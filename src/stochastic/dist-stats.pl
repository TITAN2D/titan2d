#!/usr/bin/perl -w
#
# A simple load balancing Perl script for serial jobs.  Basically
#  works as a daemon for farming out NJOBS jobs to NCPUS CPUs,
#  where (hopefully) NJOBS >> NCPUS.
#  Should print timing statistics at end of the script (provided
#  it finishes cleanly).
#
# Usage:
#        Modify the first section to create job lists (note that
#         you will need to modify the task launching commands as
#         well if they are job dependent (e.g. require unique input
#         filenames).
#
# N.B.
#        deletion of array (non-hash) elements requires Perl >= 5.6
#          (so use hashes for dynamic arrays)
#
# Written: Winter 2004 (mdj) (modified by krd)
# Modified: 2005-01 (mdj) stats version
#
# Requires one argument - name of the input file containing the run
#  parameters (needs to be one of "stat_ctl.[bed,vol,loc]")
#
# Uses rsh to launch remote jobs
# Requires PBS (relaxed 2005-02-01, can use $TMPDIR rather than $PBSTMPDIR,
#                $MACHLIST rather than $PBS_NODEFILE)
#
use strict;
use Benchmark;
use POSIX ":sys_wait_h";
use Getopt::Long;
#
my $USAGE = "\nUsage: $0 --ctlfile stat_ctl.[bed,vol,loc] [--help]\n\n" .
  "    A script for farming out TITAN jobs to a newtork of workstations\n" .
  "      for statistical analysis.\n" .
  "      Requires:  PBS (uses \$PBSTMPDIR for temporary directories for\n" .
  "        running jobs and \$PBS_NODEFILE for a machine list of available\n".
  "        CPUs).  If not using PBS, the env variable \$TMPDIR and\n" .
  "        \$MACHFILE should be set.\$MACHFILE should point to a plain \n" .
  "        text file that contains a list of machines to which you have \n" .
  "        unrestricted (i.e. password-less) rsh access.  If you choose \n" .
  "        not to use PBS and not to set \$MACHLIST and \$TMPDIR, the \n" .
  "        default behavior is to run as many tasks on your local machine \n" .
  "        as you have CPUs, and use a directory tmp.SAMPLE as the scratch\n" .
  "        space.\n " .
  "      Current limitations:\n" .
  "       1. LHS (Latin Hypercube Sampling) is assumed.\n" .
  "       2. Each independent run will be on a single processor\n" .
  "          (i.e., non-parallel).\n" .
  "       3. Restarts must be handled by hand - you must edit the stat_ctl\n" .
  "          file and insert a minus sign before the unique id (2nd column)\n".
  "          of the jobs you do not wish to be recalculated.\n";
#
# debug Benchmark 1;
#
my $POLL_INTERVAL = 2; # polling interval to check for job completion
#
# Job description section - customize to suit application
#
my $script_start = new Benchmark;
#
#  ctlfile = master control file
#  help = print help output (default)
#
my ($ctlfile,$help);
my $ok = GetOptions("ctlfile=s" => \$ctlfile,
		    "help" => \$help);
die "$USAGE" unless $ok;
if ($ctlfile) {
  open(IN,"<$ctlfile") or die "Unable to open control file $ctlfile\n";
} else {
  die "$USAGE";
}
if ($help) {
  die "$USAGE";
}
# save the control file suffix - it specifies the type of statistics run
my $stats_type = $ctlfile;
$stats_type =~ s/stat_ctl\.//;
print "Stats run, type = ",$stats_type,"\n";
#
# Process control file; negative unique job IDs (2nd column) are ignored
#
# JOB_type = [lhs,pcq]
# JOB_level = refinement level
# JOB_inputs = actual variable entries to be passed to TITAN
#
my $JOB_type="lhs"; # only supported option at this point
my (%JOB_inputs,%JOB_level);
my ($ijob,$NJOBS) = (-1,0);
while (<IN>) {
  # grab the line, parse/store it it for job info
  chomp;
  my $ctl_line = $_;
  #
  # 1st column is refinement level (LHS) or weight (PCQ)
  # 2nd column is unique ID - if <0 we toss this one
  # 3rd column+ is job run info for TITAN - _entire line_ must be stashed
  #
  my @words = split " ",$_;
  print "control file line = ",$_,"\n";
  my $level = $words[0];
  my $uniqueID = $words[1];
  print "level = ",$level," uniqueID = ",$uniqueID,"\n";
  if ($uniqueID >= 0) {
    $ijob++;
    $JOB_level{$ijob} = sprintf("%02d",$level);
    $JOB_inputs{$ijob} = $ctl_line;
    #$JOB_inputs{$ijob} = $words[2];
    #for (my $i=3; $i<=scalar(@words)-1; $i++) {
    #  $JOB_inputs{$ijob} = $JOB_inputs{$ijob} . " " . $words[$i];
    #}
  }
}
close (IN);
$NJOBS = $ijob+1;
printf "Number of Jobs = %d\n",$NJOBS;
for (my $i=0; $i<=$NJOBS-1; $i++) {
  printf "Job entry %d, %s, %s\n",$i,$JOB_level{$i},$JOB_inputs{$i};
}
#my $NJOBS = $maxrunno - $minrunno + 1;
my $NCPUS_per_JOB = 1; # Assumption
#
# Creates hash of actual jobs (all the jobs)
#
# JOB_list{$ijob} - actual run invocation, including setup and post-processing
#                   file movement
# JOB_inouts{$ijob} - input file overrides for TITAN, from the
#                   stat_ctl.[bed,loc,vol] control file
# local_run - indicates a non-PBS run on only the local machine (default if
#                   no $MACHFILE present)
# ncpus_local - tries to get the number of cpus from /proc filesystem
#
my ($ncpus_local,$local_run) = (0,0);
my %JOB_list;
my ($WORKDIR,$TMPDIR,$NODEFILE);
if ($ENV{PBS_ENVIRONMENT}) {
  $WORKDIR = $ENV{PBS_O_WORKDIR};
  $TMPDIR = $ENV{PBSTMPDIR};
  $NODEFILE = $ENV{PBS_NODEFILE};
} else { # not PBS
  $WORKDIR = `pwd`;
  chomp $WORKDIR;
  $TMPDIR = $ENV{TMPDIR};
  if (! $TMPDIR) { #set default TMPDIR
    $TMPDIR = `pwd`;
    chomp $TMPDIR;
    $TMPDIR = $TMPDIR . "/tmp.SAMPLE";
    system("mkdir $TMPDIR");
    print "TMPDIR = $TMPDIR","\n";
  }
  $NODEFILE = $ENV{MACHFILE};
  if (! $NODEFILE) { #
    $NODEFILE = "tmp.machfile";
    printf "NODEFILE = %s\n",$NODEFILE;
    # try to get number of cpus from /proc filesystem
    open(CPUINFOFILE,"</proc/cpuinfo");
    open(TMPNODEFILE,">$NODEFILE");
    while (<CPUINFOFILE>) {
      if ($_ =~ /processor/) { $ncpus_local++; }
    }
    close(CPUINFOFILE);
    printf "Number of Local CPUs (this machine): %d\n",$ncpus_local;
    for (my $i=1; $i<=$ncpus_local; $i++) {
        print TMPNODEFILE `hostname`;
    }
    close(TMPNODEFILE);
    $local_run = 1;
  }
}
if (!-e $NODEFILE) {
  die "Unable to find \$PBS_NODEFILE or \$MACHFILE; check your environment.\n";
}
if (!-e $TMPDIR) {
  die "Unable to find TMPDIR; you must have either \$TMPDIR or \$PBSTMPDIR\n" .
    "   set in your environment\n";
}
open(IN,"<$NODEFILE");
my $icpu = -1;
my @CPU_list;
#my %CPU_mpiexec_conf;
while (<IN>) {
  $icpu ++;
  chomp;
  $CPU_list[$icpu] = $_;
#  $CPU_mpiexec_conf{$icpu} = "rm -f conf ; echo \"$CPU_list[$icpu] : titan\" > conf ;"
  printf "CPU_list entry %d = %s\n",$icpu,$CPU_list[$icpu];
}
my $NCPUS = scalar(@CPU_list);
printf "Number of CPUs on which to spawn jobs: %d\n",$NCPUS;
close(IN);

my $sampleno;
my $runno =0;
my $dest;
my $orig;
for (my $ijob=0; $ijob <= $NJOBS-1; $ijob++) {
  #$runno = $ijob+$minrunno;
  $runno = $ijob;
  $sampleno = sprintf("%06d",$runno);
  #$dest =$TMPDIR/RUN$runno;
  #$orig =$WORKDIR/RUN$runno;
  $JOB_list{$ijob} = "mkdir $TMPDIR/SAMPLE.$ijob ; cd $TMPDIR/SAMPLE.$ijob ;" .
    "cd $WORKDIR ; cp -f funky0000.inp simulation.data titan " .
    " $TMPDIR/SAMPLE.$ijob ; cd $TMPDIR/SAMPLE.$ijob ; " .
      " echo $JOB_inputs{$ijob} > statin.$stats_type ;  ./titan >&/dev/null;" .
	"cat statout_lhs.$JOB_level{$ijob} >> $WORKDIR/statout_lhs.$JOB_level{$ijob} ";
#  $JOB_list{$ijob} = "sleep 5"; # dummy job for testing purposes
}
for (my $i=0; $i<=$NJOBS-1; $i++) {
  printf "Job entry %d, %s\n",$i,$JOB_list{$i};
}

#
# Available CPUs section
#
if ($NCPUS > $NJOBS) { die "You know, it is not nice to leave idle CPUs (NCPUS>NJOBS)"; }
#
# 1st set of Jobs Launching Section
#
my $repeats_over_CPUs = int( $NJOBS/$NCPUS );
printf "Number of repeats over CPUs: %d\n",$repeats_over_CPUs;
$ijob=-1;
my (%wall_start,%wall_end,%wall_elapsed);        # arrays for walltime stats
my (%pids_running,%cpus_running,%jobs_running);
for (my $i=$NCPUS-1; $i>=0; $i--) {
  $ijob ++;
  #
  # Construct command string:
  #  rsh - remote machines
  #  sh  - local spawns (note -c option to force command line
  #        processing; how portable is this flag?)
  #
  my $command; # actual command string for each job
  if ($local_run) { # just background right number of tasks
    #$command = " '" . $JOB_list{$ijob} . " & '";
    $command = "sh -c " . " '" . $JOB_list{$ijob} . " '";
  } else { # rsh required
    $command = "rsh " . " " . $CPU_list[$i] . " '" . $JOB_list{$ijob} . " '";
  }
  #
  # fork tasks
  #
  my $pid_rsh = fork;
  if ($pid_rsh == 0) { # I am a kid
    exec "$command" or die "Unable to rsh $CPU_list[$i]";
  } else { # I am the master
    $wall_start{$ijob} = new Benchmark;
    $pids_running{$i} = $pid_rsh;
    $cpus_running{$i} = $CPU_list[$i];
    $jobs_running{$i} = $ijob;
    my $start_time = `date -Iseconds`;
    chomp $start_time;
    printf "Job %4d started   on CPU %s at %s, Command=%s\n",$ijob,
      $CPU_list[$i],$start_time,$command;
  }
}
#
#
my @left = keys %JOB_list;
my $todo = scalar(@left);
printf "Number of jobs to do   = %d\n",$todo;
printf "Number of jobs spawned = %d\n",scalar(keys(%pids_running));
#
#while (defined($todo)) {
while ($todo > 0) {
  sleep $POLL_INTERVAL;
  foreach my $key (keys %jobs_running) {
    my $kid = waitpid($pids_running{$key},&WNOHANG);
    if ($kid == -1) { # reap completed job, launch new one
      my $end_time = `date -Iseconds`;
      chomp $end_time;
      $wall_end{$jobs_running{$key}} = new Benchmark;
      printf "Job %4d completed on CPU %s at %s\n", $jobs_running{$key},
	$cpus_running{$key},$end_time;
      #delete $jobids_running{$key};
      #delete $cpus_running{$key};
      delete $JOB_list{$jobs_running{$key}};
      #delete $JOB_list2{$jobs_running{$key}};
      if ($ijob < $NJOBS-1) { # still work to be allocated
        $ijob ++;
	#my $command = "rsh " . " " . $cpus_running{$key} . " '" . 
	#  $JOB_list{$ijob} . " '";
	my $command;
	if ($local_run) { # just background right number of tasks
	  #$command = " '" . $JOB_list{$ijob} . " & '";
	  $command = "sh -c " . " '" . $JOB_list{$ijob} . " '";
	} else { # rsh required
	  $command = "rsh " . " " . $cpus_running{$key} . " '" . 
	    $JOB_list{$ijob} . " '";
	}
        my $pid_rsh = fork;
        if ($pid_rsh == 0) { # I am a kid
          exec "$command" or die "Unable to execute command: ",$command;
        } else { # I am the master - update lists with new info
          $wall_start{$ijob} = new Benchmark;
          $pids_running{$key} = $pid_rsh;
          # $cpus_running{$key} = $CPU_list[$i];
          $jobs_running{$key} = $ijob;
          my $start_time = `date -Iseconds`;
          chomp $start_time;
          printf "Job %4d started   on CPU %s at %s, Command=%s\n",$ijob,
	    $CPU_list[$key],
	    $start_time,$command;
        }
      } else { # no remaining work to farm out, just clean up
        delete $jobs_running{$key};
        delete $cpus_running{$key};
        delete $pids_running{$key};
      }
    } # kid == -1
  } # key loop over running jobs
  @left = keys %JOB_list;
  $todo = scalar(@left);
  #printf "Loop: Number of jobs running = %d, in pool = %d\n",scalar(keys %jobs_running),$todo;
} # defined($todo)

#
# Final stats
my $sum_walltimes = 0;
foreach my $key (keys %wall_end) {
  if (!defined($wall_start{$key})) {
    die "Could not find start time for Job ",$key;
  }
#  print "Job = ",$key, "start = ",timestr($wall_start{$key})," end = ",timestr($wall_end{$key}),"\n";
  $wall_elapsed{$key} = timediff($wall_end{$key},$wall_start{$key});
  print "Job = ",$key," elapsed time = ",timestr($wall_elapsed{$key}),"\n";
  my @timefields = split(" ",timestr($wall_elapsed{$key}));
  my $wtime=$timefields[0];
#  print "Wtime = ",$wtime,"\n";
  $sum_walltimes += $wtime;
}
my $script_end = new Benchmark;
my $script_elapsed = timediff($script_end,$script_start);
print "Runtime for Job Balancing Script = ",timestr($script_elapsed),"\n";
my @timefields = split(" ",timestr($script_elapsed));
my $wtime_script = $timefields[0];
print "Sum of Execution walltimes = ",$sum_walltimes,"\n";
printf "Scaling factor: %9.4f = %5.2f%%\n",$sum_walltimes/$wtime_script,
  100.0*($sum_walltimes/$wtime_script)/$NCPUS;
