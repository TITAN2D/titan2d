#!/usr/bin/perl -w
#
#
# Originally Written: 2005-01 (mdj)
#
use strict;
use Getopt::Long; # not currently being used (for later)
use POSIX;        # ceil and other math functions
#
#
my $stats_file = "statout.plot";
my $NFIELDS = 20; # number of expected fields in stats output file
my @FIELDFILES = ("nsamples","prob","vave_mean","vave_stddev","vave_skew",
  "hmax_mean","have_stddev","have_skew","xc_mean","xc_stddev","xc_skew",
  "yc_mean","yc_stddev","yc_skew","xspr_mean","xspr_stddev","xspr_skew",
  "yspr_mean","yspr_stddev","yspr_skew");
my @FIELDLABELS = ("N_{samples}","Prob","\$v_{ave}:mean\$","vave_stddev","vave_skew",
  "hmax_mean","have_stddev","have_skew","xc_mean","xc_stddev","xc_skew",
  "yc_mean","yc_stddev","yc_skew","xspr_mean","xspr_stddev","xspr_skew",
  "yspr_mean","yspr_stddev","yspr_skew");

my $NFIELDFILES = scalar(@FIELDFILES);

open(STATSDATA,"<$stats_file") or die "Unable to open stats input file\n";
#
# parse the consolidated data file into a form better suited for gnuplot
#  (why can't gnuplot handle multiple columns?)
#
my $stats_file_line = -1;
my (@nsamples,@prob,@vave_mean,@vave_stddev,@vave_skew,@hmax_mean,@hmax_stddev,@hmax_skew,
    @xc_mean,@xc_stddev,@xc_skew,@yc_mean,@yc_stddev,@yc_skew,@xspr_mean,@xspr_stddev,@xspr_skew,
    @yspr_mean,@yspr_stddev,@yspr_skew);
my @outputdata;

while(<STATSDATA>) {
  chomp; #get rid of EOL
  my @words = split " ",$_;
  my $nwords = scalar(@words);
  #
  #
  # fixed format of stats file
  #  c1 - number of samples (ordinate)
  #  c2 - probability
  #  c3-5   vave
  #  c6-8   hmax
  #  c9-11  xc
  #  c12-14 yc
  #  c15-17 xspread
  #  c18-20 yspread
  #
  #
  if ($nwords != $NFIELDS) {
    die "Found only ",$nwords," columns in ",$stats_file,"\n";
  }
  $stats_file_line ++;
  $nsamples[$stats_file_line] = $words[0];
  $prob[$stats_file_line] = $words[1];
  $vave_mean[$stats_file_line] = $words[2];
  $vave_stddev[$stats_file_line] = $words[3];
  $vave_skew[$stats_file_line] = $words[4];
  $hmax_mean[$stats_file_line] = $words[5];
  $hmax_stddev[$stats_file_line] = $words[6];
  $hmax_skew[$stats_file_line] = $words[7];
  $xc_mean[$stats_file_line] = $words[8];
  $xc_stddev[$stats_file_line] = $words[9];
  $xc_skew[$stats_file_line] = $words[10];
  $yc_mean[$stats_file_line] = $words[11];
  $yc_stddev[$stats_file_line] = $words[12];
  $yc_skew[$stats_file_line] = $words[13];
  $xspr_mean[$stats_file_line] = $words[14];
  $xspr_stddev[$stats_file_line] = $words[15];
  $xspr_skew[$stats_file_line] = $words[16];
  $yspr_mean[$stats_file_line] = $words[17];
  $yspr_stddev[$stats_file_line] = $words[18];
  $yspr_skew[$stats_file_line] = $words[19];
  for (my $j=0;$j<$NFIELDS;$j++) {
    $outputdata[$stats_file_line][$j] = $words[$j];
  }
  print "Opened stats file: ",$stats_file," found ",$nwords," columns on line ",
    $stats_file_line,".\n";
  print " nsamples = ",$nsamples[$stats_file_line]," prob = ",$outputdata[$stats_file_line][1],"\n";
}
close(STATSDATA);
my $total_nsamples = scalar(@nsamples);
#
# Output temporary gnuplot files
#
my @STATSOUTFILE;
for (my $i=1; $i<$NFIELDS;$i++) {
  $STATSOUTFILE[$i] = "tmp_" . $FIELDFILES[$i];
  open(OUTFILE,">$STATSOUTFILE[$i]") or die "Unable to open output file: ",$STATSOUTFILE[$i],"\n";
  for (my $j=0; $j<$total_nsamples; $j++) {
    printf OUTFILE "%d %s\n",$nsamples[$j],$outputdata[$j][$i];
    #printf "%d %s\n",$nsamples[$j],$outputdata[$j][$i];
  }
  close(OUTFILE);
}
#
# Construct gnuplot script file
#
my $max_screen_x = 1.0;
my $max_screen_y = 1.0;
my $nx = 2;
my $ny = 4;
my $plots_per_page = $nx*$ny;
my $num_pages = ceil( ($NFIELDS-1)/$plots_per_page );
print "Number of output pages: ",$num_pages," plots per page: ",$plots_per_page,"\n";
my $dx = $max_screen_x/$nx;
my $dy = $max_screen_y/$ny;

#my ($GNUFILENM,$PSFILENM,$pos_x,$pos_y);
#my ($GNUFILENM,$PSFILENM);
my $iplot = 0;
for (my $ipage=1; $ipage<=$num_pages; $ipage++) {
  #
  # construct gnuplot file name for this page
  #
  #my $string_ipage = sprintf "%s",$ipage;
  my $GNUFILENM = "plot_stats_p" . $ipage . ".gnuplot";
  #$GNUFILENM = "plot_stats.gnuplot";
  open(GNUFILE,">$GNUFILENM") or die "Unable to open gnuplot file: ",$GNUFILENM,"\n";
  my $PSFILENM = "plot_stats_p" . $ipage . ".ps";
  #$PSFILENM = "plot_stats_1.ps";
  printf GNUFILE "set terminal postscript portrait\nset output \"%s\"\n",$PSFILENM;
  printf GNUFILE "set size %f, %f\nset origin 0.0, 0.0\n",$max_screen_x,$max_screen_y;
  printf GNUFILE "set multiplot\n";
  printf GNUFILE "set key left\n";
  #printf GNUFILE "set size %f, %f\n",$dx,$dy;
  #
  # Fill plots from top, left->right
  #
  my $pos_x = 0.0;
  my $pos_y = $max_screen_y-$dy;
  for (my $i=1; $i<=$plots_per_page;$i++) {
    $iplot++;
    if ($iplot < $NFIELDS) {
      printf GNUFILE "set origin %f, %f\n",$pos_x,$pos_y;
      printf GNUFILE "set size %f, %f\n",$dx,$dy;
      printf GNUFILE "set logscale x\n";
      printf GNUFILE "plot \"%s\" title \"%s\" with linespoints\n",$STATSOUTFILE[$iplot],
	$FIELDLABELS[$iplot];
    }

    $pos_x += $dx;
    if ($pos_x >= ($max_screen_x-0.1*$dx) ) {
      # decrement y, go again
      $pos_y -= $dy;
      $pos_x = 0.0;
    }
    if ($pos_y < 0.0) { # next page
    }
  }
  #printf GNUFILE "unset multiplot\n";
  close(GNUFILE);
  system("/usr/bin/gnuplot $GNUFILENM");
} # ipage
#
# remove temporary files
#
unlink(glob "tmp_*");
#system("/Projects/CCR/jonesm/gnuplot/3.8k.2/bin/gnuplot $GNUFILENM");
