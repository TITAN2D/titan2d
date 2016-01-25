#!/usr/bin/perl

my $infoDir = shift;
my $outDir = shift;
my $hemi = shift;
chomp($hemi);
open(ZONEFILE, ">$outDir"."/zone.txt");
my $grepString = `grep -ir "zone:" $infoDir`; 
if($grepString =~ m/zone:\s*(\d\d)$/){
   print ZONEFILE "$1 \n";
   print ZONEFILE "$hemi \n";
}
