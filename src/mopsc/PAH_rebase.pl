#!/usr/bin/perl

# This script is designed to stip out the first column of output from a premix calculation
# and use it as the vertices of a grid that can be specified for brush.

use strict;
use warnings;

my $fileName;
foreach $fileName (@ARGV) {
  my $dataFile;
  open($dataFile, $fileName) or die "Failed to open $fileName\n";

  # Read the first line to find when the profile begins
  my $firstLine = <$dataFile>;

  # The start time is the floating point number at the beginning of the line
  $firstLine =~ /^([-+]?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?),(.*)/;
  my $startTime = $1;
  print "0.0,$4\n";

  while(<$dataFile>){
    if(/^([-+]?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?),(.*)/) {
      my $time = $1 - $startTime;
      print "$time,$4\n";
    }
  }
}
