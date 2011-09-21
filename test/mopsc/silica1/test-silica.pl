#!/usr/bin/perl

#  Copyright (C) 2009 Robert I A Patterson.
#	modified by wjm34, 2011
#
# Licence:
#    This file is part of "mops".
#
#    brush is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 2
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#  Contact:
#    Prof Markus Kraft
#    Dept of Chemical Engineering
#    University of Cambridge
#    New Museums Site
#    Pembroke Street
#    Cambridge
#    CB2 3RA
#    UK
#
#    Email:       mk306@cam.ac.uk
#    Website:     http://como.cheng.cam.ac.uk

use strict;
use warnings;

# this test is designed for checking whether PAH-PP model work correctly by comparing to precalculated value

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Parse the moments file
my $momentFile;
open($momentFile, "<silica-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $dcol = 0;
my $dpri = 0;
my $sl = 0;
my $num_si = 0;
my $num_o = 0;
my $si_to_o = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begina with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.8
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.8) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], \n";

      $dcol = $fields[8];
      #print "8: $fields[8] \n";
      
      $dpri = $fields[44];
      #print "44: $fields[44] \n";
      
      $sl = $fields[46];
      #print "46: $fields[46] \n";
      
      $num_si = $fields[36];
      
      $num_o = $fields[38];
      last;
  }
}

$dcol = 1.0e9 * $dcol;
$si_to_o = $num_si / $num_o;
#print "si/o: $si_to_o\n";

####################################################
# Comparison with test values (Boost 1.47)
####################################################
# 99.9% CI shown for 10 runs
####################################################
# m0 = 2.10e14 +/- 3.3e12
# dcol = 78.8 +/- 1.26 nm
# dpri = 7.90 +/- 0.3 nm
# sl = 0.079 +/- 0.010
# si_to_o = 0.549 =/- 0.020
####################################################

if(abs($m0 -  2.10e14) > 3.3e12) {
  print "Simulated mean M0 was $m0, when 2.10e14m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}
if(abs($dcol -  78.8) > 1.26) {
  print "Simulated mean dcol was $dcol, when 78.8 nm expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}
if(abs($dpri -  7.90) > 0.30) {
  print "Simulated mean dpri was $dpri, when 7.90 nm expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}
if(abs($sl -  0.079) > 0.010) {
  print "Simulated mean sl was $sl, when 0.079 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 4;
}
if(abs($si_to_o -  0.549) > 0.020) {
  print "Simulated mean si_to_o was $si_to_o, when 0.549 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 5;
}

print "All tests passed\n";
exit 0;
