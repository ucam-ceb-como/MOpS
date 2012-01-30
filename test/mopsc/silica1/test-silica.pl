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

# this test is designed for checking the silica model under direct simulation

# Path of executable should be supplied as first argument to this script
#my $program = $ARGV[0];

# Parse the moments file
my $momentFile;
open($momentFile, "<silica-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $fv = 0;
my $dcol = 0;
my $dpri = 0;
my $sl = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begina with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.8
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.8) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], \n";
      
      $fv = $fields[16];

      $dcol = $fields[8];
      #print "8: $fields[8] \n";
      
      $dpri = $fields[44];
      #print "44: $fields[44] \n";
      
      $sl = $fields[46];
      #print "46: $fields[46] \n";
      
      last;
  }
}

print "$m0 $fv $dcol $dpri $sl\n";

####################################################
# Comparison with test values (Boost 1.47)
####################################################
# 99.9% CI shown for Nsp = 512, L = 10
# actual test uses Nsp = 512, L = 4 fo efficiency
# Commit  6ed6f4b5f4b3670ebe25f2344bb86c5900197b94
####################################################
# m0: 2.13E+014	7.00E+012
# fv: 7.22E-009	1.38E-010
# dc: 8.56E-008	2.74E-009
# dp: 8.71433	0.464715
# sl: 0.0601418	0.0144396
####################################################

if(abs($m0 -  2.13e14) > 7.00E+012) {
  print "Simulated mean M0 was $m0, when 2.26e14 1/m3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}
if(abs($fv -  7.22E-009) > 1.38E-010) {
  print "Simulated mean Fv was $fv, when 7.22e-09 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 5;
}
if(abs($dcol -  8.56E-008) > 2.74E-009) {
  print "Simulated mean dcol was $dcol, when 7.92e-8 m expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}
if(abs($dpri -  8.71433) > 0.464715) {
  print "Simulated mean dpri was $dpri, when 8.47 nm expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}
if(abs($sl -  0.0601418) > 0.0144396) {
  print "Simulated mean sl was $sl, when 0.079 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 4;
}

print "All tests passed\n";
exit 0;
