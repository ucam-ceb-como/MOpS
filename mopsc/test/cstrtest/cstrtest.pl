#!/usr/bin/perl

#  Copyright (C) 2010 Rebecca C Riehl.
#
#
# Licence:
#    This file is part of "brush".
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

# Parse the moments file
my $ChemFile;
open($ChemFile, "<J3-chem.csv") or die "ERR: failed to open chemistry file: $!";


my $TiCl4 = 0;
my $TiO2 = 0;

while(<$ChemFile>) {
  my @fields = split /,/;

  # Look for a line that begins with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.0001
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.05) < 1e-6 )) {
      # Second field should be TiCl4 composition
      $TiCl4 = $fields[2];
      #print "2: $fields[2], ";

      $TiO2 = $fields[58];
      #print "58: $fields[58] \n";

      last;
  }
}


my $SensFile;
open($SensFile, "<J3-sensi.csv") or die "ERR: failed to open sensitivity file: $!";

my $TiCl4Sens = 2;

while(<$SensFile>) {
  my @fields2 = split /,/;

  # Look for a line that begins with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.0001
  if(($fields2[0] =~ /^\d+/) && (abs($fields2[0] - 0.05) < 1e-6 )) {
      # Second field should be TiCl4 composition
      $TiCl4Sens = $fields2[3];
      #print "2: $fields[2], ";


      last;
  }
}

print "$TiCl4, $TiO2, $TiCl4Sens\n";
if(abs($TiCl4 -  8.408e-7) > 1e-8) {
  print "Final TiCl4 composition was $TiCl4, when  8.408e-7 mol/cm3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($TiO2 - 3.265e-6) > 1e-7) {
  print "Final TiO2 composition was $TiO2, when 3.265e-6 mol/cm3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
  }
  
  if(abs($TiCl4Sens - 0.0223) > 1e-3) {
  print "Final TiCl4 Sensitivity was $TiCl4Sens, when 0.0223 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}

close $ChemFile;
close $SensFile;
print "All tests passed\n";

exit 0;