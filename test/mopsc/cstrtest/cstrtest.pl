#!/usr/bin/perl

#  Copyright (C) 2011 William Menz
#
#  cstrtest
#  Tests the gas-phase chemistry for the CSTR (PSR) solver.
#  The cantera source is also included in cantera-cstr.py
#  Uses a test example for the decomposition of silane:
#     SiH4 -> Si + 2 H2
#
#  Licence:
#    This file is part of "mops".
#
#    mops is free software; you can redistribute it and/or
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
#    Dr Markus Kraft
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
open($ChemFile, "<silane-chem.csv") or die "ERR: failed to open chemistry file: $!";

my $SIH4 = 0;
my $SI = 0;
my $H2 = 0;

while(<$ChemFile>) {
  my @fields = split /,/;

  # Look for a line that begins with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.0001
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 20) < 1e-6 )) {
      $SIH4 = $fields[6];
      $SI = $fields[8];
      $H2 = $fields[4];
      last;
  }
}

print "$SIH4, $SI, $H2\n";
if((abs($SIH4 -  3.80E-09) / 3.80E-09) > 0.5) {
  print "Final SiH4 composition was $SIH4, when  3.80E-09 mol/cm3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if((abs($SI -  4.60E-09) / 4.60E-09) > 0.5) {
  print "Final Si composition was $SI, when  3.80E-09 mol/cm3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}
if((abs($H2 -  9.19E-09) / 9.19E-09) > 0.5) {
  print "Final H2 composition was $H2, when  3.80E-09 mol/cm3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

close $ChemFile;
print "All tests passed\n";

exit 0;
