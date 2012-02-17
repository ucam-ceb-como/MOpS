#!/usr/bin/perl

#  Copyright (C) 2010 William J Menz
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

# this test is designed for checking the silica model with WPMs

# Path of executable should be supplied as first argument to this script
#my $program = $ARGV[0];

# Parse the moments file
my $momentFileFM;
my $momentFileSF;
open($momentFileFM, "<silica-fm-part.csv") or die "ERR: failed to open moment file: $!";
open($momentFileSF, "<silica-sf-part.csv") or die "ERR: failed to open moment file: $!";

my @fm_test;
my @sf_test;

while(<$momentFileFM>) {
  my @fields = split /,/;
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 1.0) < 1e-6 )) {
      @fm_test = ($fields[4], $fields[16], $fields[8], $fields[44], $fields[46]);
      last;
  }
}

while(<$momentFileSF>) {
  my @fields = split /,/;
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 1.0) < 1e-6 )) {
      @sf_test = ($fields[4], $fields[16], $fields[8], $fields[44], $fields[46]);
      last;
  }
}


my @names = ("M0", "Fv", "dcol", "dpri", "sint level");
###################################################################
# FM comparison values, generated with N=4096, L=10 and rounded-up
# Commit  d841e368ed254bf2382fe7de5fcef5950e763a46
##################################################################
#           M0            Fv         dcol       dpri   sint level
my @fm_true = (2.85E+014, 9.66E-011, 8.31E-009, 7.537, 0.0472);
my @fm_errs = (5.20E+012, 1.0E-013,  5.0E-011,  0.03,  0.0020);
##################################################################
# SF comparison values, generated with N=512, L=10 and rounded-up
# Commit  d841e368ed254bf2382fe7de5fcef5950e763a46
##################################################################
#              M0         Fv         dcol       dpri   sint level
my @sf_true = (1.67E+014, 4.8E-09,   7.3E-08,   10.2,  0.081);
my @sf_errs = (5.19E+012, 1.0E-11,   1.3E-09,   0.05,  0.006);

my $i = 0;
for($i = 0; $i < 5; $i++) {
    if(abs($fm_test[$i] -  $fm_true[$i]) > $fm_errs[$i]) {
        print "TEST FAILED: $names[$i] expected $fm_true[$i], got $fm_test[$i]\n";
        exit ($i+1);
    } else {
        print "TEST PASSED: $names[$i] got $fm_test[$i], within $fm_true[$i] +/- $fm_errs[$i]\n";
    }
}

for($i = 0; $i < 5; $i++) {
    if(abs($sf_test[$i] -  $sf_true[$i]) > $sf_errs[$i]) {
        print "TEST FAILED: $names[$i] expected $sf_true[$i], got $sf_test[$i]\n";
        exit ($i+1);
    } else {
        print "TEST PASSED: $names[$i] got $sf_test[$i], within $sf_true[$i] +/- $sf_errs[$i]\n";
    }
}

print "All tests passed\n";
exit 0;
