#!/usr/bin/perl

#  Copyright (C) 2011 Laurence R. McGlashan.
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
my $profileFile;
open($profileFile, "<profile.dat") or die "ERR: failed to open profile.dat file: $!";

# Parse the moments file
my $originalProfileFile;
open($originalProfileFile, "<profileOriginal.dat") or die "ERR: failed to open profileOriginal.dat file: $!";

my @timeNew = ();
my @timeOld = ();
my @tempNew = ();
my @tempOld = ();

while(<$profileFile>) {
    my @fields = split(" +");
    push(@timeNew, $fields[0]);
    push(@tempNew, $fields[3]);
}
while(<$originalProfileFile>) {
    my @fields = split(" +");
    push(@timeOld, $fields[0]);
    push(@tempOld, $fields[3]);
}

#Compare temperature in old and new output files, skip the first line (line 0)
# because it contains the column headings.
for(my $i=1;$i<@tempOld;++$i) {
    if (abs($tempOld[$i] - $tempNew[$i])/$tempOld[$i] > 1e-4) {
      print "Old temp was $tempOld[$i], but new temp is $tempNew[$i] ($i)\n";
      print "**************************\n";
      print "****** TEST FAILURE ******\n";
      print "**************************\n";
      exit 1;
    }
}

#Compare time in old and new output files, skip the title line and the first line
#of data which should be zero.
for(my $i=2;$i<@timeOld;++$i) {
    if (abs($timeOld[$i] - $timeNew[$i])/$timeOld[$i] > 1e-4) {
      print "Old time was $timeOld[$i], but new time is $timeNew[$i] ($i)\n";
      print "**************************\n";
      print "****** TEST FAILURE ******\n";
      print "**************************\n";
      exit 1;
    }
}
exit 0;
