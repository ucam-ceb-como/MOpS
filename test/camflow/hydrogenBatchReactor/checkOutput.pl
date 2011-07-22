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
    push(@timeNew, "Time: $fields[0] \n");
    push(@tempNew, "Time: $fields[3] \n");
}
while(<$originalProfileFile>) {
    my @fields = split(" +");
    push(@timeOld, "Time: $fields[0] \n");
    push(@tempOld, "Time: $fields[3] \n");
}

#Compare temperature in old and new output files.
for(my $i=0;$i<@tempOld;++$i) {
    if ($tempOld[$i] ne $tempNew[$i]) {
      print "**************************\n";
      print "****** TEST FAILURE ******\n";
      print "**************************\n";
      exit 1;
    }
}

#Compare time in old and new output files.
for(my $i=0;$i<@timeOld;++$i) {
    if ($timeOld[$i] ne $timeNew[$i]) {
      print "**************************\n";
      print "****** TEST FAILURE ******\n";
      print "**************************\n";
      exit 1;
    }
}
exit 0;