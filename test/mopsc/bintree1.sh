#!/bin/bash

# MOPS BINTREE 1 TEST
# (C) WILLIAM MENZ (wjm34) 2012
#
# File purpose:
#       This test is designed to check the generic particle model,
#       BintreePrimary. It is also demonstrated how a model similar
#       to SilicaPrimary can be implemented in the generic format.
#       Full particle ensembles are also written to verify the
#       BintreeSerializer class.

#Path to executable should be supplied as first argument to
#this script.  Script will fail and return a non-zero value
#if no executable specified.
program=$1

if test -z "$program"
  then
    echo "No executable supplied to $0"
    exit 255
fi

# An optional second argument may specify the working directory
if test -n "$2"
  then
    cd "$2"
    echo "changed directory to $2"
fi
cd bintree1

# Declare check function
function CheckErr {
    if [ "$1" -gt "0" ]; then
        echo "#############################"
        echo "FAILED POSTPROCESSING."
        echo "ERROR CODE $1"
        echo "#############################"
        exit $1
    fi
}

# Run the executable
"$program" -p --strang --ensemble  > /dev/null
CheckErr $?

# Run the perl script as a here document (laziness)
perl << EOF
##################################################################
# Comparison values, generated with N=4096, L=4 and rounded-up
# Commit 3a041cbb227c543d76ddd911c75e07e677bbe31c
##################################################################
my @names = (  "M0",        "Fv",     "dcol",   "dpri",    "sint level");
my @true = (4.3391E+015, 1.7752E-09, 8.6742E-09, 6.9869E-09,  0.84218);
my @errs = (0.1526E+015, 0.0022E-09, 0.0787E-09, 0.0821E-09,  0.0050);

# Open file and parse results
my @results;
my \$ifile;
open(\$ifile, "<bintree-part.csv") or die "Can't open: \$!";

while(<\$ifile>) {
  my @fields = split /,/;
  if((@fields[0] =~ /^\d+/) && (abs(@fields[1] - 0.1) < 1e-6 )) {
      @results = (@fields[4], @fields[16], @fields[8], @fields[46], @fields[48]);
      last;
  }
}

# Loop over array to check if it matches the solution above
my \$i;
for(\$i = 0; \$i < 5; \$i++) {
    if(abs(@results[\$i] - @true[\$i]) > @errs[\$i]) {
        print "TEST FAILED: @names[\$i] expected @true[\$i], got @results[\$i]\n";
        exit (\$i+1);
    } else {
        print "TEST PASSED: @names[\$i] got \$results[\$i], within @true[\$i] +/- @errs[\$i]\n";
    }
}

exit 0
EOF
CheckErr $?

# If we've made it to here, all tests have passed!
echo "All tests passed. :D"
rm -f bintree*
