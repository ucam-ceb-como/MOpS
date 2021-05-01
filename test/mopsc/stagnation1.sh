#!/bin/bash

# MOPS STAGNATION 1 TEST
# (C) CASPER LINDBERG 2019
#
# File purpose:
#       Test flamepp with thermophoretic correction for a stagnation flame

#Path to executable should be supplied as first argument to
#this script.  Script will fail and return a non-zero value
#if no executable specified.

# First argument is executable
if [ -z "$1" ]; then
    echo "No executable supplied to $0"
    exit 255
else
    program=$1
fi

# Second argument is working directory (required)
if [ -z "$2" ]; then
    echo "No working dir supplied to $0"
    exit 254
else
    wkdir="$2"
fi
cd "$wkdir"

######################################################################
# DEFINE FUNCTIONS

# Checks for failed calculations
function CheckErr {
    if [ "$1" -gt "0" ]; then
        echo "FAILED CALCULATION WITH ERR $1."
        exit $1
    fi
}

# Uses perl to compare two floats.
function CheckValues {
perl << EOF
my \$oldval=$1;
my \$newval=$2;
my \$abserr=$3;
my \$err=0;
my \$retval=0;
if(\$oldval != 0) {
    \$err=abs(\$oldval-\$newval);
    if (\$err > \$abserr) {
            print "TEST FAILURE!!! GOT \$newval, wanted \$oldval (err \$err).\n";
            \$retval=2
    }
} 
exit \$retval
EOF

CheckErr $?
}

#####################################################################
# Reference values from 8 runs with 1024 particles 
# SHA-1: a52d39424e90f6e532d174db1a1c0cf4dab14c00
#####################################################################
fname="Z1-part.csv"

m0True="1.294E+17"
m0err="3.0E+15"
massTrue="3.008E-4"	
masserr="3.2E-6"
mass2True="1.436E-24"
mass2err="8.1E-26"

# Run MOPS
echo "Running MOPS for stagnation flame 1"
"$program" -p --flamepp -g "gasphase.csv" > /dev/null
CheckErr $?

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
CheckValues $m0True $m0 $m0err
mass=${line1[20]}
echo "Checking Mass..."
CheckValues $massTrue $mass $masserr
mass2=${line1[24]}
echo "Checking Mass2..."
CheckValues $mass2True $mass2 $mass2err

rm -f Z1*

# If we've made it here, everything has been going well!
echo "Stagnation flame 1 test passed."

cd ..

exit 0