#!/bin/bash

# MOPS TITANIA 3 TEST
# (C) ASTRID BOJE 2019
#
# File purpose:
#       Test hybrid particle-number/particle model (preprint 211) 

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
# Reference values from 20 runs with 1024 particles 
# SHA-1: aa72af443a1a0de64f4ddb65595ead7221add564
#####################################################################
# Run MOPS PSR
echo "Running MOPS for titania hybrid 1: PSR"
"$program" -p --strang -s "sweep-fo-spherical.xml" -r "mops-psr.xml" --ensemble > /dev/null
CheckErr $?

fname="PSR-part.csv"
m0True="4.2518E+15"
m0Err="2.1259E+14"
massTrue="8.6889E-3"	
massErr="4.3445E-4"
dpTrue="4.9767E-8"
dpErr="2.4884E-9"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
echo $m0
CheckValues $m0True $m0 $m0Err
mass=${line1[20]}
echo "Checking Mass..."
echo $mass
CheckValues $massTrue $mass $massErr
dp=${line1[8]}
echo "Checking Avg. diameter..."
echo $dp
CheckValues $dpTrue $dp $dpErr

echo "--------"

fname="PSR-total-particle-number.csv"
npTrue="222.15"
npErr="11.108"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
np=${line1[2]}
echo "Checking temperature..."
echo $np
CheckValues $npTrue $np $npErr

echo "--------"

rm -f PSR*

# If we've made it here, everything has been going well!
echo "Titaniahybrid 1 test 1/2 passed: PSR."

#####################################################################
# Run MOPS Network
echo "Running MOPS for titania hybrid 1: Network"
"$program" -p --strang -s "sweep-fo-detailed.xml" -r "mops-network.xml" -w --ensemble > /dev/null
CheckErr $?

fname="Network(stage1)-part.csv"
m0True="2.4783E+17"
m0Err="1.2392E+16"
massTrue="0.41164"	
massErr="0.020582"
dpTrue="3.9853E-8"
dpErr="1.9927E-9"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
echo $m0
CheckValues $m0True $m0 $m0Err
mass=${line1[20]}
echo "Checking Mass..."
echo $mass
CheckValues $massTrue $mass $massErr
dp=${line1[38]}
echo "Checking Avg. primary diameter..."
echo $dp
CheckValues $dpTrue $dp $dpErr

echo "--------"

fname="Network(stage2)-part.csv"
m0True="1.4543E+17"
m0Err="7.2715E+15"
massTrue="0.48191"	
massErr="0.024096"
dpTrue="7.3569E-8"
dpErr="3.6785E-9"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
echo $m0
CheckValues $m0True $m0 $m0Err
mass=${line1[20]}
echo "Checking Mass..."
echo $mass
CheckValues $massTrue $mass $massErr
dp=${line1[38]}
echo "Checking Avg. primary diameter..."
echo $dp
CheckValues $dpTrue $dp $dpErr

echo "--------"

rm -f Network*

# If we've made it here, everything has been going well!
echo "Titaniahybrid 1 test 2/2 passed: Network."

#####################################################################
# If we've made it here, everything has been going well!
echo "Titaniahybrid 1 test passed."

cd ..

exit 0
