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
# SHA-1: 
#####################################################################
# Run MOPS PSR
echo "Running MOPS for titania 3: PSR"
"$program" -p --strang -s "sweep-fo-spherical.xml" -r "mops-psr.xml" --ensemble > /dev/null
CheckErr $?

fname="PSR-part.csv"
m0True="9.655E+16"
m0Err="3.8E+15"
massTrue="1.659E-4"	
massErr="2.4E-6"
dpTrue="7.786E-9"
dpErr="9.0E-11"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
CheckValues $m0True $m0 $m0Err
mass=${line1[20]}
echo "Checking Mass..."
CheckValues $massTrue $mass $massErr
dp=${line1[58]}
echo "Checking Avg. primary diameter..."
CheckValues $dpTrue $dp $dpErr

fname="PSR-chem.csv"
tempTrue="9.655E+16"
tempErr="3.8E+15"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
temp=${line1[60]}
echo "Checking temperature..."
CheckValues $tempTrue $temp $tempErr

rm -f PSR*

# If we've made it here, everything has been going well!
echo "Titaniahybrid 1 test 1/2 passed: PSR."

#####################################################################
# Run MOPS Network
echo "Running MOPS for titania 3: Network"
"$program" -p --strang -s "sweep-fo-detailed.xml" -r "mops-network.xml" -w --ensemble > /dev/null
CheckErr $?

fname="Network(stage1)-part.csv"
m0True="9.655E+16"
m0Err="3.8E+15"
massTrue="1.659E-4"	
massErr="2.4E-6"
dpTrue="7.786E-9"
dpErr="9.0E-11"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
CheckValues $m0True $m0 $m0Err
mass=${line1[20]}
echo "Checking Mass..."
CheckValues $massTrue $mass $massErr
dp=${line1[58]}
echo "Checking Avg. primary diameter..."
CheckValues $dpTrue $dp $dpErr

fname="Network(stage2)-part.csv"
m0True="9.655E+16"
m0Err="3.8E+15"
massTrue="1.659E-4"	
massErr="2.4E-6"
dpTrue="7.786E-9"
dpErr="9.0E-11"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
CheckValues $m0True $m0 $m0Err
mass=${line1[20]}
echo "Checking Mass..."
CheckValues $massTrue $mass $massErr
dp=${line1[58]}
echo "Checking Avg. primary diameter..."
CheckValues $dpTrue $dp $dpErr

rm -f Network*

# If we've made it here, everything has been going well!
echo "Titaniahybrid 1 test 2/2 passed: Network."

#####################################################################
# If we've made it here, everything has been going well!
echo "Titaniahybrid 1 test passed."

cd ..

exit 0