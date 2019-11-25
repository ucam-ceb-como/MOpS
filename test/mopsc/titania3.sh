#!/bin/bash

# MOPS TITANIA 3 TEST
# (C) ASTRID BOJE 2019
#
# File purpose:
#       Test energy balance with particle contributions for titania (preprint 242)

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
"$program" -p --strang -s "sweep-fo-detailed.xml" -r "mops-psr.xml" --ensemble > /dev/null
CheckErr $?

fname="PSR-part.csv"
m0True="1.7738E+17"
m0Err="8.8690E+15"
massTrue="0.019853"	
massErr="9.9265E-4"
dpTrue="2.6238E-8"
dpErr="1.3119E-9"

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

echo "---------"

fname="PSR-chem.csv"
tempTrue="1230.2"
tempErr="61.510"

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
temp=${line1[60]}
echo "Checking temperature..."
echo $temp
CheckValues $tempTrue $temp $tempErr

echo "---------"

rm -f PSR*

# If we've made it here, everything has been going well!
echo "Titania 3 test 1/2 passed: PSR."

#####################################################################
# Run MOPS Network
echo "Running MOPS for titania 3: Network"
"$program" -p --strang -s "sweep-fo-detailed.xml" -r "mops-network.xml" -w --ensemble > /dev/null
CheckErr $?

fname="Network(stage1)-part.csv"
m0True="7.1754E+17"
m0Err="3.5877E+16"
massTrue="0.37676"	
massErr="1.8838E-2"
dpTrue="3.1998E-8"
dpErr="1.5999E-9"

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

echo "---------"

fname="Network(stage2)-part.csv"
m0True="3.5562E+17"
m0Err="1.7781E+16"
massTrue="0.39044"	
massErr="1.9522E-2"
dpTrue="4.7397E-8"
dpErr="2.3699E-9"

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

echo "---------"

rm -f Network*

# If we've made it here, everything has been going well!
echo "Titania 3 test 2/2 passed: Network."

#####################################################################
# If we've made it here, everything has been going well!
echo "Titania 3 test passed."

cd ..

exit 0
