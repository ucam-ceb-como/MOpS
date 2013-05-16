#!/bin/bash

# TITANIA 1 TEST FOR SWEEP
# (C) WILLIAM J MENZ (WJM34) 2012

# TESTS KEY PROCESSES FOR TITANIA SURFACE REACTION
# CASE 1: FIRST ORDER IN OXYGEN
#   Same surface reaction as in Preprint 99
#
# CASE 2: SIMPLE ELEY-RIDEAL MODEL
#   Pseudo steady state concentration of surface active sites
#   (similar to ABF soot model implementation)
#
# CASE 3: DETAILED ELEY-RIDEAL MODEL
#   Use a bigger chemical type space to track active sites

######################################################################
# GET ARGUMENTS

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

######################################################################
# 1: FIRST ORDER IN OXYGEN
trueDcol="2.49e-7"
absErrDcol="0.01e-7"
trueTiCl4="6.39e-9"
absErrTiCl4="0.01e-9"

# Run calculation
echo "Running MOPS for FIRST ORDER IN OXYGEN case."
"$program" -p --strang -s "sweep-fo.xml" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 Z1-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
val=${line1[8]}
CheckValues $trueDcol $val $absErrDcol

csvline1=`tail -1 Z1-chem.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
val=${line1[2]}
CheckValues $trueTiCl4 $val $absErrTiCl4

# Clean files
echo "FIRST ORDER IN OXYGEN passes."
rm Z1*

######################################################################
# 2: SIMPLE ELEY-RIDEAL MODEL
trueDcol="3.02e-7"
absErrDcol="0.01e-7"
trueTiCl4="4.34e-7"
absErrTiCl4="0.01e-7"

# Run calculation
echo "Running MOPS for SIMPLE ELEY-RIDEAL MODEL case."
"$program" -p --strang -s "sweep-eley-simple.xml" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 Z1-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
val=${line1[8]}
CheckValues $trueDcol $val $absErrDcol

csvline1=`tail -1 Z1-chem.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
val=${line1[2]}
CheckValues $trueTiCl4 $val $absErrTiCl4

# Clean files
echo "SIMPLE ELEY-RIDEAL MODEL passes."
rm Z1*


######################################################################
# 3: DETAILED ELEY-RIDEAL MODEL
trueDcol="4.52e-9"
absErrDcol="0.01e-9"
trueTiCl4="4.37e-7"
absErrTiCl4="0.01e-7"

# Run calculation
echo "Running MOPS for DETAILED ELEY-RIDEAL MODEL case."
"$program" -p --strang -s "sweep-detailed.xml" -r "mops-detailed.inx" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 Z1-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
val=${line1[8]}
CheckValues $trueDcol $val $absErrDcol

csvline1=`tail -1 Z1-chem.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
val=${line1[2]}
CheckValues $trueTiCl4 $val $absErrTiCl4

# Clean files
echo "DETAILED ELEY-RIDEAL MODEL passes."
rm Z1*

echo "All tests passed! :D"
