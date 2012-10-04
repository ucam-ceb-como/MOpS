#!/bin/bash

# SILICA 1 TEST FOR SWEEP
# (C) WILLIAM J MENZ (WJM34) 2012

# TESTS KEY PROCESSES FOR THE SILICA PARTICLE MODEL.
# CASE 1: INCEPTION ONLY
#   M0 checked
#
# CASE 2: INTERPARTICLE REACTION ONLY
#   Mass, number of oxygen and water concentration checked
#
# CASE 3: SURFACE REACTION ONLY
#   Dcol, Mass, number of oxygen and water concentration checked
#
# Note that coagulation and sintering are taken care of by the test
# mops.silica2

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
cwd=`pwd`
if [ -z "$2" ]; then
    echo "No working dir supplied to $0"
    exit 254
else
    wkdir=$2
fi
cd $wkdir

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
# 1: INCEPTION ONLY
# Obtained for M0 with 10 runs (512 SPs)
# hash 76f42c72b130ea0dccd1ff6d7abbd37ca7b37686
true="7.82119E+019"
absErr="4.76906E+016"

# Run calculation
echo "Running MOPS for inception case."
$program -p -strang -rr "mops-incep.inx" -s "sweep-incep.xml" -c "chem.inp" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 silica-incep-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
CheckValues $true $m0 $absErr

# Clean files
echo "Inception passes."
rm silica-incep*

######################################################################
# 2: INTERPARTICLE ONLY
# Obtained with 20 runs (128 SPs)
# hash 76f42c72b130ea0dccd1ff6d7abbd37ca7b37686
trueM="3.32E-008"
errM="0.01E-008"
trueOxygen="9.00E11"
errOxygen="0.009E11"
trueWater="1.66E-013"
errWater="0.01E-13"

# Run calculation
echo "Running MOPS for interparticle case."
$program -p -strang -rr "mops-intp-sr.inx" -s "sweep-intp.xml" -c "chem-intp.inp" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 silica-intp-sr-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
mass=${line1[20]}
CheckValues $trueM $mass $errM
oxygen=${line1[36]}
CheckValues $trueOxygen $oxygen $errOxygen

csvline1=`tail -1 silica-intp-sr-chem.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
water=${line1[20]}
CheckValues $trueWater $water $errWater

# Clean files
echo "Interparticle passes."
rm silica-intp-sr*

######################################################################
# 3: SURFACE REACTION ONLY
# Obtained with 20 runs (128 SPs)
# hash 76f42c72b130ea0dccd1ff6d7abbd37ca7b37686
trueDcol="2.05E-008"
errDcol="0.01E-8"
trueM="9.88E-007"
errM="0.01E-007"
trueOxygen="8.14E12"
errOxygen="0.01E12"
trueWater="1.22E-011"
errWater="0.01E-11"

# Run calculation
echo "Running MOPS for surface reaction case."
$program -p -strang -rr "mops-intp-sr.inx" -s "sweep-sr.xml" -c "chem-intp.inp" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 silica-intp-sr-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
dcol=${line1[8]}
CheckValues $trueDcol $dcol $errDcol
mass=${line1[20]}
CheckValues $trueM $mass $errM
oxygen=${line1[36]}
CheckValues $trueOxygen $oxygen $errOxygen

csvline1=`tail -1 silica-intp-sr-chem.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
water=${line1[20]}
CheckValues $trueWater $water $errWater

# Clean files
echo "SR passes."
rm silica-intp-sr*

echo "All tests passed! :D"
cd $cwd
