#!/bin/bash

# MOPS BINTREE 2 TEST
# (C) CASPER LINDBERG 2019
#
# File purpose:
#       This test is designed to check the generic particle model,
#       BintreePrimary with coordinate tracking. 

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

###################################################################
# CASE A: Test ballistic cluster cluster aggregation (BCCA)
# Obtained with 4 runs 1024 particles 
# SHA-1: 308216377734f9f4f2e655a7cd65243e3bc7cddb
# #################################################################
fname="bintreeA-part.csv"

m0True="1.688E+14"
m0err="8.4E+12"
DcTrue="1.124E-08"	
Dcerr="6.9E-10"
npTrue="59.3"
nperr="3.0"
dpTrue="1.81165148E-09"
dperr="5E-17"

# Run MOPS
echo "Running MOPS for CASE A: BCCA coagulation..."
"$program" -p --flamepp -s "sweep-A.xml" -r "mops-A.inx" -g "gasphase-A.inp" > /dev/null
CheckErr $?

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
CheckValues $m0True $m0 $m0err
Dc=${line1[8]}
echo "Checking Dc..."
CheckValues $DcTrue $Dc $Dcerr
np=${line1[36]}
echo "Checking np..."
CheckValues $npTrue $np $nperr
dp=${line1[38]}
echo "Checking dp..."
CheckValues $dpTrue $dp $dperr

rm -f bintreeA*

# If we've made it here, everything has been going well!
echo "BCCA coagulation passed."

cd ..

exit 0