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
# CASE B: Test surface growth (with coagulation)
# Obtained with 2 runs 1024 particles 
# SHA-1: babede7cc254609b55bd760bb33826992630621a
# #################################################################
fname="bintreeB-part.csv"

m0True="3.122E+14"	
m0err="2.15e+13"
DcTrue="1.114E-08"
Dcerr="2.9E-10"
npTrue="11.34"	
nperr="0.56"
dpTrue="5.877E-09"
dperr="3E-12"

# Run MOPS
echo "Running MOPS for CASE B: surface growth..."
"$program" -p --flamepp -s "sweep-B.xml" -r "mops-B.inx" -g "gasphase-B.inp" > /dev/null
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

rm -f bintreeB*

# If we've made it here, everything has been going well!
echo "Surface growth and coagulation passed."

cd ..

exit 0