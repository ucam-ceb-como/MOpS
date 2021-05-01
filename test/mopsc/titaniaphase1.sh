#!/bin/bash

# MOPS TITANIA PHASE 1 TEST
# (C) CASPER LINDBERG 2019
#
# File purpose:
#       Test Gibbs energy model in titania_melting_model 
#	The test case is flame 3 in c4e 238.

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
# Reference values from 20 runs with 512 particles 
# SHA-1: 3531c17e7eed8a8512dbf857a4045083a04d21f8
#####################################################################
fname="Z1-part.csv"

m0True="9.655E+16"
m0err="3.8E+15"
massTrue="1.659E-4"	
masserr="2.4E-6"
TiAnTrue="5.37E+14"
TiAnerr="3.2E+13"
TiRuTrue="7.14E+14"
TiRuerr="2.8E+13"
dpTrue="7.786E-9"
dperr="9.0E-11"

# Run MOPS
echo "Running MOPS for titania phase 1"
"$program" -p --flamepp --endpoint -g "gasphase.csv" > /dev/null
CheckErr $?

csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
echo "Checking M0..."
CheckValues $m0True $m0 $m0err
mass=${line1[20]}
echo "Checking Mass..."
CheckValues $massTrue $mass $masserr
TiAn=${line1[40]}
echo "Checking Ti_An..."
CheckValues $TiAnTrue $TiAn $TiAnerr
TiRu=${line1[48]}
echo "Checking Ti_Ru..."
CheckValues $TiRuTrue $TiRu $TiRuerr
dp=${line1[58]}
echo "Checking Avg. primary diameter..."
CheckValues $dpTrue $dp $dperr

rm -f Z1*

# If we've made it here, everything has been going well!
echo "Titania phase 1 test passed."

cd ..

exit 0