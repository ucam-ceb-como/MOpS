#!/bin/bash

# MOPS NETWORK2 TEST
# (C) WILLIAM MENZ (wjm34) 2013
#
# This regression test compares the gas-phase solver part of the MOPS.
# Solutions are compared to those generated using Cantera. A basic
# gas-phase mechanism is used for simplicity.

# First argument is executable path
exe=$1
if test -z "$exe"
  then
    echo "No executable supplied to $0"
    exit 255
fi

# An optional second argument may specify the working directory
if test -n "$2"
  then
    rundir="$2"
    echo "changed directory to $2"
  else
    rundir="network2"
fi

# Set the global relative error tolerance (3%)
err_threshold="0.01"

# Declare check function
function CheckErr {
    if [ "$1" -gt "0" ]; then
        echo "#############################"
        echo "FAILED POSTPROCESSING."
        echo "#############################"
        exit 1
    fi
}

# Declare check value function (call to Perl)
function CheckValues {
perl << EOF
my \$oldval=$1;
my \$newval=$2;
my \$relerr=$err_threshold;
my \$err=0;
my \$retval=0;
if(\$oldval != 0) {
    \$err=abs(\$oldval-\$newval)/\$oldval;
    if (\$err > \$relerr) {
            print "TEST FAILURE!!! GOT \$newval, wanted \$oldval (relerr \$err).\n";
            \$retval=2
    }
} 
exit \$retval
EOF
CheckErr $?
}

function CheckTest {
csvline1=`tail -1 "$fname"`
line1=(`echo $csvline1 | tr ',' '\n'`)
ar=${line1[2]}
echo "Checking Ar concentration..."
CheckValues $arTrue $ar
siH4=${line1[6]}
echo "Checking SiH4 concentration..."
CheckValues $siH4True $siH4
si=${line1[8]}
echo "Checking Si concentration..."
CheckValues $siTrue $si
}

cd $rundir

###################################################################
# CASE A: 4 REACTOR CHAIN, ALL INITIALISED WITH ARGON, SIH4 INFLOW
#         TO FIRST REACTOR. ALL HAVE RESTIME=0.5s.
# #################################################################

# Run MOPS
echo "Running MOPS for GAS PHASE..!"
$exe -p -w > /dev/null
CheckErr $?

# True values obtained from Cantera solution
fname="chain(r1)-chem.csv"
arTrue=9.36E-011
siH4True=6.09E-006
siTrue=2.03E-006
CheckTest

# True values obtained from Cantera solution
fname="chain(r2)-chem.csv"
arTrue=2.65E-009
siH4True=3.42E-006
siTrue=2.92E-006
CheckTest

# True values obtained from Cantera solution
fname="chain(r3)-chem.csv"
arTrue=2.38E-008
siH4True=2.05E-006
siTrue=3.37E-006
CheckTest

# True values obtained from Cantera solution
fname="chain(r4)-chem.csv"
arTrue=1.17E-007
siH4True=1.28E-006
siTrue=3.60E-006
CheckTest

# If we've made it here, everything has been going well!
echo "All tests passed. :D"
# Remove temporary files
rm -f chain*

cd ..

exit 0
