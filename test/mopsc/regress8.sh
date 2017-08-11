#!/bin/bash

# MOPS REGRESS 8 TEST
# (C) WILLIAM MENZ (wjm34) 2013
#
# File purpose:
#       This test is designed to check the adiabatic PFR gas-phase only 
#       solver. A simple (non-physical!) test problem is used for both
#       a constant pressure and constant volume reactor. Solutions are
#       checked from Cantera results.
#

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
    rundir="regress8"
fi

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
csvline1=`tail -1 "adiabatic-chem.csv"`
line1=(`echo $csvline1 | tr ',' '\n'`)
ar=${line1[2]}
echo "Checking Ar concentration..."
CheckValues $arTrue $ar
siH4=${line1[6]}
echo "Checking SiH4 concentration..."
CheckValues $siH4True $siH4
h2=${line1[4]}
echo "Checking H2 concentration..."
CheckValues $h2True $h2
T=${line1[10]}
echo "Checking Temperature..."
CheckValues $tempTrue $T
rho=${line1[12]}
echo "Checking Density..."
CheckValues $rhoTrue $rho
}

# Change to test directory
cd "$rundir"

# relative error theshold tolerable
err_threshold="0.005"

###################################################################
# CASE A: CONSTANT PRESSURE REACTOR
# #################################################################
tempTrue=465.96
rhoTrue=2.62E-005
siH4True=1.64E-006
arTrue=2.21E-005
h2True=1.62E-006
pwd
# Run mops
echo "Running MOPS for constant pressure..!"
$exe -p -r mops-constp.inx > /dev/null
CheckErr $?

# Check results
CheckTest

###################################################################
# CASE B: CONSTANT VOLUME REACTOR
# #################################################################
tempTrue=197.26
rhoTrue=1.30E-005
siH4True=8.17E-007
arTrue=1.10E-005
h2True=8.03E-007

# Run mops
echo "Running MOPS for constant volume..!"
$exe -p -r mops-constv.inx > /dev/null
CheckErr $?

# Check results
CheckTest

# If we've made it here, everything has been going well!
echo "All tests passed. :D"
# Remove temporary files
echo "Cleaning files.."
rm -f adiabatic*

cd ..

exit 0
