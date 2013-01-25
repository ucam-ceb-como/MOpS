#!/bin/bash

# MOPS PSR2 TEST
# (C) WILLIAM MENZ (wjm34) 2013
#
# This test includes a Python script for the ODE solution of a very 
# basic population balance system. The script has the following 
# requirements: Python 2.6-2.7, Numpy, Scipy.
# 
# The ODE solutions can be obtained by manipulating the input file
# pinput.py. The conditions for each of the test cases are explained 
# below. The regression test does not actually run the test Python
# scripts, but they are included to demonstrate how a simple ODE 
# solution to the population balance equations could be generated.
#
# Note that the particle process rates are written as short functions of
# time. To keep them as a constant rate, leave them as t*0 + const. The
# tests also assume a residence time of 1.0.
# 
# CASE A: INITIALISED PARTICLES AND OUTFLOW
# All process rates at 0, e.g. K = lambda t: t*0.0 + 0.0
# NINIT = [1.0]
# 
# CASE B: REACTOR FILLED WITH PARTICLES, OUTFLOW + COAG ONLY.
# NINIT = [1.0]
# K = lambda t: t*0.0 + 1.0
# 
# CASE C: REACTOR EMPTY, INCEPTION, SR AND OUTFLOW ONLY.
# NINIT = [0.0]
# S = lambda t: t*0.0 + 1.0
# I = lambda t: t*0.0 + 1.0
# 
# CASE D: REACTOR EMPTY; INFLOW AND OUTFLOW ONLY
# NIN = [1.0]
# 
# CASE E: REACTOR EMPTY; INFLOW, OUTFLOW AND WEIGHTED COAG
# NINIT = [0.0]
# NIN = [1.0]
# K = lambda t: t*0.0 + 1.0


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
    rundir="psr2"
fi

# Set the global relative error tolerance (3%)
err_threshold="0.03"

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
m0=${line1[4]}
echo "Checking M0..."
CheckValues $m0True $m0
fv=${line1[16]}
echo "Checking FV..."
CheckValues $fvTrue $fv
m2=${line1[24]}
echo "Checking M2(Mass)..."
CheckValues $m2True $m2
m3=${line1[26]}
echo "Checking M3(Mass)..."
CheckValues $m3True $m3
}

cd $rundir

###################################################################
# CASE A: REACTOR FILLED WITH PARTICLES, OUTFLOW ONLY.
# #################################################################
# True values obtained from ODE solution.
fname="psrtest-part.csv"
m0True=4.54009277986e-05
fvTrue=4.54009277986e-05
m2True=45.4009277986
m3True=45400.9277986

# Run MOPS
echo "Running MOPS for CASE A..."
$exe -p --strang -s "a-sweep.xml" -r "a-mops.inx" > /dev/null
CheckErr $?

CheckTest

###################################################################
# CASE B: REACTOR FILLED WITH PARTICLES, OUTFLOW + COAG ONLY.
# #################################################################
# True values obtained from ODE solution.
m0True=3.0267588829e-05
fvTrue=4.54007105228e-05
m2True=90.7993867522
m3True=249691.654296

# Run MOPS
echo "Running MOPS for CASE B..."
$exe -p --strang -s "b-sweep.xml" -r "a-mops.inx" > /dev/null
CheckErr $?

CheckTest

###################################################################
# CASE C: INCEPTION + SURFACE REACTION + OUTFLOW
# #################################################################
# True values obtained from ODE solution.
m0True=0.99995460647
fvTrue=1.99945525823
m2True=5992918.1065
m3True=25901213879.0

# Run MOPS
echo "Running MOPS for CASE C..."
$exe -p --strang -s "c-sweep.xml" -r "c-mops.inx" > /dev/null
CheckErr $?

CheckTest

###################################################################
# CASE D: INFLOW AND OUTFLOW ONLY
# #################################################################
# True values obtained from ODE solution.
m0True=1.0
fvTrue=1.0
m2True=1000000.0
m3True=1000000000.0

# Run MOPS
echo "Running MOPS for CASE D..."
$exe -p --strang -s "a-sweep.xml" -r "d-mops.inx" > /dev/null
CheckErr $?

CheckTest

###################################################################
# CASE E: INFLOW, SWA COAG (W3) AND OUTFLOW ONLY
# #################################################################
# True values obtained from ODE solution.
m0True=0.732050779366
fvTrue=1.0
m2True=1999046.08692
m3True=6982310054.91

# Run MOPS
echo "Running MOPS for CASE E..."
$exe -p --strang -s "e-sweep.xml" -r "d-mops.inx" > /dev/null
CheckErr $?

CheckTest

# If we've made it here, everything has been going well!
echo "All tests passed. :D"
# Remove temporary files
rm psrtest*

cd ..

exit 0
