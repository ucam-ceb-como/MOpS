#!/bin/bash

# MOPS NETWORK1 TEST
# (C) WILLIAM MENZ (wjm34) 2013
#
# This regression tests checks for integrity of the population balance
# part of the MOPS networking code. ODE solutions can be obtained by
# hacking mops.psr2's Python source.
# 
# CASE A: INITIALISED PARTICLES AND IN/OUTFLOW
# All process rates at 0, e.g. K = lambda t: t*0.0 + 0.0
# NINIT = [1.0]

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
    rundir="network1"
fi

# An optional 3rd argument to specify the cases to run.
if test -n "$3"
  then
    cases="$3"
  else
    cases="a b c d"
fi

# Set the global relative error tolerance (3%)
err_threshold="0.03"

# 1:match, 0:failed
function contains {
    [[ $1 =~ $2 ]] && echo 1 || echo 0
}

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
# CASE A: INITIALISED PARTICLES, PARTICLE INFLOW, OUTFLOW, SWA COAG
#         (ONE REACTOR ONLY)
# #################################################################
# True values obtained from ODE solution.
if [ $(contains "$cases" "a") -eq 1 ]; then
    fname="casea(r1)-part.csv"

    m0True=0.732050815536
    fvTrue=2.7999704928
    m2True=20037817.4839
    m3True=233396567499.0

    # Run MOPS
    echo "Running MOPS for CASE A..."
    $exe -p --strang -s "sweep-coag.xml" -r "a-mops.inx" -w > /dev/null
    CheckErr $?

    CheckTest
    
    rm casea*
fi

###################################################################
# CASE B: PARTICLE INFLOW, OUTFLOW, NO COAG
#         (TWO REACTORS: Inflow -> R1 -> R2)
# #################################################################
# True values obtained from ODE solution (of the second reactor)
if [ $(contains "$cases" "b") -eq 1 ]; then
    fname="caseb(r2)-part.csv"

    m0True=0.999500562156
    fvTrue=2.79860157404
    m2True=12193906.8583
    m3True=65167436652.6

    # Run MOPS
    echo "Running MOPS for CASE B..."
    $exe -p --strang -s "sweep-nocoag.xml" -r "b-mops.inx" -w > /dev/null
    CheckErr $?

    CheckTest
    
    rm caseb*
fi

###################################################################
# CASE C: INERT INFLOW, OUTFLOW, INCEP, SWA COAG, SURFACE GROWTH
#         (TWO REACTORS: Inflow -> R1 -> R2)
# #################################################################
# True values obtained from ODE solution (of the second reactor)
if [ $(contains "$cases" "c") -eq 1 ]; then
    fname="casec(r2)-part.csv"

    m0True=1.10029111766
    fvTrue=3.18909564028
    m2True=17300205.2357
    m3True=150797280051.0

    # Run MOPS
    echo "Running MOPS for CASE C..."
    $exe -p --strang -s "sweep.xml" -r "c-mops.inx" -w > /dev/null
    CheckErr $?

    CheckTest
    
    rm cased*
fi

###################################################################
# CASE D: 'TRIANGULAR' NETWORK WITH NO PARTICLE PROCESSES
#
#   inf--->R1--->|‾‾‾‾|
#                | R3 |---->
#   inf--->R2--->|____|
#
#   with M0(inf) = 1.0 #/m3, all with RT = 1.0s
# #################################################################
# True values obtained from ODE solution (of the second reactor)
if [ $(contains "$cases" "d") -eq 1 ]; then
    fname="cased(r3)-part.csv"
    m0True=0.9595715909
    fvTrue=0.9595715909
    m2True=959571.590851
    m3True=959571590.851

    # Run MOPS
    echo "Running MOPS for CASE D..."
    $exe -p --strang -s "sweep-nocoag.xml" -r "d-mops.inx" -w > /dev/null
    CheckErr $?

    CheckTest
    
    rm cased*
fi

# If we've made it here, everything has been going well!
echo "All tests passed. :D"

cd ..

exit 0
