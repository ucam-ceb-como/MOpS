#!/bin/bash

# MOPS NETWORK3 TEST
# (C) WILLIAM MENZ (wjm34) 2013
#
# This regression test checks the validity of MOPS' gas-phase solver in
# reactor networks.
#
# Overview:
#   All reactors are initialised filled with argon. The only inflow
#   stream composition is 10% SiH4, which reacts according to:
#       SiH4 --> Si + 2H2       (k = 1.0)
#
#   Solutions were generated from analytic solutions, a generic ODE
#   solver and Cantera where possible.
#
# Case A: A four reactor 'triangular' network
# Case B: A four reactor linear netowrk
# Case C: A two reactor network with a recycle loop

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

headers=(AR SIH4 SI TEMP DENSITY)

function CheckTest {
csvline1=`tail -1 "$fname"`
CheckErr $?
line1=(`echo $csvline1 | tr ',' '\n'`)

# Get the mops results in an array
trueVals=($trueList)
mopsVals=(${line1[2]} ${line1[6]} ${line1[8]} ${line1[10]} ${line1[12]})

# loop over each and check
count=${#mopsVals[@]}
for i in `seq 1 $count`
do
    echo "Checking ${headers[$i-1]}..."
    CheckValues ${trueVals[$i-1]} ${mopsVals[$i-1]}
done
}
cd $rundir

# relative error theshold tolerable
err_threshold="0.005"

###################################################################
# NETWORK CASE A: DOUBLE INFLOW, NO LOOPS
#
# |‾‾‾‾|
# | R1 |
# |____|---->|‾‾‾‾|     |‾‾‾‾|
#            | R3 |---->| R4 |---->
# |‾‾‾‾|---->|____|     |____|
# | R2 |
# |____|
#
#  The results in R4 are evaluated.
#  Residence times are : R1: 0.5, R2: 1.0, R3: 0.6, R4: 1.0
# #################################################################
mfile="mops-case-a.inx"
tempfile="temp.inx"
fname="casea(r4)-chem.csv"

# A1: NO REACTIONS
# A1a: CONST TEMP, CONST P
trueList="1.21E-005 1.13E-007 0.0 1000 1.22E-005"
echo "Running case A1a: Double inflow network with no reactions (constT, constP)"
$exe -p -w -r $mfile -c blank.inp > /dev/null
CheckErr $?
CheckTest

# A1b: ADIABATIC, CONST P
sed -e 's/constt="true"/constt="false"/g' $mfile > $tempfile
trueList="1.21E-005 1.13E-007 0.0 1000 1.22E-005"
echo "Running case A1b: Double inflow network with no reactions (adiabatic, constP)"
$exe -p -w -r $tempfile -c blank.inp > /dev/null
CheckErr $?
CheckTest

# A2: SIMPLE REACTIONS, CONST TEMP
# A2a: CONST P
trueList="1.20E-005 6.43E-008 4.69E-008 1000 1.22E-005"
echo "Running case A2a: Double inflow network with reactions (constT, constP)"
$exe -p -w -r $mfile -c chem.inp > /dev/null
CheckErr $?
CheckTest

# A2b: CONST V
sed -e 's/constv="false"/constv="true"/g' $mfile > $tempfile
trueList="1.21E-005 6.57E-008 4.80E-008 1000 1.23E-005"
echo "Running case A2b: Double inflow network with reactions (constT, constV)"
$exe -p -w -r $tempfile -c chem.inp > /dev/null
CheckErr $?
CheckTest

# A3: SIMPLE REACTIONS, ADIABATIC
# A3a: CONST P
sed -e 's/constt="true"/constt="false"/g' $mfile > $tempfile
trueList="1.31E-005 7.79E-008 5.73E-008 913.96 1.33E-005"
echo "Running case A3a: Double inflow network with reactions (adiabatic, constP)"
$exe -p -w -r $tempfile -c chem.inp > /dev/null
CheckErr $?
CheckTest

# A3b: CONST V
sed -i -e 's/constv="false"/constv="true"/g' $tempfile
trueList="1.21E-005 6.56E-008 4.79E-008 838.7 1.23E-005"
echo "Running case A3b: Double inflow network with reactions (adiabatic, constV)"
$exe -p -w -r $tempfile -c chem.inp > /dev/null
CheckErr $?
CheckTest

rm -f casea*

###################################################################
# NETWORK CASE B: Linear network
# Inflow is pure SiH4, reactors initially filled with Ar
#       |‾‾‾‾|    |‾‾‾‾|    |‾‾‾‾|    |‾‾‾‾|
#  ---->| R1 |--->| R2 |--->| R3 |--->| R4 |--->
#       |____|    |____|    |____|    |____|
# 
# All have residence time 0.5
# #################################################################
mfile="mops-case-b.inx"
echo "Running case B: Linear network (constT, constP)"
$exe -p -w -r $mfile -c chem.inp > /dev/null
CheckErr $?

fname="caseb(r1)-chem.csv"
trueList="3.78E-008 6.09E-006 2.02E-006 1000 1.22E-005"
CheckTest

fname="caseb(r2)-chem.csv"
trueList="3.44E-007 3.42E-006 2.81E-006 1000 1.22E-005"
CheckTest

fname="caseb(r3)-chem.csv"
trueList="1.32E-006 2.03E-006 2.95E-006 1000 1.22E-005"
CheckTest

fname="caseb(r4)-chem.csv"
trueList="3.17E-006 1.22E-006 2.60E-006 1000 1.22E-005"
CheckTest

rm -f caseb*

###################################################################
# NETWORK CASE C: RECYCLE LOOPS, AAARGH!
#   A recycle stream of fraction 0.4 is sent back from R2 to R1.
#
#       ----------------------- fR=0.4
#       |                     |
#       -->|‾‾‾‾| f=1 |‾‾‾‾|---
# fIn=0.6  | R1 |---->| R2 |     fOut=0.6
#     ---->|____|     |____|----->
#
#  The results in R2 are evaluated.
#  Residence times are : R1: 1.0, R2: 1.0
# #################################################################
mfile="mops-case-c.inx"
fname="casec(r2)-chem.csv"

# C1: NO REACTIONS
# C1a: CONST TEMP, CONST V
# Note that an example ODE solution is provided in caseb-odesl.py
trueList="1.15E-005 6.99E-007 0.0 1000 1.22E-005"
echo "Running case C1a: Recycle network, no reactions (constT, constV)"
$exe -p -w -r $mfile -c blank.inp > /dev/null
CheckErr $?
CheckTest

# C2: SIMPLE REACTIONS
# C2a: CONST TEMP, CONST P
sed -e 's/constv="true"/constv="false"/g' $mfile > $tempfile
trueList="1.06E-005 1.87E-007 4.66E-007 1000 1.22E-005"
echo "Running case C2a: Recycle network with reactions (constT, constP)"
$exe -p -w -r $tempfile -c chem.inp > /dev/null
CheckErr $?
CheckTest

# C2b: ADIABATIC, CONST V
sed -e 's/constt="true"/constt="false"/g' $mfile > $tempfile
trueList="1.15E-005 1.96E-007 5.02E-007 36.75 1.32E-005"
echo "Running case C2b: Recycle network with reactions (adiabatic, constV)"
$exe -p -w -r $tempfile -c chem.inp > /dev/null
CheckErr $?
CheckTest

# If we've made it here, everything has been going well!
echo "All tests passed. :D"
# Remove temporary files
echo "Cleaning files.."

rm -f casec*
rm -f $tempfile

cd ..
