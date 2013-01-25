#!/bin/bash

# SILICA 2 TEST FOR SWEEP
# (C) WILLIAM J MENZ (WJM34) 2012

# TESTS SILICA MODEL + TRANSITION KERNEL WITH WEIGHTS
# CASE 1: FREE-MOLECULAR KERNEL
#   M0 checked
#
# CASE 2: INTERPARTICLE REACTION ONLY
#   Mass, Si:O ratio and water concentration checked

######################################################################
# GET ARGUMENTS

# First argument is executable
if [ -z "$1" ]; then
    echo "No executable supplied to $0"
    exit 255
else
    program="$1"
fi

# Second argument is working directory (required)
cwd=`pwd`
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
# 1: FREE-MOLECULAR KERNEL
# Rebased 3 Oct 2012, migration to BinTreeSilicaPrimary
# on commit ca583a307d2c3f5c72f940c7ebb3bab9861cf7ac, 512 SPs & 20 runs

# Define true vales
name="M0 Fv dcol dpri sintlevel"
true="2.86E+013 9.66E-012 8.24E-009 7.66E-009 0.919267"
errs="5.86E+011 9.71E-016 6.94E-011 5.38E-011 0.003847"

# Run calculation
echo "Running MOPS for free-molecular kernel case."
"$program" -p --strang -r "mops-fm.inx" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 silica-fm-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
fv=${line1[16]}
dc=${line1[8]}
dp=${line1[46]}
sl=${line1[48]}
results="$m0 $fv $dc $dp $sl"

# Check values
array0=($name)
array1=($true)
array2=($errs)
array3=($results)
count=${#array1[@]}
for i in `seq 1 $count`
do
    echo "Checking for ${array0[$i-1]}."
    CheckValues ${array1[$i-1]} ${array3[$i-1]} ${array2[$i-1]}
    echo "    ..passed"
done

# Clean files
echo "Free-molecular passes."
rm silica-fm*

######################################################################
# 2: SLIP-FLOW (CONTINUUM) KERNEL
# Rebased 3 Oct 2012, migration to BinTreeSilicaPrimary
# on commit ca583a307d2c3f5c72f940c7ebb3bab9861cf7ac, 512 SPs & 20 runs

# Define true vales
name="M0 Fv dcol dpri sintlevel"
true="1.62E+014 4.80E-009 7.50E-008 1.02E-008 0.0790885"
errs="1.09E+013 7.64E-012 3.11E-009 5.97E-011 0.0077061"

# Run calculation
echo "Running MOPS for slip flow kernel case."
"$program" -p --strang -r "mops-sf.inx" > /dev/null
CheckErr $?

# Parse CSV to check values.
csvline1=`tail -1 silica-sf-part.csv`
line1=(`echo $csvline1 | tr ',' '\n'`)
m0=${line1[4]}
fv=${line1[16]}
dc=${line1[8]}
dp=${line1[46]}
sl=${line1[48]}
results="$m0 $fv $dc $dp $sl"

# Check values
array0=($name)
array1=($true)
array2=($errs)
array3=($results)
count=${#array1[@]}
for i in `seq 1 $count`
do
    echo "Checking for ${array0[$i-1]}."
    CheckValues ${array1[$i-1]} ${array3[$i-1]} ${array2[$i-1]}
    echo "    ..passed"
done

# Clean files
echo "Free-molecular passes."
rm silica-sf*

######################################################################
# FINISH UP
echo "All tests passed! :D"
cd "$cwd"
