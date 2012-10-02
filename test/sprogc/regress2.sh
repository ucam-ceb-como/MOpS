#!/bin/bash

# REGRESSION TEST 2 FOR SPROG
# (C) WILLIAM J MENZ (WJM34) 2012

# TESTS THE FALLOFF AND THIRD-BODY REACTIONS EXPRESSIONS IN SPROG
# 3 BODY REACTIONS ARE OF THE FORM:
#       SIH4 + M =  H2 + SIH2 + M
#       
#       R = A T^n exp(E/RT) [SiH4] [M]
#       A has units cm3/mol.s
# 
# FALLOFF REACTIONS ARE OF THE FORM
#       SIH4 (+M) =  H2 + SIH2 (+M)
#       
#       MAY USE LINDEMANN FORMULA (LOW) OR TROE FORMULA (TROE)
#       R = k_falloff T^n exp(E/RT) [SiH4]
#
#       k_falloff = k_infty * (P_r/(1+P_r)) * F
#       P_r = k_0 * [M] / k_infty (Reduced pressure)
#       k_0: low pressure rate constant, A has units cm3/mol.s
#       k_infty: high pressure rate constant, A has units 1/s
#       F: Troe fitting formula (F=1 for Lindemann)
#       
#       NOTE THAT AS [M] -> 0, k_falloff -> k_0 * [M]
#                    [M] -> infty, k_falloff -> k_infty

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
    wkdir="$2"
fi

# Define functions for collecting compare data
function GetDataFromFileAtTime {
    fname=$1
    time=",$2,"
    ind=$3
    ans=`grep ${time} $fname | cut -d, -f ${ind}`
    echo $ans
}

# Checks for failed calculations
function CheckErr {
    if [ "$1" -gt "0" ]; then
        echo "FAILED CALCULATION WITH ERR $1."
        exit $1
    fi
}

# Declare check value function (call to Perl)
# relative error theshold tolerable
err_threshold="0.004"   # 0.4% relative deviation allowed
function CheckValues {
perl << EOF
my \$oldval=$1;
my \$newval=$2;
my \$err=0;
my \$threshold=$err_threshold;
my \$retval=0;
if(\$oldval != 0) {
    \$err=abs(\$oldval-\$newval)/\$oldval;
    if (\$err > \$threshold) {
            print "TEST FAILURE!!! ";
            print "WANTED \$oldval, GOT \$newval (ERROR \$err).\n";
            \$retval=2
    }
} 
exit \$retval
EOF
}

# Compares all results for each of the cases against their true values
function CompareResults {
    true=($1)   # True solutions
    this=($2)   # Current solutions
    count=${#true[@]}
    for i in `seq 1 $count`; do
        CheckValues ${true[$i-1]} ${this[$i-1]}
        CheckErr $?
    done
}

#####################################################################
############# MAIN PROGRAM BODY #####################################
#####################################################################
cd "$wkdir"

# Define chem files for the test
chemFiles="chem.3body.inp chem.lindemann.inp chem.troe.inp"

# Define true data for each test
# format: 
#       gp rxn 0 net rate at t=10s, 20s
#       sih4 conc at t=10s, 20s
#       sih2 conc at t=10s, 20s
#          rate@10s rate@20s sih4@10s sih4@20s sih2@10s sih2@20s
true3body="7.55E-04 2.79E-04 9.30E-09 4.48E-09 1.22E-08 1.70E-08"
 trueLind="7.47E-04 2.55E-04 8.73E-09 4.13E-09 1.27E-08 1.73E-08"
 trueTroe="7.38E-04 2.66E-04 9.02E-09 4.38E-09 1.25E-08 1.71E-08"

# Loop over the files and collect their data
for f in $chemFiles; do
    
    echo "Running case $f."
    
    # Run calculation
    "$program" -p -c "$f" > /dev/null
    CheckErr $?
    
    # Use grep and cut to parse cSV
    r10=`GetDataFromFileAtTime "silicon-gp-rates.csv" 10 3`
    r20=`GetDataFromFileAtTime "silicon-gp-rates.csv" 20 3`
    s10=`GetDataFromFileAtTime "silicon-chem.csv" 10 3`
    s20=`GetDataFromFileAtTime "silicon-chem.csv" 20 3`
    l10=`GetDataFromFileAtTime "silicon-chem.csv" 10 5`
    l20=`GetDataFromFileAtTime "silicon-chem.csv" 20 5`
    result="$r10 $r20 $s10 $s20 $l10 $l20"
    #echo $result
    
    # Compare against true solutions
    if [ "$f" == "chem.3body.inp" ]; then
        CompareResults "$true3body" "$result"
        CheckErr $?
    elif [ "$f" == "chem.lindemann.inp" ]; then
        CompareResults "$trueLind" "$result"
        CheckErr $?
    elif [ "$f" == "chem.troe.inp" ]; then
        CompareResults "$trueTroe" "$result"
        CheckErr $?
    else
        echo "TEST FAILURE: Unrecognised file found!"
        exit 1
    fi
done

# If we get here then all tests have passed.
rm silicon*
cd "$cwd"
echo "All tests passed; your falloff reactions are working well! :D"
exit 0
