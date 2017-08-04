#!/bin/bash

# MOPS REGRESS 7 TEST
# (C) WILLIAM MENZ (wjm34) 2012
#
# File purpose:
#       This test is designed to check the reading and writing of 
#       particle ensemble (*.ens) binary files, as well as the
#       serialisation/deserialisation of binary trees for complex
#       particle models (silica in this case).
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
    cd "$2"
    echo "changed directory to $2"
fi

rundir="regress7"
# relative error theshold tolerable
err_threshold="1.0e-5"

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
my \$err=0;
my \$threshold=$err_threshold;
my \$retval=0;
if(\$oldval != 0) {
    \$err=abs(\$oldval-\$newval)/\$oldval;
    if (\$err > \$threshold) {
            print "TEST FAILURE!!! GOT ERROR OF \$err.\n";
            \$retval=2
    }
} 
exit \$retval
EOF
}

cd $rundir

# First run the executables:
echo "Running 1-silica.."
"$exe" -p --strang --ensemble -r "mops-1-silica.inx" -s "sweep-1-silica.xml" > /dev/null
CheckErr $?

# Get M0 from output file
csvline=`tail -1 1-silica-part.csv`
line=`echo $csvline | tr ',' '\n'`
i=0
for l in $line
do
    if [ $i -eq 4 ]; then
    m0=$l
    fi
    i=`expr $i + 1`
done
echo "Found M0 of $m0 #/m3."

# Prepare temporary mops.inx file
inp="mops-1-silica.inx"
temp1="mops-2-silica.inx"
temp2="mops-2-silica-temp.inx"
ifile="1-silica(0)-SP(100).ens"
filenames="2-silica"
timesteps='<time steps="1" splits="1">0.00001<\/time>'
eval "sed '21s/.*/        <m0>$m0<\/m0>/' $inp" > $temp1
eval "sed '22s/.*/        <file>$ifile<\/file>/' $temp1" > $temp2
eval "sed '35s/.*/    $timesteps/' $temp2" > $temp1
eval "sed '62s/.*/    <filename>$filenames<\/filename>/' $temp1" > $temp2
eval "sed '23,27d' $temp2" > "$temp1"
rm -f $temp2

# Run the calculation
echo "Running second calculation..."
"$exe" -p --strang -r $temp1 -s "sweep-2-silica.xml" > /dev/null
CheckErr $?

# Compare -part.csv files... they *should* be identical, except 
# for rounding errors.
echo "Comparing -part.csv files.."
csvline1=`tail -1 1-silica-part.csv`
line1=`echo $csvline1 | tr ',' '\n'`
csvline2=`tail -1 2-silica-part.csv`
line2=`echo $csvline2 | tr ',' '\n'`
array1=($line1)
array2=($line2)

count=${#array1[@]}
for i in `seq 1 $count`
do
    if [ $i -gt 2 ]; then
        a=${array1[$i-1]} 
        b=${array2[$i-1]}
        CheckValues $a $b
        CheckErr $?
    fi
done

# Defer to python to check PSL files (only where available)
type python > /dev/null 2>&1
if [ $? -gt 0 ]; then
    echo "Python not installed.. not a problem."
else
echo "Checking PSL files now..."
python << EOF
f1=open("1-silica-psl(1s).csv", "r"); f1.readline()
f2=open("2-silica-psl(1e-05s).csv", "r"); f2.readline()
data1=[]
data2=[]
for csvline in f1:
    line=csvline.split(",")
    newline=[]
    for l in line:
        newline.append(float(l.strip()))
    data1.append(newline)
f1.close()
for csvline in f2:
    line=csvline.split(",")
    newline=[]
    for l in line:
        newline.append(float(l.strip()))
    data2.append(newline)
f2.close()

# CALCULATE ERRORS
flag=0
err=[]
for l1, l2 in zip(data1, data2):
    errline=[]
    i=0
    for x1, x2 in zip(l1, l2):
        e=0
        if (i!=7):
            if (x1!=0):
                e=abs(x1-x2)/x1
                if (e > $err_threshold):
                    flag+=1
        errline.append(e)
        i+=1
    err.append(errline)
exit(flag)
EOF
CheckErr $?
fi

# If we've made it here, everything has been going well!
echo "All tests passed. :D"
# Remove temporary files
echo "Cleaning files.."
rm -f 1-silica* 2-silica*
rm -f $temp1

cd ..

exit 0
