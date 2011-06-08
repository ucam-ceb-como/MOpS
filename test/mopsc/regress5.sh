#!/bin/bash


#Path to executable should be supplied as first argument to
#this script.  Script will fail and return a non-zero value
#if no executable specified.
program=$1

# An optional second argument may specify the working directory
if test -n "$2"
  then
    cd $2
    echo "changed directory to $2"
fi

#Get rid of any results from earlier runs of this test
rm regression5a*


echo "DSA for constant coagulation kernel"
$program -flamepp -p -gp ./regress5/regress5.inp -rr ./regress5/regress5.inx -s ./regress5/regress5.xml -c ./regress5/chem.inp -t ./regress5/therm.dat
echo "Finished simulation"
echo "========================"

# Array of particle numbers - these should be the count of particles of sizes 1-5 in the psl file
# Put a negative values at the start so that the number of particles of size 1 comes at index 1
# These numbers are for a seed of 123 in the Mersenne Twister random number generator
# Analytic solution is 1458 908 565 352 219
testValues=(-1 1428 849 563 343 234)

# Grep seems to require the file in unix format, even under cygwin
if((windows==1))
then
    dos2unix "regression5a-psl(3s).csv"
fi

i=1
while ((i <= 5))
do
  count=`grep ",$i$" "regression5a-psl(3s).csv" | wc -l`
  #echo "$i $count"
  if((count != testValues[i])) 
    then
      # Regression test has failed; print explanatory message and exit with non zero
      # value showing the size class of the first difference
      echo "Found $count particles of size $i, when ${testValues[i]} expected"
      echo "**************************"
      echo "****** TEST FAILURE ******"
      echo "**************************"
      exit $i
  fi
  ((i+=1))
done

if((i==6))
then
  # All tests passed
  echo "All tests passed"
  rm regression5a*
  exit 0
else
  exit 1
fi

