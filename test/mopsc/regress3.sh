#!/bin/bash

#Path to executable should be supplied as first argument to
#this script.  Script will fail and return a non-zero value
#if no executable specified.
program=$1

#Get rid of any results from earlier runs of this test
rm regression3a*


# run mops on a very simple problem with weighted particles and additive kernel
echo "Weighted particle simulation of additive coagulation kernel"
$program -flamepp -p -gp ./regress3/regress3.inp -rr ./regress3/regress3.inx -s ./regress3/regress3.xml -c ./regress3/chem.inp -t ./regress3/therm.dat
echo "Finished simulation"
echo "========================"

# Array of particle numbers - these should be the count of particles of sizes 1-5 in the psl file
# Put a negative values at the start so that the number of particles of size 1 comes at index 1
# These numbers are for a seed of 123 in the Mersenne Twister random number generator
# Analytic solution is 1427 888 622 459 349
testValues=(-1 1416 903 623 457 323)

# Grep seems to require the file in unix format, even under cygwin
if((windows==1))
then
    dos2unix "regression3a-psl(0.5s).csv"
fi

i=1
while ((i <= 5))
do
  count=`grep ",$i$" "regression3a-psl(0.5s).csv" | wc -l`
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
  rm regression3a*
  exit 0
else
  exit 1
fi

