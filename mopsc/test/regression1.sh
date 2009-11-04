#!/bin/bash

# run mops on a very simple problem
../bin/mops_d.x -flamepp -p -gp regress1.inp -rr regress1.inx -s regress1.xml -c chem.inp -t therm.dat
echo "Finished simulation"
echo "========================"

# Array of particle numbers - these should be the count of particles of sizes 1,2 and 3 in the psl file
# Put a negative values at the start so that the number of particles of size 1 comes at index 1
# These numbers are for a seed of 123 in the Mersenne Twister random number generator
# Analytic solution is 1053 91 12 2
testValues=(-1 1022 116 6 2 0)

i=1
while ((i <= 5))
do
  count=`grep ",$i$" "regression1-128-10-psl(0.1s).csv" | wc -l`
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

# All tests passed
echo "All tests passed"
rm regression1-128-10*
exit 0

