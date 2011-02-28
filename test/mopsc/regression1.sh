#!/bin/bash

#Choose the windows or linux names for the executable
uname -s | grep --ignore-case CYGWIN 
if(($?==0))
then
    windows=1
else
    windows=0
fi

if((windows==1))
then
	program="../../bin/debug/mops.exe"
else 
	program="../../bin/debug/mops" 
fi

# run mops on a very simple problem
$program -flamepp -p -gp regress1.inp -rr regress1.inx -s regress1.xml -c chem.inp -t therm.dat
echo "Finished simulation"
echo "========================"

# Array of particle numbers - these should be the count of particles of sizes 1,2 and 3 in the psl file
# Put a negative values at the start so that the number of particles of size 1 comes at index 1
# These numbers are for a seed of 123 in the Mersenne Twister random number generator
# Analytic solution is 1053 91 12 2
testValues=(-1 1031 105 9 3 0)

# Grep seems to require the file in unix format, even under cygwin
if((windows==1))
then
    dos2unix "regression1-128-10-psl(0.1s).csv"
fi

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
  fi
  ((i+=1))
done

if((i==6))
then
  # All tests passed
  echo "All tests passed"
  rm regression1-128-10*
  exit 0
else
  exit 1
fi

