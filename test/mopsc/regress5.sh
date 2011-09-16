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
rm -f regression5a*


echo "DSA for constant coagulation kernel"
$program -flamepp -p -gp ./regress5/regress5.inp -rr ./regress5/regress5.inx -s ./regress5/regress5.xml -c ./regress5/chem.inp -t ./regress5/therm.dat
echo "Finished simulation"
echo "========================"

# Array of particle numbers - these should be the count of particles of sizes 1-5 in the psl file
# Put a negative values at the start so that the number of particles of size 1 comes at index 1
# These numbers are for a seed of 123 in the Mersenne Twister random number generator
# Analytic solution is 1458 908 565 352 219
referenceValues=(-1 1483 899 531 374 232)
simulatedValues=(-1 -1   -1  -1 -1 -1)

# Grep seems to require the file in unix format, even under cygwin
if((windows==1))
then
    dos2unix "regression5a-psl(3s).csv"
fi

i=1
while ((i <= 5))
do
  simulatedValues[i]=`grep ",$i$" "regression5a-psl(3s).csv" | wc -l`
  #echo "$i $count"
  ((i+=1))
done

echo "Analytic  solution is 1458 908 565 352 219"
echo "Simulated solution is ${simulatedValues[1]} ${simulatedValues[2]} ${simulatedValues[3]} ${simulatedValues[4]} ${simulatedValues[5]}"

i=1
while ((i <= 5 ))
do
  if((simulatedValues[i] != referenceValues[i])) 
    then
      # Regression test has failed; print explanatory message and exit with non zero
      # value showing the size class of the first difference
      echo "Found ${simulatedValues[i]} particles of size $i, when ${referenceValues[i]} expected"
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
  rm -f regression5a*
  exit 0
else
  exit 1
fi

