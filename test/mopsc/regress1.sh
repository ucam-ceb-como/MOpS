#!/bin/bash


#Absolute path to executable should be supplied as first argument to
#this script.  Script will fail and return a non-zero value
#if no executable specified.
program=$1

if test -z "$program"
  then
    echo "No executable supplied to $0"
    exit 255
fi

# An optional second argument may specify the working directory
if test -n "$2"
  then
    cd "$2"
fi

# run mops on a very simple problem
"$program" --flamepp -p -g regress1/regress1.inp -r regress1/regress1.inx -s regress1/regress1.xml -c regress1/chem.inp -t regress1/therm.dat

# capture exit value of simulation
simulationResult=$?

if((simulationResult==0))
  then
    echo "Finished simulation"
    echo "========================"
else
  echo "****** Simulation failed ******"
  exit 255
fi

# Array of particle numbers - these should be the count of particles of sizes 1,2 and 3 in the psl file
# Put a negative values at the start so that the number of particles of size 1 comes at index 1
# Analytic solution is 1053 91 12 2
referenceValues=(-1 1064 84 9  4  1)
simulatedValues=(-1 -1   -1 -1 -1 -1)

# Grep seems to require the file in unix format, even under cygwin
if((windows==1))
then
    dos2unix "regression1-128-10-psl(0.1s).csv"
fi

i=1
while ((i <= 5))
do
  simulatedValues[i]=`grep ",$i$" "regression1-128-10-psl(0.1s).csv" | wc -l`
  #echo "$i $count"
  ((i+=1))
done

echo "Analytic  solution is 1053 91 12 2"
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
  rm -f regression1-128-10*
  exit 0
else
  exit 1
fi

