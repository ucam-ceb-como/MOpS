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
    cd $2
fi

echo "sweep.FixedChem"
# run the test
$program

# capture exit value of simulation
simulationResult=$?

if((simulationResult==0))
  then
    echo "Finished sweep.FixedChem"
    echo "========================"
else
  echo "****** sweep.FixedChem failed ******"
  exit $simulationResult
fi
