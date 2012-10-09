#!/bin/bash

# This script regression tests the Camflow Batch Reactor model. (H2-O2 flame)

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

# run camflow
"$program"

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

./checkOutput.pl
result=$?
if (($result!=0))
then
  exit 1
else
  # All tests passed
  echo "All tests passed"
  rm -rf reactionsParsed speciesParsed
  exit 0
fi

