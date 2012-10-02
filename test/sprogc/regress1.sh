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

# run sprogc on a very simple problem
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

if (diff logIN logOUT >/dev/null \
    && diff serializeTest serializeTestOriginal >/dev/null \
    && diff reactionsParsed reactionsParsedOriginal >/dev/null \
    && diff speciesParsed speciesParsedOriginal >/dev/null);
then
  # All tests passed
  echo "All tests passed"
  rm -rf logIN logOUT reactionsParsed serializeTest speciesParsed
  exit 0
else
  echo "*** Test failure ***"
  exit 1
fi

