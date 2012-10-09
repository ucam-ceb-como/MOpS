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
"$program" chemkinReader/chem.inp chemkinReader/therm.dat chemkinReader/tran.dat
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

#Compare the output files to the reference output
diff --brief speciesParsed chemkinReader/speciesParsedOriginal
speciesCompResult=$?
diff --brief reactionsParsed chemkinReader/reactionsParsedOriginal
reactionsCompResult=$?

if((speciesCompResult != 0))
  then
    exit 1
fi

if((reactionsCompResult != 0))
  then
    exit 2
fi

# All tests passed if we get this far
echo "All tests passed"
rm -rf reactionsParsed speciesParsed
exit 0
