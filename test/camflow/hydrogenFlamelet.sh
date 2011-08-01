#!/bin/bash

# This script regression tests the Camflow flamelet model. (H2-O2 flame)

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

# run Camflow on the H2-O2 Flamelet (no arguments needed - just have to be in the correct directory)
cd hydrogenFlamelet
$program

# Now do a diff between new profile.dat and originalProfile.dat
# and output the results to file
diff profile.dat profileOriginal.dat > diffResults.txt

# Check to see if the diffResults file is empty
if [-s diffResults.txt]
then
      echo "**************************"
      echo "****** TEST FAILURE ******"
      echo "*****Differences found****"
      echo "**************************"
      exit 1
else
      echo "Test passed"
      exit 0
fi 


# Now change back up to original directory
cd .. 




