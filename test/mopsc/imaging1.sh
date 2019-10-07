#!/bin/bash

# IMAGING TEST 1 FOR SWEEP
# (C) WILLIAM J MENZ (WJM34) 2012

# TEST THE ABILITY OF SWEEP TO CORRECTLY WRITE OUT POVRAY FILES FOR 
# TEM IMAGING OF PARTICLES
# 
# CASE 1:
#   SPHERICAL PARTICLES
#   initialised particle with dx = 1000000
#   tests that this is correctly written to .pov
# CASE 2:
#   SURFVOL PARTICLES
#   initialised particle with dx = 1000000, surf = 10e-15 m2
#   tests that an aggregate with 22 primaries is written
# CASE 3:
#   BINTREE PRIMARY PARTICLES
#   particles created via inception (not initialised)
#   tests that the correct number (4) of -track.csv and .pov files
#      are written
#   tests povray running the output, where povray is installed

######################################################################
# GET ARGUMENTS

# First argument is executable
if [ -z "$1" ]; then
    echo "No executable supplied to $0"
    exit 255
else
    program=$1
fi

# Second argument is working directory (required)
if [ -z "$2" ]; then
    echo "No working dir supplied to $0"
    exit 254
else
    wkdir="$2"
fi
cd "$wkdir"

######################################################################
# DEFINE FUNCTIONS

# Checks for failed calculations
function CheckErr {
    if [ "$1" -gt "0" ]; then
        echo "FAILED CALCULATION WITH ERR $1."
        exit $1
    fi
}

######################################################################
# 1: SPHERICAL PARTICLES
true="sphere {<0, 0, 0>, 16.8459, 1.0}"

# Run calculation
echo "Running MOPS for spherical primary case."
"$program" -p --strang -r "mops.spherical.inx" -s "sweep.spherical.xml" > /dev/null
CheckErr $?

# Check that the particle given by $true is exactly the same.
ans=`grep -i "sphere" spherical-tem\(0.1s\,\ 0\).pov` 
if [ "$ans" == "$true" ]; then
    echo "Correct POVray particle: $true."
    rm -f spherical*
else
    echo "ERROR: WRONG PARTICLE FROM SPHERICAL PRIMARY."
    exit 1
fi

######################################################################
# 2: SURFVOL PARTICLES
true="22"   # 22 particles obtained with specified surface area & diam

# Run calculation
echo "Running MOPS for spherical primary case."
"$program" -p --strang -r "mops.surfvol.inx" -s "sweep.surfvol.xml" > /dev/null
CheckErr $?

# Check that the particle given by $true is exactly the same.
ans=`grep -i "sphere" surfvol-tem\(0.1s\,\ 0\).pov | wc -l`
if [ "$ans" -eq "$true" ]; then
    echo "Correct POVray particle: $true."
    rm -f surfvol*
else
    echo "ERROR: WRONG PARTICLE FROM SURFVOL PRIMARY."
    exit 2
fi



######################################################################
# 3: BINTREE PRIMARY PARTICLES

# Run calculation
echo "Running MOPS for bintree primary case."
"$program" -p --strang -r "mops.bintree.inx" -s "sweep.bintree.xml" > /dev/null
CheckErr $?

# Check that there are four POV and -track CSV files.
numPovFiles=`ls bintree*.pov -l | wc -l`
numCsvFiles=`ls bintree\(0\)-track*.csv -l | wc -l`
if [ $numPovFiles -eq 4 ] && [ $numCsvFiles -eq 4 ]; then
    echo "Correct number of files found for bintree primary."
else
    echo "ERROR: COULDN'T FIND REQUIRED FILES FROM BINTREE PRIMARY."
    exit 3
fi

# Now try running POVray, if installed.
type povray > /dev/null 2>&1
if [ $? -gt 0 ]; then
    echo "POVray not installed.. skipping test!"
else
    echo "Trying to run POVray with supplied files."
    povray +H600 +W800 -D +iparticle.pov +obintree-output.png > /dev/null 2> /dev/null
    CheckErr $?
    echo "POVray passed!"
fi

# Clean bintree files
rm -f bintree*

echo "All tests passed! :D"
