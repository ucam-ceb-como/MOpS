#!/bin/bash
# This script creates the folder structure required by MOpS
# on a node scratch space.
# It is being used by a wrapper script.
# DO NOT EXECUTE THIS DIRECTLY ON THE COMMAND LINE!

cd /scratches/cares005
# cd /scratch

# Create sub-folder named after the user.
crsid=$(id -u -n)
mkdir $crsid 2> /dev/null
cd $crsid

# Create sub-folder named after the job (first argument).
mkdir $1
cd $1

# Copy MOpS bin folder (second argument, full path) into current location.
cp $2 ./ -rf
