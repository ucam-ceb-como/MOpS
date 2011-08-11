#!/bin/bash

#Script takes one argument which is the base of the output filenames

# Get a list of files that contain particle lists from the named simulation
filePattern=$1'[0-9]*_psl.csv'
pslFiles=(`ls $filePattern`)

# Copy the header row from one file (the first is the most convenient) into
# the merged file
mergedFile=$1'Merged_psl.csv'
head -n 1 ${pslFiles[0]} > $mergedFile

# Now copy the data rows from the files
for pslFile in ${pslFiles[@]}
do
  echo $pslFile
  tail -n +2 $pslFile >> $mergedFile
done

