#!/bin/bash

#Script takes one argument which is the base of the output filenames

# Get a list of files that contain particle statistics from the named simulation
filePattern=$1'[0-9]*_partstats.csv'
momentFiles=(`ls $filePattern`)

# Copy the header row from one file (the first is the most convenient) into
# the merged file
mergedFile=$1'Merged_partstats.csv'
head -n 1 ${momentFiles[0]} > $mergedFile

# Now copy the data rows from the files
for momentFile in ${momentFiles[@]}
do
  echo $momentFile
  grep "^[[:digit:]]" $momentFile >> $mergedFile
done

