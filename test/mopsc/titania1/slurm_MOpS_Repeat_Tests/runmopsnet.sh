#!/bin/bash

mopsapp=$1
mopsbuild=$2
rv=$3
nreac=$4
file=RunList.txt
swpfile=sweep-fo-detailed-w3-1.xml
 
echo "Seed is: " > $file
echo $rv >> $file
echo "Stage: " >> $file
echo "1" >> $file

$mopsapp -p --strang --ensemble -t ../therm.dat -w -r ../mops-hm-s1.xml -c ../chem.inp -s ../$swpfile -e $rv 

echo "Finished." >> $file

if [ $nreac -gt 1 ]; then
    for i in {1..$nreac} #Number of reactors to do after number 1 above
    do
	    $mopsbuild $i
	    j=`expr $i + 1`
            echo "Stage: " >> $file
            echo $j >> $file
	    $mopsapp -p --strang --ensemble -t ../therm.dat -w -r mops-hm-s$j.xml -c ../chem.inp -s ../$swpfile -e $rv
    	    echo "Finished." >> $file
    done 
fi



