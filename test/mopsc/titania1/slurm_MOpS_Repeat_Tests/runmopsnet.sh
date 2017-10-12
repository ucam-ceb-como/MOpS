#!/bin/bash

mopsapp=$1
mopsbuild=$2
rv=$3
file=RunList.txt

$mopsapp -p --strang --ensemble -t ../therm.dat -w -r ../mops-hm-s1.xml -c ../chem.inp -s ../sweep-fo-detailed.xml
 -e $rv 

echo "Finished stage:" > $file
echo "1" >> $file 

for i in {1..1..1} #Number of reactors to do after number 1 above
do
	$mopsbuild $i
	j=`expr $i + 1`
	$mopsapp -p --strang --ensemble -t ../therm.dat -w -r mops-hm-s$j.xml -c ../chem.inp -s ../sweep-fo-detailed.xml -e $rv
	echo $j >> $file
done 




