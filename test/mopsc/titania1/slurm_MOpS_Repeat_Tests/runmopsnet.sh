#!/bin/bash

mopsapp=$1
nreac=$2
rv=$3
file=RunList.txt
swpfile=sweep-fo-detailed-w3-1.xml

ncstr=10
 
echo "Seed is: " > $file
echo $rv >> $file
echo "DZ, stage: " >> $file
echo "1" >> $file
$mopsapp -p --strang --ensemble -t ../therm.dat -w -r ../mops-hm-s1.xml -c ../chem.inp -s ../$swpfile -e $rv 
echo "Finished." >> $file

if [ $nreac -gt 1 ]; then
        for i in `eval echo {1..$nreac}`  #Number of reactors to do after number 1 above
        do
		j=`expr $i + 1`
		if [ $i -lt $ncstr ]; then 
                echo "DZ, stage: " >> $file
                echo $j >> $file
		        mopsbuild='build-mopsxml-dz.sh'
				cp ../$mopsbuild .
			    $mopsbuild $i
        	    $mopsapp -p --strang --ensemble -t ../therm.dat -r ../mops-hm-s$j.xml -c ../chem.inp -s ../$swpfile -e $rv -w
    	        echo "Finished." >> $file
				rm -f $mopsbuild
		elif [ $i -eq $ncstr ]; then
                echo "WZ, stage: " >> $file
                echo $j >> $file
	            mopsbuild='build-mopsxml-wz.sh'
				cp ../$mopsbuild .
			    $mopsbuild $i
        	    $mopsapp -p --strang --ensemble -t ../therm.dat -r ../mops-hm-s$j.xml -c ../chem.inp -s ../$swpfile -e $rv
    	        echo "Finished." >> $file
				rm -f $mopsbuild
		else
                echo "CZ, stage: " >> $file
                echo $j >> $file
	            mopsbuild='build-mopsxml-cz.sh'
				cp ../$mopsbuild .
			    $mopsbuild $i
        	    $mopsapp -p --strang --ensemble -t ../therm.dat -r ../mops-hm-s$j.xml -c ../chem.inp -s ../$swpfile -e $rv
    	        echo "Finished." >> $file
				rm -f $mopsbuild
		fi
        done 
fi



