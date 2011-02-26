#!/bin/bash
#Run brush with the arguments given to this script
#and collect some information on the performance.

#All output will be written to out.txt in the directory
#from which the command is issued.

date > out.txt
echo $HOST >> out.txt
echo $@
echo $@ >> out.txt
echo " "
tail -n 25 /proc/cpuinfo >>out.txt

#Use GNU time rather than the shell built in
/usr/bin/time -a $@ >> out.txt 2>&1

if (($? != 0))
  then
    echo "*********************************"
    echo "******** command failed *********"
    echo "*********************************"
fi
 
