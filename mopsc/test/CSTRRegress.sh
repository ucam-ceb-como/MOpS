#!/bin/bash

cd cstrtest
# run mops on a very simple CSTR problem
../../bin/mops_d.exe -gpc -p -diag4 -rr mops.inx -s CSTRRegress.xml -c chem.inp -t therm.dat
if(($?==0)) 
then
  echo "Finished simulation"
else
  echo "Simulation failed"
  exit $?
fi
echo "========================"

# Check to see if the same outlet concentrations for 
# Put in standard values for this type of reaction.
# If the functionality has been changed, the CSTR will not 
#give the same outputs.


dos2unix "J3-chem.csv"
dos2unix "J3-sensi.csv"
./cstrtest.pl
if(($?!=0)) 
  then
    cd ..
    exit $?
fi


rm J3*
# All tests passed
#echo "All tests passed"



exit 0

