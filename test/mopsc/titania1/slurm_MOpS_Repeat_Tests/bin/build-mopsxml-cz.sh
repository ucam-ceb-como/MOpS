#!/bin/bash

prev=$1
curr=`expr $prev + 1`

# Simulation files
stm="Network(stage"$prev")"

# Run python script to produce params file
python << EOF
from numpy import genfromtxt

# Get data from previous run
chem_txt = genfromtxt('$stm-chem.csv', dtype=None, delimiter=',')
chem_dat = genfromtxt('$stm-chem.csv', dtype=float, delimiter=',')
part_dat = genfromtxt('$stm-part.csv', dtype=float, delimiter=',')

# Assign species concs, op conditions and m0 to variables
num_pnts = len(chem_dat)-1
num_spcs = (len(chem_dat[1])-8)/2
sp_names = chem_txt[0,2:2*num_spcs+1:2]
sp_concs = chem_dat[num_pnts,2:2*num_spcs+1:2]
tot_conc = sum(sp_concs) 
temp_out = chem_dat[num_pnts,2*num_spcs+2]
dens_out = chem_dat[num_pnts,2*num_spcs+4]
pres_out = chem_dat[num_pnts,2*num_spcs+6]
mom0_out = part_dat[num_pnts,4]

# Compute species mol fractions
sp_fracs = sp_concs/tot_conc 

# Write to simple text file
file = open('mops_params.txt','w')
for j in range (0,len(sp_fracs)):
	file.write(str(sp_fracs[j]))
	file.write(',')
file.write(str(temp_out))
file.write(',')
file.write(str(dens_out))
file.write(',')
file.write(str(pres_out))
file.write(',')
file.write(str(num_pnts-1))
file.write(',')
file.write(str(mom0_out))
file.close()


# Write headers to file
file = open('mops_heads.txt','w')
for j in range(0,len(sp_fracs)):
	file.write(sp_names[j])
	file.write(',')
file.close()
EOF

# Read the values from the params file
textline=`tail mops_params.txt`
line1=`echo $textline | tr ',' '\n'`
headline=`tail mops_heads.txt`
line2=`echo $headline | tr ',' '\n'`

# Loop variables
array1=($line2)
i=0
j=0

inp="mops-hm-temp-cz.xml"
temp1="mops-hm-s"$curr".xml"
temp2="tmops-hm-temp.xml"

cp $inp $temp1

for l in $line1
do
    k=$(echo "2*$i;" | bc)
    sloc=`expr $j + 12`
    tloc=41
    ploc=42
    mloc=44

    # Gas species mole fractions from previous reactor
    if [ $i -lt 29 ]; then
        ci=$l
	ni=${array1[$k]}
        eval "sed $sloc's/.*/	    <component id=\""$ni"\">$ci<\/component>/' $temp1" > $temp2	
	eval "sed ' ' $temp2" > "$temp1"
    fi
    # Temperature
    if [ $i -eq 29 ]; then
        t=$l
	eval "sed $tloc's/.*/	    <temperature units=\"K\">$t<\/temperature>/' $temp1" > $temp2
	eval "sed ' ' $temp2" > "$temp1"	
    fi
    # Pressure
    if [ $i -eq 31 ]; then
        p=$l 
	eval "sed $ploc's/.*/	    <pressure units=\"Pa\">$p<\/pressure>/' $temp1" > $temp2
	eval "sed ' ' $temp2" > "$temp1"
    fi
    # Number of timesteps
    if [ $i -eq 32 ]; then
	ntsteps=$l
    fi
    # Number density and ensemble
    if [ $i -eq 33 ]; then 
        m0=$l
	ifile="$stm(0)-SP($ntsteps).ens"
        eval "sed $mloc's/.*/            <file>$ifile<\/file>/' $temp1" > $temp2
        eval "sed ' ' $temp2" > "$temp1"
	mloc=`expr $mloc + 1`
        eval "sed $mloc's/.*/            <m0>$m0<\/m0>/' $temp1" > $temp2
	eval "sed ' ' $temp2" > "$temp1"
    fi
    i=`expr $i + 1`
    j=`expr $j + 1`
done
	
rm -f $temp2 "mops_heads.txt" "mops_params.txt"
