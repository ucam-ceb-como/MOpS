12 October 2017
Astrid Boje, aab64@cam.ac.uk

Contents of this folder
========================================================================
This folder contains the scripts used to run MOpS with slurm on the 
CARES cluster, in different configurations. 

1) The files with batch/cstr in their name are set up to 
test different settings in the sweep input files for a single reactor. 

2) The files tagged network_repeats can be used to 
- run for multiple reactors in a network, constructing the input files 
for downstream reactors using the buildmops-x.sh scripts and the 
mops-hm-temp-x.xml input templates where x=dz, wz, or cz for the dosing
zone (CSTR), working zone (batch) and cooling zone (batch) respectively. 
- run for multiple random seeds. The python script in the bin can be 
used to compute errors and variances from these runs and combine the 
PSL files.

3) The bin contains the input files for all cases described above. 

4) The scripts with darwin in their names are set up to run the 
tests in (1) and (2) above on the Cambridge Darwin cluster. 


For building MOpS on Darwin (can add to bashrc)
========================================================================
# Load consistent gcc and boost (it also works to build boost v1.47 
# from vienna with gcc 4.8.1
#module load gcc/5.2.0
#module load boost/1.62-python2.7.10-gcc5.2.0

# Add gcc to path
export CC=/usr/local/Cluster-Apps/gcc/5.2.0/bin/gcc
export CXX=/usr/local/Cluster-Apps/gcc/5.2.0/bin/g++

# Add boost libraries to path
export LD_LIBRARY_PATH=/usr/local/Cluster-Apps/boost/1.62/python2.7.10-gcc5.2.0/stage/lib:$LD_LIBRARY_PATH
export INCLUDE=//usr/local/Cluster-Apps/boost/1.62/python2.7.10-gcc5.2.0:$INCLUDE
export LIBRARY_PATH=/usr/local/Cluster-Apps/boost/1.62/python2.7.10-gcc5.2.0/stage/lib:$LIBRARY_PATH
