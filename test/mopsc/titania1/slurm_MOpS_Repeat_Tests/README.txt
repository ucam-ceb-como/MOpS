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
for downstream reactors using the buildmops script and the 
mops-hm-temp.xml input template. 
- run for multiple random seeds. The python script in the bin can be 
used to compute errors and variances from these runs (probably with 
some tweaking).

3) The bin contains the input files for all cases described above. 


