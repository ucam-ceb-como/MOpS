#!/bin/bash
# Slurm job submission script for MOpS.
# It is being used by a wrapper script.
# DO NOT EXECUTE THIS DIRECTLY ON THE COMMAND LINE!
# DO NOT SUBMIT THIS DIRECTLY TO SLURM!

#SBATCH -p cares.cluster                     # partition (queue)
#SBATCH --mem 32768                          # real memory required per node [MB]
#SBATCH --time 7-0:0                         # wall-time limit [D-HH:MM]
#SBATCH --output slurm.%u.%N.%j.stdout.txt   # (%u,%N,%j)=(user, node, job allocation #)  
#SBATCH --error slurm.%u.%N.%j.errout.txt    #
#SBATCH --mail-type=END,FAIL                 # notifications for job done & fail

echo 'Preparing scratch spaces of all nodes involved in the job...'
srun --ntasks-per-node=1 -l ./scratchprep.sh Job$SLURM_JOBID $(pwd)/$1
echo

TheFolder="/scratches/cares005/$(id -u -n)/Job$SLURM_JOBID/$1"

echo "Changing to folder: $(hostname)$TheFolder"
cd $TheFolder
echo

echo "Launching MOpS ($2)..."
echo
for t in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'
do
	dirname=swa$t
        mkdir $dirname
        cd $dirname
	srun -n1 --exclusive $2 -p --strang --ensemble -t ../therm.dat -r ../mops-batch.xml -c ../chem.inp -s ../sweep-fo-detailed-w3-$t.xml &
        cd ..
done
for t in '1' '2' '3' '4' '5'
do
	dirname=dsa$t
        mkdir $dirname
        cd $dirname
	srun -n1 --exclusive $2 -p --strang --ensemble -t ../therm.dat -r ../mops-batch.xml -c ../chem.inp -s ../sweep-fo-detailed-$t.xml &
        cd ..
done
wait
cd ../../

jobfolder="Job$SLURM_JOBID"
homedir="/home/userspace/aab64/MOpS-suite/EB-branch/mops-c-Git/test/mopsc/titania1/"
mv $jobfolder $homedir

echo
echo 'Slurm job diagnostics:'
sacct --job $SLURM_JOBID --format "JobName,Submit,Elapsed,AveCPU,CPUTime,UserCPU,TotalCPU,NodeList,NTasks,AveDiskRead,AveDiskWrite"
