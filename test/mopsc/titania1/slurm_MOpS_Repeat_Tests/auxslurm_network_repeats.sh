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
#SBATCH --cpus-per-task=20

echo 'Preparing scratch spaces of all nodes involved in the job...'
srun --ntasks-per-node=1 -l ./scratchprep.sh Job$SLURM_JOBID $(pwd)/$1
echo

TheFolder="/scratches/cares005/$(id -u -n)/Job$SLURM_JOBID/$1"

echo "Changing to folder: $(hostname)$TheFolder"
cd $TheFolder
echo

#Uncomment to do OpenMP run
#export OMP_NUM_THREADS=20

echo "Launching MOpS ($2)..."
echo
for t in {1..1..1} #Number of repeats to do
do
	RANDOM=$t
	rv=$RANDOM
	echo "Seed "$t " is " $rv
	dirname=$t
        mkdir $dirname
        cd $dirname
	srun -n1 $3 $2 $4 $rv $t & 
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
