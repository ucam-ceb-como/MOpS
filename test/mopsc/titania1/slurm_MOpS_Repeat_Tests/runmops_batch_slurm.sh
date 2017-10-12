#!/bin/bash
# Author: A. Boje, modified from MoDS script of S. Mosbach
echo 'This wrapper script launches MOpS'
echo 'on the CARES cluster using 20 cores per node.'
echo

# Function to print usage instructions.
function usage {
    echo "Usage: $0 [<args>]"
    echo "Arguments:"
    echo '    -b,--bin <folder>    The MoDS bin folder to be used, default "./bin".'
    echo "    -N,--nodes <number>  The job will be run on <number> nodes, default 1."
}

# Default bin folder.
BinFolder='./bin'

# Default node number.
NumNodes=1

# Default MOpS executable.
MOpSExe='mops-app'

# Parse command line arguments.
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -b|--bin)
	BinFolder="$2"
	shift # past argument
	;;
    -N|--nodes)
	NumNodes=$2
	shift # past argument
	;;
    --help|*)
	usage
	exit
	;;
esac
shift # past argument or value
done

if [[ -d $BinFolder && ! -z $BinFolder ]]
then

usremailadr=$(git config user.email)
echo "Notification emails will be sent to: $usremailadr"
echo '(NB Edit your git config in order to change this.)'
echo

echo "Bin folder:      $BinFolder"
echo "Number of nodes: $NumNodes"
echo

echo 'Submitting job to Slurm...'
sbatch --mail-user=$usremailadr --job-name='MOpS_pfr' --time 20-00:00:00 --nodes=$NumNodes ./auxslurm_batch.sh $BinFolder $(pwd)/$MOpSExe
echo 'Type squeue to watch it.'
echo
echo 'Done.'

else

echo "Bin folder \"$BinFolder\" not found."
usage
exit

fi
