#!/bin/bash
#SBATCH -J DOS
#SBATCH -o ./DOS.stdout%j
#SBATCH -e ./DOS.stderr%j
##SBATCH --mail-type=end
##SBATCH --mail-type=fail
#SBATCH --mem=256000
#SBATCH --ntasks-per-node=48
##SBATCH --reservation=co2
#SBATCH --partition=6342
#### Load the VOTCA module
module purge
module load app/votca
module load app/gromacs/2019.6

#### Define name of the job and working directory (current folder)
job=$SLURM_JOB_NAME
jobid=$SLURM_JOB_ID
workdir=$PWD

ulimit -s unlimited
ulimit -t unlimited

echo "Job execution start: $(date)" >>  $workdir/$job.out
echo "Shared library path: $LD_LIBRARY_PATH" >>  $workdir/$job.out
echo "Slurm Job ID is: ${jobid}" >>  $workdir/$job.out
echo "Slurm Job name is: ${job}" >>  $workdir/$job.out
echo $SLURM_JOB_NODELIST  >> $workdir/$job.out
echo $workdir >> $workdir/$job.out
echo $JOBDIR >> $workdir/$job.out

python calc_IE_parallel.py

joberror=$?

exit $joberror
