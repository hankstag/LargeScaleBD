#!/bin/bash
# first we ensure a clean running environment:
module purge
#SBATCH --time=1:00:00
#SBATCH --mem=10000
#SBATCH --mail-type=END
#SBATCH --mail-user=hs3048@nyu.edu
#SBATCH --output=slurm_%j.out  
# and ensure we can find the executable:
SRCDIR=$HOME/LargeScaleBD
  
# create a unique directory to run this job in, as per the script above
#RUNDIR=$SCRATCH/local_injective_parameterization/run-${SLURM_JOB_ID/.*}
#mkdir $RUNDIR
  
# By default the script will have started running in the directory we ran sbatch from.
# Let's assume our input file is in the same directory in this example. SLURM
# sets some environment variables with information about the job, including
# SLURM_SUBMIT_DIR which is the directory the job was submitted from. So lets
# go there and copy the input file to the run directory on /scratch:
#cd $SLURM_SUBMIT_DIR
#cp models/$1 $RUNDIR
  
# go to the run directory to begin the run:
cd $SCRATCH/LargeScaleBD/

# load whatever environment modules the executable needs:
module load /share/apps/matlab/2017b
  
# run the executable (sending the contents of my_input_params.inp to stdin)
bash $SRCDIR/run_bd_map.sh /share/apps/matlab/2017b $1
