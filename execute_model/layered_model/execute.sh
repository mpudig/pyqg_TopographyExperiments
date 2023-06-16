#!/bin/bash

#SBATCH --job-name=<JOBNAME>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --time=0:30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=mp6191@nyu.edu
#SBATCH --output=slurm_%j.out

# Purge modules to be safe
module purge

# Activate singuality and conda environments, then run script
singularity exec \
	--overlay /scratch/mp6191/pyqg_dev_env/overlay-15GB-500K.ext3:ro /scratch/work/public/singularity/cuda11.2.2-cudnn8-devel-ubuntu20.04.sif \
	/bin/bash -c "source /ext3/env.sh;
python -u /scratch/mp6191/pyqg_expts/<EXPTNAME>/setup/execute_model/layered_model/execute_layered_model.py"
