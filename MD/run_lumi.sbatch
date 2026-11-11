#!/bin/bash
# LUMI-G job for ONE replica: nvt -> npt1..npt4 -> 100 ns prod (scripts/run_md.sh logic).
# One MI250X GCD + 7 CPU cores on small-g. Submit from MD/ with the replica dir as --chdir:
#   sbatch --chdir=state1_lumen-open/rep1 --job-name=s1_rep1 run_lumi.sbatch
# or simply ./submit_lumi.sh for all 6. The run is restartable: resubmitting the same
# replica skips finished stages and continues an interrupted one from its checkpoint.
#SBATCH --account=project_<number>
#SBATCH --partition=small-g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%x-%j.out

export EBU_USER_PREFIX=/project/project_<number>/EasyBuild
module load LUMI/25.03 partition/G
module load GROMACS/2026.0-cpeAMD-25.03-VkFFT-rocm-hip

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

echo "job $SLURM_JOB_ID on $(hostname) in $PWD, GPU(s): ${ROCR_VISIBLE_DEVICES:-?}"
module list 2>&1 | grep -i gromacs || true

export GMX="srun gmx_mpi"
export MDRUN_FLAGS="-ntomp $SLURM_CPUS_PER_TASK -pin on -nb gpu -pme gpu -bonded gpu -update gpu"

./run_md.sh
