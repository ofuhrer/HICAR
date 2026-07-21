#!/bin/bash
#SBATCH --job-name=hypre-device-smoke
#SBATCH --partition=debug
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/mch/olifu/icon_hicar/deps/hypre-device-smoke-%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/deps/hypre-device-smoke-%j.err

# Four-rank, four-GPU device-memory validation using HYPRE's maintained
# CUDA-native IJ test driver.  This is deliberately independent of HICAR.
set -euo pipefail

[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc

build=/scratch/mch/olifu/icon_hicar/deps/hypre-build-nvhpc-cuda
test -x "$build/test/ij"

export MPICH_GPU_SUPPORT_ENABLED=1
srun --distribution=block:block --cpu-bind=none "$build/test/ij" -laplacian -solver 0
