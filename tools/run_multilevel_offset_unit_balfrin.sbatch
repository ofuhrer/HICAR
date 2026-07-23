#!/bin/bash -l
#SBATCH --job-name=hicar-ml-offset-unit
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00

set -euo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc \
  netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc \
  hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=${HICAR_MULTILEVEL_ROOT:-$SCRATCH/icon_hicar/HICAR-wind-galerkin-multilevel}
build="$root/build_gpu_nccl"
result=${HICAR_MULTILEVEL_RESULT:-$root/validation/multilevel_offset_unit_${SLURM_JOB_ID:?}.out}
test -x "$build/tests/HICAR-tester"

export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export HICAR_WIND_MULTILEVEL=0

srun -N1 -n1 --ntasks-per-node=1 --gres=gpu:1 --cpu-bind=none \
  bash -c 'export ACC_DEVICE_TYPE=nvidia ACC_DEVICE_NUM=0 CUDA_VISIBLE_DEVICES=$SLURM_LOCALID; cd "$1"; exec ./tests/HICAR-tester wind_iterative' _ "$build" \
  2>&1 | tee "$result"

grep -q 'multilevel_device_offset \[PASSED\]' "$result"
