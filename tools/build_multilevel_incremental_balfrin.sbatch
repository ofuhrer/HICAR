#!/bin/bash -l
# Rebuild only sources changed since the last clean target-stack configure.
# Use build_wind_operator_audit_balfrin.sbatch for the final clean gate.
#SBATCH --job-name=hicar-ml-incremental
#SBATCH --partition=pp-short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00

set -euo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc \
  cmake/3.24.4-gcc gmake/4.4.1-gcc \
  netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc \
  hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=${HICAR_MULTILEVEL_ROOT:-$SCRATCH/icon_hicar/HICAR-wind-galerkin-multilevel}
build=${HICAR_MULTILEVEL_BUILD:-$root/build_gpu_nccl}
cmake --build "$build" --target HICAR HICAR-tester --parallel "$SLURM_CPUS_PER_TASK"
test -x "$build/HICAR_gpu"
test -x "$build/tests/HICAR-tester"
