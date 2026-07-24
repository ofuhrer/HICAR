#!/bin/bash -l
# One-node target-stack gate for the discretely adjoint projection.  The
# projection flag itself must be sufficient to activate the exact hierarchy.
#SBATCH --job-name=hicar-adjoint-unit
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --time=00:10:00
#SBATCH --exclusive

set -euo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc \
  netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc \
  hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=${HICAR_MULTILEVEL_ROOT:-$SCRATCH/icon_hicar/HICAR-wind-galerkin-multilevel}
build=${HICAR_MULTILEVEL_BUILD:-$root/build_gpu_nccl}
tester="$build/tests/HICAR-tester"
validation=${HICAR_ADJOINT_VALIDATION:-$root/validation}
result="$validation/adjoint_projection_unit_${SLURM_JOB_ID:?}.out"
test -x "$tester"
mkdir -p "$validation"

export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export HICAR_WIND_ADJOINT_PROJECTION=1
unset HICAR_WIND_MULTILEVEL
unset HICAR_WIND_OPERATOR_AUDIT
unset HICAR_WIND_OPERATOR_AUDIT_MAX_ITERS

srun -N1 -n4 --ntasks-per-node=4 --gres=gpu:4 --distribution=block:block --cpu-bind=none \
  bash -c 'export ACC_DEVICE_TYPE=nvidia ACC_DEVICE_NUM=0 CUDA_VISIBLE_DEVICES=$SLURM_LOCALID; cd "$1"; exec ./tests/HICAR-tester -v wind_iterative' _ "$build" \
  2>&1 | tee "$result"

for test_name in \
  multilevel_transfer multilevel_device multilevel_device_offset multilevel_halo \
  collective_coarse harmonic_ritz adjoint_projection iter_wind_solve_decomp
do
  grep -q "$test_name \\[PASSED\\]" "$result"
done
grep -q 'HICAR exact Galerkin hierarchy ready:' "$result"
grep -q 'HICAR terminal physical solve:.*status=0' "$result"
diagnostic=$(grep 'adjoint projection diagnostic:' "$result" | tail -n 1)
test -n "$diagnostic"
matrix_relative=$(sed -n 's/.*matrix_relative=[[:space:]]*\\([^[:space:]]*\\).*/\\1/p' <<<"$diagnostic")
constraint_relative=$(sed -n 's/.*constraint_relative=[[:space:]]*\\([^[:space:]]*\\).*/\\1/p' <<<"$diagnostic")
awk -v matrix="$matrix_relative" -v constraint="$constraint_relative" \
  'BEGIN { exit !((matrix + 0) <= 1.0e-5 && (constraint + 0) <= 2.0e-5) }'
if grep -Eq 'wind solve rejected|adjoint projection rejected|native FGMRES\\+line failed' "$result"; then
  echo "Rejected solver state found in $result" >&2
  exit 1
fi

echo "Adjoint GPU unit gate passed: matrix_relative=$matrix_relative constraint_relative=$constraint_relative"
