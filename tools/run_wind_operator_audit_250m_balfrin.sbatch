#!/bin/bash -l
#SBATCH --job-name=hicar-wind-audit-250m
#SBATCH --partition=debug
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --time=00:25:00
#SBATCH --exclusive
#SBATCH --output=/scratch/mch/olifu/icon_hicar/HICAR-wind-operator-audit/validation/audit_250m_%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/HICAR-wind-operator-audit/validation/audit_250m_%j.err

set -euo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc \
  netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=${HICAR_AUDIT_ROOT:-$SCRATCH/icon_hicar/HICAR-wind-operator-audit}
exe="$root/build_gpu_nccl/HICAR_gpu"
tester="$root/build_gpu_nccl/tests/HICAR-tester"
nml_src="$SCRATCH/icon_hicar/case_studies/icon_ch1_eps_20260710T18_alps_250m/hicar_wind_design_sensitivity_24h/v2_auto1_n80_top12_s26/input/icon_hicar_250m_v2_auto1_n80_top12_s26.nml"
wrapper="$SCRATCH/icon_hicar/case_studies/swiss_200m/scripts/gpu_rank_wrapper.sh"
run="$root/validation/audit_250m_${SLURM_JOB_ID:?}"

test -x "$exe" && test -x "$tester" && test -f "$nml_src" && test -x "$wrapper"
export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0

srun -N1 -n4 --ntasks-per-node=4 --gres=gpu:4 --distribution=block:block --cpu-bind=none \
  bash -c 'export ACC_DEVICE_TYPE=nvidia ACC_DEVICE_NUM=0 CUDA_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$1" -v halo_exch' _ "$tester" \
  2>&1 | tee "$run.halo.out"
grep -q 'PASSED' "$run.halo.out"

mkdir -p "$run/output" "$run/restart"
sed -e "/^[[:space:]]*restart_folder[[:space:]]*=/c\  restart_folder = '$run/restart/'" \
    -e "/^[[:space:]]*output_folder[[:space:]]*=/c\  output_folder = '$run/output/'" \
    "$nml_src" > "$run/input.nml"
for asset in NoahmpTable.TBL rrtmgp_support mp_support; do cp -a "$root/run/$asset" "$run/"; done

export HICAR_WIND_OPERATOR_AUDIT=1
export HICAR_WIND_OPERATOR_AUDIT_FILE="$run/operator_audit.csv"
cd "$run"
srun --distribution=block:block --cpu-bind=none "$wrapper" "$exe" "$run/input.nml" | tee "$run/model.out"
grep -q 'Projection operator equivalence: relative L2=' "$run/model.out"
grep -q 'HICAR calibrated wind operator structural audit:' "$run/model.out"
grep -q 'HICAR wind Krylov audit written to' "$run/model.out"
grep -q 'HICAR native FGMRES+line:' "$run/model.out"
grep -q 'Timing across all compute images:' "$run/model.out"
test -s "$run/operator_audit.csv"
test -n "$(find "$run/output" -name '*.nc' -print -quit)"
