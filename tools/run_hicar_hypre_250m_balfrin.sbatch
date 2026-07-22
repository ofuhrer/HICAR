#!/bin/bash -l
# Target-stack 250 m regression for the exact HICAR operator with
# AMG-preconditioned FGMRES.  Each invocation owns its run directory.
#SBATCH --job-name=hicar-hypre-250m
#SBATCH --partition=debug
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --time=00:25:00
#SBATCH --exclusive
#SBATCH --output=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/hypre_250m_%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/hypre_250m_%j.err
set -euo pipefail

[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc \
  netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate
# The HYPRE validation build is not installed into bin/.  Run the executable
# produced by the matching CMake build so this regression cannot silently use
# an older native-solver binary.
exe="$root/build_gpu_nccl_nocopy/HICAR_gpu"
tester="$root/build_gpu_nccl_nocopy/tests/HICAR-tester"
nml_src="$SCRATCH/icon_hicar/case_studies/icon_ch1_eps_20260710T18_alps_250m/hicar_wind_design_sensitivity_24h/v2_auto1_n80_top12_s26/input/icon_hicar_250m_v2_auto1_n80_top12_s26.nml"
wrapper="$SCRATCH/icon_hicar/case_studies/swiss_200m/scripts/gpu_rank_wrapper.sh"
run="$root/revalidation/hypre_250m_${SLURM_JOB_ID:?}"

test -x "$exe" && test -x "$tester" && test -f "$nml_src" && test -x "$wrapper"
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=0 MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0

# Transport regression gate, independent of the solver regression below.
srun -N1 -n4 --ntasks-per-node=4 --gres=gpu:4 --distribution=block:block --cpu-bind=none \
  bash -c 'export ACC_DEVICE_TYPE=nvidia ACC_DEVICE_NUM=0 CUDA_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$1" -v halo_exch' _ "$tester" \
  2>&1 | tee "$run.halo.out"
grep -q 'PASSED' "$run.halo.out"

mkdir -p "$run/output" "$run/restart"
sed -e "/^[[:space:]]*restart_folder[[:space:]]*=/c\  restart_folder = '$run/restart/'" \
    -e "/^[[:space:]]*output_folder[[:space:]]*=/c\  output_folder = '$run/output/'" \
    "$nml_src" > "$run/input.nml"
for asset in NoahmpTable.TBL rrtmgp_support mp_support; do cp -a "$root/run/$asset" "$run/"; done

cd "$run"
srun --distribution=block:block --cpu-bind=none "$wrapper" "$exe" "$run/input.nml" \
  | tee "$run/model.out"
grep -q 'HICAR HYPRE FGMRES+vertical-line:' "$run/model.out"
grep -q 'Timing across all compute images:' "$run/model.out"
test -n "$(find "$run/output" -name '*.nc' -print -quit)"
