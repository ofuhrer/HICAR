#!/bin/bash -l
# Four-node Swiss 200 m acceptance gate for the exact HYPRE wind backend.
# The runner owns an isolated namelist/output directory and never alters the
# historical native-solver experiments under case_studies/swiss_200m/runs.
#SBATCH --job-name=hicar-hypre-swiss-200m
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --output=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/hypre_swiss_%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/hypre_swiss_%j.err
set -euo pipefail

[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc \
  netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate
exe="$root/build_gpu_nccl_nocopy/HICAR_gpu"
tester="$root/build_gpu_nccl_nocopy/tests/HICAR-tester"
nml_src="${NML_SRC:-$SCRATCH/icon_hicar/case_studies/swiss_200m/runs/integration_20100101_0030_dynamic_prec4/hicar.nml}"
wrapper="$SCRATCH/icon_hicar/case_studies/swiss_200m/scripts/gpu_rank_wrapper.sh"
run="$root/revalidation/hypre_swiss_${SLURM_JOB_ID:?}"

test -x "$exe" && test -x "$tester" && test -f "$nml_src" && test -x "$wrapper"
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=0 MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
# HYPRE level 2 provides rank-0 Krylov residual progress for this expensive
# national diagnostic; regular 250 m regressions retain the silent default.
export HICAR_HYPRE_PRINT_LEVEL="${HICAR_HYPRE_PRINT_LEVEL:-2}"

# Retain the independent four-GPU halo gate alongside the national solve.
srun -N1 -n4 --ntasks-per-node=4 --gres=gpu:4 --distribution=block:block --cpu-bind=none \
  bash -c 'export ACC_DEVICE_TYPE=nvidia ACC_DEVICE_NUM=0 CUDA_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$1" -v halo_exch' _ "$tester" \
  2>&1 | tee "$run.halo.out"
grep -q 'PASSED' "$run.halo.out"

mkdir -p "$run/output" "$run/restart"
forcing_list_src=$(awk -F"'" '/^[[:space:]]*forcing_file_list[[:space:]]*=/{print $2; exit}' "$nml_src")
test -f "$forcing_list_src"
# HICAR reads this list with Fortran list-directed input: slash-containing
# absolute paths must be quoted, otherwise it can wait for a corrupt filename.
forcing_list="$run/forcing_files.txt"
awk 'NF && $1 !~ /^#/ {gsub(/^"|"$/, "", $0); printf "\"%s\"\n", $0}' "$forcing_list_src" > "$forcing_list"
test -s "$forcing_list"
sed -e "/^[[:space:]]*restart_folder[[:space:]]*=/c\  restart_folder = '$run/restart/'" \
    -e "/^[[:space:]]*output_folder[[:space:]]*=/c\  output_folder = '$run/output/'" \
    -e "/^[[:space:]]*forcing_file_list[[:space:]]*=/c\  forcing_file_list = '$forcing_list'" \
    "$nml_src" > "$run/input.nml"
for asset in NoahmpTable.TBL rrtmgp_support mp_support; do cp -a "$root/run/$asset" "$run/"; done

static_file=$(awk -F"'" '/^[[:space:]]*init_conditions_file[[:space:]]*=/{print $2; exit}' "$run/input.nml")
test -f "$static_file" && test -f "$static_file.ready" && test -f "$forcing_list"
while IFS= read -r forcing; do
  forcing=${forcing#\"}; forcing=${forcing%\"}
  test -f "$forcing" && test -f "$forcing.ready"
done < "$forcing_list"

cd "$run"
srun --distribution=block:block --cpu-bind=none "$wrapper" "$exe" "$run/input.nml" | tee "$run/model.out"
grep -q 'HICAR HYPRE/native operator relative L2=' "$run/model.out"
grep -q 'HICAR HYPRE FGMRES+MGR/AMG:' "$run/model.out"
! grep -q 'HICAR HYPRE .*rejected\|HICAR HYPRE FGMRES+MGR/AMG failed' "$run/model.out"
grep -q 'Timing across all compute images:' "$run/model.out"
test -n "$(find "$run/output" -name '*.nc' -print -quit)"
