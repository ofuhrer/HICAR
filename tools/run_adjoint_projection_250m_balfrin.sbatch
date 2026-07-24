#!/bin/bash -l
# Small 250 m regression gate.  This protects the established 80-level
# SLEVE-2/6 configuration but is not a national-scale conditioning proxy.
#SBATCH --job-name=hicar-adjoint-250m
#SBATCH --partition=normal
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --exclusive
#SBATCH --output=/scratch/mch/olifu/icon_hicar/case_studies/icon_ch1_eps_20260710T18_alps_250m/validation/adjoint_250m_%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/case_studies/icon_ch1_eps_20260710T18_alps_250m/validation/adjoint_250m_%j.err

set -euo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc \
  netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc \
  hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=${HICAR_MULTILEVEL_ROOT:-$SCRATCH/icon_hicar/HICAR-wind-galerkin-multilevel}
case_root=${HICAR_250M_CASE:-$SCRATCH/icon_hicar/case_studies/icon_ch1_eps_20260710T18_alps_250m}
build=${HICAR_MULTILEVEL_BUILD:-$root/build_gpu_nccl}
exe="$build/HICAR_gpu"
nml_src=${HICAR_250M_NML_SOURCE:-$case_root/hicar_wind_design_sensitivity_24h/v2_auto1_n80_top12_s26/input/icon_hicar_250m_v2_auto1_n80_top12_s26.nml}
wrapper=${HICAR_GPU_WRAPPER:-$SCRATCH/icon_hicar/case_studies/swiss_200m/scripts/gpu_rank_wrapper.sh}
run="$case_root/validation/adjoint_250m_${SLURM_JOB_ID:?}"

test -x "$exe" && test -f "$nml_src" && test -x "$wrapper"
grep -Eq '^[[:space:]]*nz[[:space:]]*=[[:space:]]*80([[:space:]]|$)' "$nml_src"
grep -Eq '^[[:space:]]*model_top_height[[:space:]]*=[[:space:]]*12000' "$nml_src"
grep -Eq '^[[:space:]]*height_lowest_level[[:space:]]*=[[:space:]]*15' "$nml_src"
grep -Eq '^[[:space:]]*stretch_fac[[:space:]]*=[[:space:]]*0\\.65' "$nml_src"
grep -Eq '^[[:space:]]*decay_rate_L_topo[[:space:]]*=[[:space:]]*2' "$nml_src"
grep -Eq '^[[:space:]]*decay_rate_S_topo[[:space:]]*=[[:space:]]*6' "$nml_src"
grep -Eq "^[[:space:]]*wind[[:space:]]*=[[:space:]]*'variational solver'" "$nml_src"
grep -Eq '^[[:space:]]*Sx[[:space:]]*=[[:space:]]*\\.True\\.' "$nml_src"
grep -Eq '^[[:space:]]*smooth_wind_distance[[:space:]]*=[[:space:]]*500' "$nml_src"
grep -Eq '^[[:space:]]*RK3[[:space:]]*=[[:space:]]*\\.True\\.' "$nml_src"

mkdir -p "$run/output" "$run/restart"
forcing_list_src=$(awk -F"'" '/^[[:space:]]*forcing_file_list[[:space:]]*=/{print $2; exit}' "$nml_src")
static_file=$(awk -F"'" '/^[[:space:]]*init_conditions_file[[:space:]]*=/{print $2; exit}' "$nml_src")
test -f "$forcing_list_src" && test -s "$forcing_list_src"
test -f "$static_file" && test -f "$static_file.ready"
forcing_list="$run/forcing_files.txt"
awk 'NF && $1 !~ /^#/ {gsub(/^"|"$/, "", $0); printf "\"%s\"\n", $0}' \
  "$forcing_list_src" > "$forcing_list"
while IFS= read -r forcing; do
  forcing=${forcing#\"}
  forcing=${forcing%\"}
  test -f "$forcing" && test -f "$forcing.ready"
done < "$forcing_list"

sed -e "/^[[:space:]]*restart_folder[[:space:]]*=/c\\  restart_folder = '$run/restart/'" \
    -e "/^[[:space:]]*output_folder[[:space:]]*=/c\\  output_folder = '$run/output/'" \
    -e "/^[[:space:]]*forcing_file_list[[:space:]]*=/c\\  forcing_file_list = '$forcing_list'" \
    "$nml_src" > "$run/input.nml"
for asset in NoahmpTable.TBL rrtmgp_support mp_support; do
  cp -a "$root/run/$asset" "$run/"
done

export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export HICAR_WIND_ADJOINT_PROJECTION=1
export HICAR_WIND_RECYCLE_DIM=0
unset HICAR_WIND_MULTILEVEL
unset HICAR_WIND_OPERATOR_AUDIT
unset HICAR_WIND_OPERATOR_AUDIT_MAX_ITERS

cd "$run"
srun --distribution=block:block --cpu-bind=none "$wrapper" "$exe" "$run/input.nml" \
  2>&1 | tee "$run/model.out"

grep -q 'HICAR discretely adjoint wind projection enabled' "$run/model.out"
grep -q 'HICAR exact Galerkin hierarchy ready:' "$run/model.out"
grep -q 'HICAR native FGMRES+line:' "$run/model.out"
grep -q 'HICAR adjoint conservation: relative_Bq=' "$run/model.out"
grep -q 'Timing across all compute images:' "$run/model.out"
awk -f "$root/tools/check_adjoint_projection_log.awk" "$run/model.out"
test -n "$(find "$run/output" -name '*.nc' -print -quit)"

echo "Adjoint 250 m gate passed."
