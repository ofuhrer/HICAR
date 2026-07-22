#!/bin/bash -l
#SBATCH --job-name=hicar-hypre-adapter
#SBATCH --partition=debug
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --output=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/hypre_adapter_%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/hypre_adapter_%j.err
set -euo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc
root=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate
hypre=${HYPRE_PREFIX:-/scratch/mch/olifu/icon_hicar/deps/hypre-nvhpc-cuda}
mpicc "$root/tools/hypre_adapter_smoke.c" "$root/src/physics/wind_hypre_adapter.c" \
  -I"$hypre/include" -L"$hypre/lib64" -Wl,-rpath,"$hypre/lib64" -lHYPRE -o "$root/revalidation/hypre_adapter_smoke"
export MPICH_GPU_SUPPORT_ENABLED=1
srun --distribution=block:block --cpu-bind=none bash -c '
  export CUDA_VISIBLE_DEVICES="$SLURM_LOCALID"
  exec /scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/hypre_adapter_smoke
'
