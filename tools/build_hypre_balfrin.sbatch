#!/bin/bash
#SBATCH --job-name=hypre-cuda-build
#SBATCH --partition=pp-short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/mch/olifu/icon_hicar/deps/hypre-build-%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/deps/hypre-build-%j.err

# Build the pinned HYPRE v2.32.0 source checkout for the Balfrin NVHPC/MPI stack.
# The resulting prefix is intentionally outside the HICAR source tree.
set -euo pipefail

[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc cmake/3.24.4-gcc gmake/4.4.1-gcc

src=/scratch/mch/olifu/icon_hicar/deps/hypre-src/src
build=/scratch/mch/olifu/icon_hicar/deps/hypre-build-nvhpc-cuda
prefix=/scratch/mch/olifu/icon_hicar/deps/hypre-nvhpc-cuda

test -f "$src/CMakeLists.txt"
rm -rf "$build"
cmake -S "$src" -B "$build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_CUDA_COMPILER=nvcc \
  -DHYPRE_WITH_MPI=ON \
  -DHYPRE_WITH_CUDA=ON \
  -DHYPRE_WITH_GPU_AWARE_MPI=ON \
  -DHYPRE_ENABLE_CUDA_STREAMS=ON \
  -DHYPRE_CUDA_SM=80 \
  -DHYPRE_ENABLE_SHARED=ON \
  -DHYPRE_BUILD_TESTS=ON
cmake --build "$build" --parallel "$SLURM_CPUS_PER_TASK"
cmake --install "$build"

test -f "$prefix/lib/libHYPRE.so" -o -f "$prefix/lib64/libHYPRE.so"
