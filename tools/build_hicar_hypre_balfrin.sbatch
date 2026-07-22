#!/bin/bash -l
#SBATCH --job-name=hicar-hypre-build
#SBATCH --partition=pp-short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/build_hypre_%j.out
#SBATCH --error=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate/revalidation/build_hypre_%j.err

# Rebuild the isolated Balfrin validation checkout with the managed HYPRE
# backend enabled.  The source files are synchronized by the caller.
set -euxo pipefail
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module purge || true
module use "${USER_ENV_ROOT:-/mch-environment/v8}/modules"
module load nvhpc/24.5 cray-mpich-nvhpc/8.1.30 cuda/12.3.0-gcc cmake/3.24.4-gcc gmake/4.4.1-gcc netcdf-c/4.9.2-nvhpc netcdf-fortran/4.6.1-nvhpc hdf5/1.14.3-nvhpc fftw/3.3.10-gcc

root=/scratch/mch/olifu/icon_hicar/HICAR-443f60f1-revalidate
build="$root/build_gpu_nccl_nocopy"
hypre=/scratch/mch/olifu/icon_hicar/deps/hypre-nvhpc-cuda
export HICAR_DEPS_SRC_CACHE="$SCRATCH/icon_hicar/HICAR-gpu-wind-prec3/.hicar-deps-cache-gpu"
ncp=$(nc-config --prefix); nfp=$(nf-config --prefix); ffp=$(pkg-config --variable=prefix fftw3)
mpip=$(dirname "$(dirname "$(command -v mpifort)")")
cuda_target="$CUDA_HOME/targets/x86_64-linux"
nvroot=$(dirname "$(dirname "$(dirname "$(command -v nvfortran)")")")
export NCCL_DIR=$(dirname "$(dirname "$(find "$nvroot/comm_libs" -name libnccl.so -print -quit)")")

test -f "$hypre/lib64/cmake/HYPRE/HYPREConfig.cmake"
cmake -S "$root" -B "$build" -DTEST_DRIVE_SUBPROJECT="$root/vendor/test-drive-src" \
  -DMODE=release -DFC=mpifort -DOPENACC=ON -DNCCL=ON -DGPU_ARCH=cc80 \
  -DHICAR_HYPRE=ON -DHYPRE_DIR="$hypre/lib64/cmake/HYPRE" \
  -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++ \
  -DNETCDF_DIR="$ncp" -DNETCDF_INCLUDES="$ncp/include;$nfp/include" -DNETCDF_INCLUDES_F90="$nfp/include" \
  -DNETCDF_LIBRARIES_C="$ncp/lib/libnetcdf.so" -DNETCDF_LIBRARIES_F90="$nfp/lib/libnetcdff.so" \
  -DFFTW_DIR="$ffp" -DFFTW_INCLUDES="$ffp/include" -DFFTW_LIBRARIES="$ffp/lib/libfftw3.so" \
  -DMPI_DIR="$mpip" -DCUDAToolkit_ROOT="$CUDA_HOME" -DCUFFT_LIBRARY="$cuda_target/lib/libcufft.so" \
  -DCMAKE_BUILD_RPATH="$hypre/lib64;$nfp/lib;$ncp/lib;$NCCL_DIR/lib"
cmake --build "$build" --target HICAR HICAR-tester --parallel "$SLURM_CPUS_PER_TASK"
test -x "$build/HICAR_gpu"
test -x "$build/tests/HICAR-tester"
