if(NOT DEFINED RRTMGP_SOURCE_DIR)
    message(FATAL_ERROR "RRTMGP_SOURCE_DIR is required")
endif()

set(kernel_file
    "${RRTMGP_SOURCE_DIR}/rte/kernels/accel/mo_rte_solver_kernels.F90")
if(NOT EXISTS "${kernel_file}")
    message(FATAL_ERROR
        "Pinned RTE-RRTMGP accelerated solver kernel not found: ${kernel_file}")
endif()

file(READ "${kernel_file}" source)
set(old_directive
"  real(wp) :: scalar ! local scalar version

  !$acc                         parallel loop gang vector collapse(2)
  !$omp target teams distribute parallel do simd          collapse(2)")
set(new_directive
"  real(wp) :: scalar ! local scalar version

  !$acc                         parallel loop gang vector collapse(2) private(scalar)
  !$omp target teams distribute parallel do simd          collapse(2)")
set(old_block
"      scalar = 0.0_wp

      do igpt = 1, ngpt
        scalar = scalar + spectral_flux(icol, ilev, igpt)")
set(new_block
"      scalar = 0.0_wp

      ! NVHPC otherwise implicitly parallelizes this ordered g-point sum as a
      ! reduction inside the already-vectorized column/level loop.  At large
      ! column counts that reduction can retain one warp of stale partial sums
      ! across repeated calls.  Keep one deterministic sum per outer iteration.
      !$acc loop seq
      do igpt = 1, ngpt
        scalar = scalar + spectral_flux(icol, ilev, igpt)")

string(FIND "${source}" "${new_directive}" directive_already_patched)
if(directive_already_patched EQUAL -1)
    string(FIND "${source}" "${old_directive}" directive_patch_site)
    if(directive_patch_site EQUAL -1)
        message(FATAL_ERROR
            "Pinned RTE-RRTMGP broadband accumulator directive changed; refusing an unverified patch")
    endif()
    string(REPLACE "${old_directive}" "${new_directive}" source "${source}")
endif()

string(FIND "${source}" "${new_block}" already_patched)
if(NOT already_patched EQUAL -1)
    message(STATUS
        "Pinned RTE-RRTMGP deterministic broadband reduction patch already applied")
else()
    string(FIND "${source}" "${old_block}" patch_site)
    if(patch_site EQUAL -1)
        message(FATAL_ERROR
            "Pinned RTE-RRTMGP broadband reduction changed; refusing an unverified patch")
    endif()
    string(REPLACE "${old_block}" "${new_block}" source "${source}")
endif()

file(WRITE "${kernel_file}" "${source}")
message(STATUS
    "Patched pinned RTE-RRTMGP broadband reduction with a private ordered accumulator")
