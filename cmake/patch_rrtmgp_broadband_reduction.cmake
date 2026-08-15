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
    file(WRITE "${kernel_file}" "${source}")
    message(STATUS
        "Patched pinned RTE-RRTMGP broadband reduction for deterministic NVHPC execution")
endif()
