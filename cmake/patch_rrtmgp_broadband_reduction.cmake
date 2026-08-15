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
set(original_block
"  integer  :: icol, ilev, igpt
  real(wp) :: scalar ! local scalar version

  !$acc                         parallel loop gang vector collapse(2)
  !$omp target teams distribute parallel do simd          collapse(2)
  do ilev = 1, nlev
    do icol = 1, ncol

      scalar = 0.0_wp

      do igpt = 1, ngpt
        scalar = scalar + spectral_flux(icol, ilev, igpt)
      end do

      broadband_flux(icol, ilev) = factor * scalar
    end do
  end do")
set(private_accumulator_block
"  integer  :: icol, ilev, igpt
  real(wp) :: scalar ! local scalar version

  !$acc                         parallel loop gang vector collapse(2) private(scalar)
  !$omp target teams distribute parallel do simd          collapse(2)
  do ilev = 1, nlev
    do icol = 1, ncol

      scalar = 0.0_wp

      ! NVHPC otherwise implicitly parallelizes this ordered g-point sum as a
      ! reduction inside the already-vectorized column/level loop.  At large
      ! column counts that reduction can retain one warp of stale partial sums
      ! across repeated calls.  Keep one deterministic sum per outer iteration.
      !$acc loop seq
      do igpt = 1, ngpt
        scalar = scalar + spectral_flux(icol, ilev, igpt)
      end do

      broadband_flux(icol, ilev) = factor * scalar
    end do
  end do")
set(deterministic_block
"  integer  :: icol, ilev, igpt

  ! Keep each level on a gang and each column on one vector lane.  A private
  ! scalar accumulator still allowed NVHPC to retain one stale warp across
  ! repeated large-domain calls, so accumulate into the uniquely-owned output
  ! element and explicitly order the spectral loop.
  !$acc                         parallel loop gang
  !$omp target teams distribute parallel do simd          collapse(2)
  do ilev = 1, nlev
    !$acc loop vector
    do icol = 1, ncol

      broadband_flux(icol, ilev) = 0.0_wp

      !$acc loop seq
      do igpt = 1, ngpt
        broadband_flux(icol, ilev) = broadband_flux(icol, ilev) + &
                                      spectral_flux(icol, ilev, igpt)
      end do

      broadband_flux(icol, ilev) = factor * broadband_flux(icol, ilev)
    end do
  end do")

string(FIND "${source}" "${deterministic_block}" already_patched)
if(NOT already_patched EQUAL -1)
    message(STATUS
        "Pinned RTE-RRTMGP deterministic broadband reduction patch already applied")
else()
    string(FIND "${source}" "${original_block}" original_patch_site)
    string(FIND "${source}" "${private_accumulator_block}" private_patch_site)
    if(NOT original_patch_site EQUAL -1)
        string(REPLACE "${original_block}" "${deterministic_block}" source "${source}")
    elseif(NOT private_patch_site EQUAL -1)
        string(REPLACE "${private_accumulator_block}" "${deterministic_block}" source "${source}")
    else()
        message(FATAL_ERROR
            "Pinned RTE-RRTMGP broadband reduction changed; refusing an unverified patch")
    endif()
    file(WRITE "${kernel_file}" "${source}")
endif()

message(STATUS
    "Patched pinned RTE-RRTMGP broadband reduction with explicit level/column ownership")
