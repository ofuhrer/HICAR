if(NOT DEFINED RRTMGP_SOURCE_DIR)
    message(FATAL_ERROR "RRTMGP_SOURCE_DIR is required")
endif()

set(kernel_file
    "${RRTMGP_SOURCE_DIR}/rrtmgp/kernels/accel/mo_gas_optics_rrtmgp_kernels.F90")
if(NOT EXISTS "${kernel_file}")
    message(FATAL_ERROR
        "Pinned RTE-RRTMGP accelerated gas-optics kernel not found: ${kernel_file}")
endif()

file(READ "${kernel_file}" source)
set(original_nvcompiler_guard
"#if ( defined(_CRAYFTN) && _RELEASE_MAJOR <= 14 ) || ( defined(_OPENMP) && defined(__NVCOMPILER) )")
set(fixed_nvcompiler_guard
"#if ( defined(_CRAYFTN) && _RELEASE_MAJOR <= 14 ) || defined(__NVCOMPILER)")
set(original_declarations
"    integer :: i, minl, maxl
    logical(wl) :: mask(:,:)
    real(wp) :: a(:,:)")
set(intent_declarations
"    integer, intent(in)    :: i
    logical(wl), intent(in) :: mask(:,:)
    real(wp), intent(in)    :: a(:,:)
    integer, intent(inout)  :: minl, maxl")

string(FIND "${source}" "${intent_declarations}" already_patched)
if(NOT already_patched EQUAL -1)
    message(STATUS
        "Pinned RTE-RRTMGP minmaxloc argument intents already patched")
else()
    string(FIND "${source}" "${original_declarations}" patch_site)
    if(patch_site EQUAL -1)
        message(FATAL_ERROR
            "Pinned RTE-RRTMGP minmaxloc declarations changed; refusing an unverified patch")
    endif()
    string(REPLACE "${original_declarations}" "${intent_declarations}" source "${source}")
    file(WRITE "${kernel_file}" "${source}")
endif()

# The accelerated source already carries a sequential min/max-location helper
# for compilers whose device MINLOC/MAXLOC implementation is unsafe.  The
# original guard enabled that helper for NVHPC only under OpenMP.  HICAR builds
# the same file with OpenACC, where NVHPC 24.5 otherwise generates a vector-loop
# live-out scalar for each intrinsic.  Select the helper for every NVHPC GPU
# build, independent of the host offload model.
file(READ "${kernel_file}" source)
string(FIND "${source}" "${fixed_nvcompiler_guard}" guard_already_patched)
if(guard_already_patched EQUAL -1)
    string(FIND "${source}" "${original_nvcompiler_guard}" guard_patch_site)
    if(guard_patch_site EQUAL -1)
        message(FATAL_ERROR
            "Pinned RTE-RRTMGP NVHPC minmaxloc guard changed; refusing an unverified patch")
    endif()
    string(REPLACE "${original_nvcompiler_guard}" "${fixed_nvcompiler_guard}" source "${source}")
    file(WRITE "${kernel_file}" "${source}")
endif()

message(STATUS
    "Patched pinned RTE-RRTMGP NVHPC tropopause search to use explicit minmaxloc")
