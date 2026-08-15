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

message(STATUS
    "Patched pinned RTE-RRTMGP minmaxloc with explicit read-only argument intents")
