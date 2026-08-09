if(NOT DEFINED NOAHMP_SOURCE_DIR)
    message(FATAL_ERROR "NOAHMP_SOURCE_DIR is required")
endif()

set(init_file "${NOAHMP_SOURCE_DIR}/drivers/hicar/NoahmpInitMainMod.F90")
if(NOT EXISTS "${init_file}")
    message(FATAL_ERROR "Pinned Noah-MP HICAR initializer not found: ${init_file}")
endif()

file(READ "${init_file}" source)
set(old_block
"                NoahmpIO%SNOW(I,J)  = 0.0!max(NoahmpIO%SNOW(I,J), 10.0)            ! set SWE to at least 10mm
                NoahmpIO%SNOWH(I,J) = NoahmpIO%SNOW(I,J) * 0.005               ! SNOW in mm and SNOWH in m")
set(new_block
"                ! Preserve caller-provided glacier SWE and snow depth.  The generic
                ! consistency checks above already repair a missing member of the pair.")

string(FIND "${source}" "${new_block}" already_patched)
if(NOT already_patched EQUAL -1)
    message(STATUS "Pinned Noah-MP glacier-snow preservation patch already applied")
else()
    string(FIND "${source}" "${old_block}" patch_site)
    if(patch_site EQUAL -1)
        message(FATAL_ERROR
            "Pinned Noah-MP glacier initialization changed; refusing an unverified patch")
    endif()
    string(REPLACE "${old_block}" "${new_block}" source "${source}")
    message(STATUS "Patched pinned Noah-MP to preserve initialized glacier snow")
endif()

set(old_geometry
"             if ( NoahmpIO%ISNOWXY(I,J) >= 0 ) then
                NoahmpIO%ZSNSOXY(I, 1, J) = -NoahmpIO%DZS(1)
             else
                NoahmpIO%ZSNSOXY(I, 1, J) = NoahmpIO%ZSNSOXY(I, 0, J) - NoahmpIO%DZS(1)
             endif
             !$acc loop seq
             do IZ = 2, NoahmpIO%NSOIL
                NoahmpIO%ZSNSOXY(I, IZ, J) = NoahmpIO%ZSNSOXY(I, IZ-1, J) - NoahmpIO%DZS(IZ)
             enddo")
set(new_geometry
"             ! Reproduce NoahmpSnowInitMain's operation order exactly.  The
             ! algebraically equivalent direct DZS subtraction can differ by
             ! one single-precision ULP and seed a restart-only soil response.
             if ( NoahmpIO%ISNOWXY(I,J) >= 0 ) then
                NoahmpIO%ZSNSOXY(I, 1, J) = NoahmpIO%ZSOIL(1)
             else
                NoahmpIO%ZSNSOXY(I, 1, J) = NoahmpIO%ZSNSOXY(I, 0, J) + NoahmpIO%ZSOIL(1)
             endif
             !$acc loop seq
             do IZ = 2, NoahmpIO%NSOIL
                NoahmpIO%ZSNSOXY(I, IZ, J) = NoahmpIO%ZSNSOXY(I, IZ-1, J) + &
                     (NoahmpIO%ZSOIL(IZ) - NoahmpIO%ZSOIL(IZ-1))
             enddo")

string(FIND "${source}" "${new_geometry}" geometry_already_patched)
if(NOT geometry_already_patched EQUAL -1)
    message(STATUS "Pinned Noah-MP restart geometry patch already applied")
else()
    string(FIND "${source}" "${old_geometry}" geometry_patch_site)
    if(geometry_patch_site EQUAL -1)
        message(FATAL_ERROR
            "Pinned Noah-MP restart geometry changed; refusing an unverified patch")
    endif()
    string(REPLACE "${old_geometry}" "${new_geometry}" source "${source}")
    message(STATUS "Patched pinned Noah-MP restart geometry operation order")
endif()

file(WRITE "${init_file}" "${source}")

# Noah-MP 9e293596 fused the vegetated surface-energy kernels into one
# long-lived OpenACC parallel region.  Its numerical scratch scalars are then
# private per gang while vector lanes use them concurrently.  Reversing only
# that execution-organization commit restores the earlier independent kernels,
# where the same scratch is private per loop iteration.  The equations and
# Noah-MP options are unchanged.
find_package(Git REQUIRED)
set(persistent_region_commit
    "9e29359640362cd5f5215e0810c3bb859b5ebb5d")

function(run_persistent_region_patch direction check_only result_variable error_variable patch_commit)
    set(apply_arguments "apply" "${direction}")
    if(check_only)
        list(APPEND apply_arguments "--check")
    endif()
    list(APPEND apply_arguments "-")
    set(show_arguments "show" "--format=" "--binary" "${patch_commit}")
    if(ARGN)
        list(APPEND show_arguments "--" ${ARGN})
    endif()
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${NOAHMP_SOURCE_DIR}"
                ${show_arguments}
        COMMAND "${GIT_EXECUTABLE}" -C "${NOAHMP_SOURCE_DIR}"
                ${apply_arguments}
        RESULT_VARIABLE patch_result
        ERROR_VARIABLE patch_error
    )
    set(${result_variable} "${patch_result}" PARENT_SCOPE)
    set(${error_variable} "${patch_error}" PARENT_SCOPE)
endfunction()

run_persistent_region_patch("--reverse" TRUE reverse_check reverse_error
                            "${persistent_region_commit}")
if(reverse_check EQUAL 0)
    run_persistent_region_patch("--reverse" FALSE reverse_result reverse_error
                                "${persistent_region_commit}")
    if(NOT reverse_result EQUAL 0)
        message(FATAL_ERROR
            "Failed to restore separate Noah-MP vegetated kernels: ${reverse_error}")
    endif()
    message(STATUS "Restored separate Noah-MP vegetated OpenACC kernels")
else()
    # If the forward patch applies cleanly, the source is already in the
    # deliberately reverted state.  Anything else is an unknown mixed tree.
    run_persistent_region_patch("" TRUE forward_check forward_error
                                "${persistent_region_commit}")
    if(NOT forward_check EQUAL 0)
        message(FATAL_ERROR
            "Pinned Noah-MP vegetated kernel source changed; refusing an unverified patch: "
            "${reverse_error}; ${forward_error}")
    endif()
    message(STATUS "Separate Noah-MP vegetated OpenACC kernels already restored")
endif()

# Noah-MP 9a252eac made the analogous unsafe fusion in the bare-ground
# surface-energy path.  Reverse only the three execution-organization files;
# the same commit also changed unrelated snow, soil-water, and radiation code
# that must remain at the pinned revision.
set(bare_persistent_region_commit
    "9a252eac8ea425fe789d4ff68bca9d3c84465938")
set(bare_persistent_region_files
    "src/SurfaceEnergyFluxBareGroundMod.F90"
    "src/ResistanceBareGroundMostMod.F90"
    "src/ResistanceBareGroundChen97Mod.F90")

run_persistent_region_patch("--reverse" TRUE bare_reverse_check bare_reverse_error
                            "${bare_persistent_region_commit}"
                            ${bare_persistent_region_files})
if(bare_reverse_check EQUAL 0)
    run_persistent_region_patch("--reverse" FALSE bare_reverse_result bare_reverse_error
                                "${bare_persistent_region_commit}"
                                ${bare_persistent_region_files})
    if(NOT bare_reverse_result EQUAL 0)
        message(FATAL_ERROR
            "Failed to restore separate Noah-MP bare-ground kernels: ${bare_reverse_error}")
    endif()
    message(STATUS "Restored separate Noah-MP bare-ground OpenACC kernels")
else()
    run_persistent_region_patch("" TRUE bare_forward_check bare_forward_error
                                "${bare_persistent_region_commit}"
                                ${bare_persistent_region_files})
    if(NOT bare_forward_check EQUAL 0)
        message(FATAL_ERROR
            "Pinned Noah-MP bare-ground kernel source changed; refusing an unverified patch: "
            "${bare_reverse_error}; ${bare_forward_error}")
    endif()
    message(STATUS "Separate Noah-MP bare-ground OpenACC kernels already restored")
endif()
