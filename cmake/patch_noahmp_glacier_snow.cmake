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
