!> Shared, model-native atmospheric initialization used by both normal runs
!! and the initialization-only preprocessing driver.
module hicar_initialization_core
    use iso_c_binding, only : c_double
    use ieee_arithmetic, only : ieee_is_finite
    use domain_interface, only : domain_t
    use options_interface, only : options_t
    use icar_constants, only : STD_OUT_PE
    use wind, only : init_winds, update_winds, get_last_projection_diagnostics
    use wind_iterative, only : get_last_wind_solve_diagnostics

    implicit none
    private
    public :: initialize_atmospheric_winds

    real(c_double), parameter :: MATRIX_RELATIVE_TOLERANCE = 1.0e-5_c_double
    real(c_double), parameter :: CONTINUITY_RELATIVE_TOLERANCE = 2.0e-5_c_double

contains

    subroutine initialize_atmospheric_winds(domain, options, project_state)
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        logical, intent(in), optional :: project_state

        logical :: apply_projection

        apply_projection = .true.
        if (present(project_state)) apply_projection = project_state

        call init_winds(domain, options)
        if (apply_projection) then
            call update_winds(domain, options)
            call write_requested_diagnostics()
        endif
    end subroutine initialize_atmospheric_winds


    subroutine write_requested_diagnostics()
        integer :: env_status, path_length
        integer :: solve_status, iterations, unit_number
        real(c_double) :: solve_initial, solve_final, matrix_relative
        real(c_double) :: constraint_initial_norm2, constraint_final_norm2
        real(c_double) :: constraint_relative
        logical :: passed
        character(len=4096) :: path
        character(len=5) :: passed_text

        if (.not. STD_OUT_PE) return
        call get_environment_variable(&
            'HICAR_INITIALIZATION_DIAGNOSTICS', path, length=path_length, status=env_status)
        if (env_status /= 0 .or. path_length <= 0) return
        if (path_length > len(path)) error stop 'HICAR initialization diagnostics path is too long'

        call get_last_wind_solve_diagnostics(&
            solve_status, iterations, solve_initial, solve_final)
        call get_last_projection_diagnostics(&
            constraint_initial_norm2, constraint_final_norm2, constraint_relative)

        matrix_relative = solve_final / max(solve_initial, tiny(1.0_c_double))
        passed = solve_status == 0 .and. ieee_is_finite(matrix_relative) .and. &
                 matrix_relative <= MATRIX_RELATIVE_TOLERANCE .and. &
                 ieee_is_finite(constraint_relative) .and. &
                 constraint_relative >= 0.0_c_double .and. &
                 constraint_relative <= CONTINUITY_RELATIVE_TOLERANCE
        if (passed) then
            passed_text = 'true '
        else
            passed_text = 'false'
        endif

        open(newunit=unit_number, file=path(:path_length), status='replace', action='write')
        write(unit_number,'(A)') '{'
        write(unit_number,'(A)') '  "schema": "hicar-initialization-diagnostics-v1",'
        write(unit_number,'(A)') '  "pressure_operator": "HICAR::domain_obj.adjust_pressure_temp",'
        write(unit_number,'(A)') '  "wind_operator": "HICAR::wind.adjoint_variational_projection",'
        write(unit_number,'(A)') '  "staggering": "HICAR_C_GRID_U_V_W_NATIVE",'
        write(unit_number,'(A,I0,A)') '  "wind_solver_status": ', solve_status, ','
        write(unit_number,'(A,I0,A)') '  "wind_solver_iterations": ', iterations, ','
        write(unit_number,'(A,ES24.16E3,A)') &
            '  "wind_matrix_initial_residual": ', solve_initial, ','
        write(unit_number,'(A,ES24.16E3,A)') &
            '  "wind_matrix_final_residual": ', solve_final, ','
        write(unit_number,'(A,ES24.16E3,A)') &
            '  "wind_matrix_relative_residual": ', matrix_relative, ','
        write(unit_number,'(A,ES24.16E3,A)') &
            '  "mass_continuity_initial_norm2": ', constraint_initial_norm2, ','
        write(unit_number,'(A,ES24.16E3,A)') &
            '  "mass_continuity_final_norm2": ', constraint_final_norm2, ','
        write(unit_number,'(A,ES24.16E3,A)') &
            '  "mass_continuity_relative_residual": ', constraint_relative, ','
        write(unit_number,'(A,A,A)') '  "passed": ', trim(passed_text), ','
        write(unit_number,'(A)') '  "producer_commit": "' // GIT_COMMIT // '"'
        write(unit_number,'(A)') '}'
        close(unit_number)
    end subroutine write_requested_diagnostics

end module hicar_initialization_core
