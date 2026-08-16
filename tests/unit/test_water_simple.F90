module test_water_simple

    use module_water_simple, only : simple_water_flux
    use mod_wrf_constants, only : cp, XLV
    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none
    private

    public :: collect_water_simple_suite

contains

    subroutine collect_water_simple_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("exchange_velocity_units", test_exchange_velocity_units), &
            new_unittest("evaporation_and_condensation", test_evaporation_and_condensation) &
        ]
    end subroutine collect_water_simple_suite

    subroutine test_exchange_velocity_units(error)
        type(error_type), allocatable, intent(out) :: error
        real, parameter :: rho = 1.0, chs = 0.02, delta_t = 2.0
        real, parameter :: qsat = 0.008, qair = 0.004
        real :: sensible, evaporation, latent

        call simple_water_flux(chs, rho, 282.0, 280.0, qsat, qair, &
                               sensible, evaporation, latent)

        call check(error, abs(sensible-rho*cp*(1.0+0.8*qair)*chs*delta_t) < 1.0e-5, &
                   "sensible flux must use rho*moist_cp*CHS*dT")
        if (allocated(error)) return
        call check(error, abs(evaporation-rho*chs*(qsat-qair)) < 1.0e-10, &
                   "moisture flux must use rho*CHS*dq without another wind factor")
        if (allocated(error)) return
        call check(error, abs(latent-XLV*evaporation) < 1.0e-4, &
                   "latent heat must be XLV times moisture flux")
    end subroutine test_exchange_velocity_units

    subroutine test_evaporation_and_condensation(error)
        type(error_type), allocatable, intent(out) :: error
        real :: sensible, evaporation, latent

        call simple_water_flux(0.01, 1.1, 280.0, 280.0, 0.007, 0.005, &
                               sensible, evaporation, latent)
        call check(error, evaporation > 0.0 .and. latent > 0.0, &
                   "undersaturated air must produce upward evaporation")
        if (allocated(error)) return

        call simple_water_flux(0.01, 1.1, 280.0, 280.0, 0.005, 0.007, &
                               sensible, evaporation, latent)
        call check(error, evaporation < 0.0 .and. latent < 0.0, &
                   "supersaturated air must permit downward condensation")
    end subroutine test_evaporation_and_condensation

end module test_water_simple
