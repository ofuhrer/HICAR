module test_lsm_state

    use land_surface, only : initialize_surface_specific_humidity
    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none
    private

    public :: collect_lsm_state_suite

contains

    subroutine collect_lsm_state_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [new_unittest("surface_humidity_cold_start_and_restart", &
                                  test_surface_humidity_cold_start_and_restart)]
    end subroutine collect_lsm_state_suite

    subroutine test_surface_humidity_cold_start_and_restart(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: ims = -1, ime = 2, kms = 0, jms = -2, jme = 1
        real :: water_vapor(ims:ime,kms:2,jms:jme)
        real :: surface_humidity(ims:ime,jms:jme)
        real, parameter :: restart_value = 0.012345
        integer :: i, j, k

        do j = jms, jme
            do k = kms, 2
                do i = ims, ime
                    water_vapor(i,k,j) = 0.001 * real(100*j + 10*k + i + 250)
                enddo
            enddo
        enddo
        surface_humidity = -1.0

        !$acc data copy(water_vapor, surface_humidity)
        call initialize_surface_specific_humidity(surface_humidity, water_vapor, .false., &
                                                  ims, ime, kms, jms, jme)
        !$acc update self(surface_humidity)
        call check(error, all(surface_humidity == water_vapor(:,kms,:)), &
                   "cold start must seed surface humidity from the lowest atmospheric level")
        if (.not. allocated(error)) then
            surface_humidity = restart_value
            !$acc update device(surface_humidity)
            call initialize_surface_specific_humidity(surface_humidity, water_vapor, .true., &
                                                      ims, ime, kms, jms, jme)
            !$acc update self(surface_humidity)
            call check(error, all(surface_humidity == restart_value), &
                       "restart must retain checkpointed surface humidity")
        endif
        !$acc end data
    end subroutine test_surface_humidity_cold_start_and_restart

end module test_lsm_state
