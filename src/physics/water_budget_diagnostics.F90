!> Restart-persistent water-budget accumulation kept outside the Noah-MP
!! driver compilation unit so the diagnostic arithmetic cannot perturb its
!! existing kernels.
module water_budget_diagnostics
    use icar_constants, only : kVARS, kLC_LAND
    use domain_interface, only : domain_t

    implicit none
    private
    public :: accumulate_water_budget

contains

    subroutine accumulate_water_budget( &
        domain, lsm_dt, land_mask, &
        ims, ime, jms, jme, its, ite, jts, jte &
    )
        type(domain_t), intent(inout) :: domain
        real, intent(in) :: lsm_dt
        integer, intent(in) :: ims, ime, jms, jme
        integer, intent(in) :: its, ite, jts, jte
        real, intent(in) :: land_mask(ims:ime, jms:jme)
        integer :: i, j

        ! Noah-MP runoff is already a water depth over the completed soil
        ! timestep. QFX is the signed net surface-water flux rate. Domain
        ! accumulators are restart variables and are never reset by output,
        ! so endpoint differences describe exactly bounded intervals.
        associate(                                                                                         &
            runoff_surface_step => domain%vars_2d(domain%var_indx(kVARS%runoff_surface)%v)%data_2d,       &
            runoff_subsurface_step => domain%vars_2d(domain%var_indx(kVARS%runoff_subsurface)%v)%data_2d, &
            moisture_flux => domain%vars_2d(domain%var_indx(kVARS%QFX)%v)%data_2d,                        &
            runoff_surface_cumulative =>                                                                   &
                domain%vars_2d(domain%var_indx(kVARS%runoff_surface_cumulative)%v)%data_2d,                &
            runoff_subsurface_cumulative =>                                                                &
                domain%vars_2d(domain%var_indx(kVARS%runoff_subsurface_cumulative)%v)%data_2d,             &
            evaporation_net_cumulative =>                                                                  &
                domain%vars_2d(domain%var_indx(kVARS%evaporation_net_cumulative)%v)%data_2d)
        !$acc parallel loop gang vector collapse(2) default(present)
        do j = jts, jte
            do i = its, ite
                ! Noah-MP retains its -9999 fill value where it did not run.
                ! QFX remains valid on active HICAR land, including cells
                ! covered by an external snow model.
                if (land_mask(i,j) == real(kLC_LAND)) then
                    if (runoff_surface_step(i,j) >= 0.0) then
                        runoff_surface_cumulative(i,j) = runoff_surface_cumulative(i,j) + &
                                                        runoff_surface_step(i,j)
                    endif
                    if (runoff_subsurface_step(i,j) >= 0.0) then
                        runoff_subsurface_cumulative(i,j) = runoff_subsurface_cumulative(i,j) + &
                                                           runoff_subsurface_step(i,j)
                    endif
                    evaporation_net_cumulative(i,j) = evaporation_net_cumulative(i,j) + &
                                                     moisture_flux(i,j) * lsm_dt
                endif
            enddo
        enddo
        end associate
    end subroutine accumulate_water_budget

end module water_budget_diagnostics
