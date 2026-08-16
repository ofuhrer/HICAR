!>----------------------------------------------------------
!!  Simple open water flux calculations
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_water_simple
    use data_structures
    use options_interface,   only : options_t
    use mod_wrf_constants,   only : XLV, gravity, KARMAN, cp
    use mod_atm_utilities,   only : sat_mr
    implicit none

contains

    !> Convert the MM5 surface-layer exchange velocity into kinematic and
    !! energetic open-water fluxes.  CHS is already an exchange velocity
    !! [m s-1] (u* kappa / stability denominator), not a dimensionless bulk
    !! coefficient, so wind speed must not be applied a second time.
    pure elemental subroutine simple_water_flux(exchange_velocity, air_density, &
                                                 sst, air_temperature, qv_surf, qv_air, &
                                                 sensible_heat, evap_flux, latent_heat)
        !$acc routine seq
        implicit none
        real, intent(in)  :: exchange_velocity, air_density
        real, intent(in)  :: sst, air_temperature, qv_surf, qv_air
        real, intent(out) :: sensible_heat, evap_flux, latent_heat

        sensible_heat = air_density * cp * (1.0+0.8*qv_air) * exchange_velocity * &
                        (sst-air_temperature)
        evap_flux      = air_density * exchange_velocity * (qv_surf-qv_air)
        latent_heat    = evap_flux * XLV
    end subroutine simple_water_flux

    subroutine water_simple(options, sst, psfc, qv, temperature,  &
                            density, sensible_heat, latent_heat, landmask, lakemask, &
                            qv_surf, evap_flux, tskin, coef_heat_exch, vegtype, ims, ime, kms, kme, jms, jme, its, ite, jts, jte)
        implicit none
        type(options_t),intent(in)    :: options        
        real,    dimension(ims:ime,jms:jme),  intent(inout) :: sensible_heat, latent_heat, qv_surf, evap_flux, tskin
        real,    dimension(ims:ime,jms:jme),  intent(in)    :: sst, psfc, coef_heat_exch, lakemask
        real,    dimension(ims:ime,jms:jme),  intent(in)    :: landmask
        real,    dimension(ims:ime,kms:kme,jms:jme),  intent(in)    :: qv, temperature, density
        integer, dimension(ims:ime,jms:jme), intent(in)    :: vegtype
        integer, intent(in)                     :: ims, ime, kms, kme, jms, jme, its, ite, jts, jte

        integer :: i, j, options_water, options_water_cat, options_lake_cat

        options_water = options%physics%watersurface
        options_water_cat = options%lsm%water_category
        options_lake_cat  = options%lsm%lake_category

        !$acc data present(sst,psfc,qv,temperature,density,sensible_heat,latent_heat,landmask,lakemask,qv_surf,evap_flux,tskin,coef_heat_exch,vegtype) copyin(ims,ime,jms,jme)
        !$acc parallel loop gang vector collapse(2)
        do j=jts,jte
            do i=its,ite
                if(                                                                         &
                    ( (options_water==kWATER_SIMPLE) .AND.                   &   ! If lakemodel is not selected, use this
                      ( (landmask(i,j)==kLC_WATER)                                          & ! ocean / generic water (kLC_WATER landmask), OR
                        .OR. (vegtype(i,j)==options_water_cat)                              & ! cell tagged with the water LU index (in case landmask is inconsistent), OR
                        .OR. ((options_lake_cat>0).AND.(vegtype(i,j)==options_lake_cat))    & ! inland-lake LU cells (typically have landmask==kLC_LAND
                      )                                                                     & !   in USGS/MODIS, so they fall through the landmask gate above)
                    )                                                                       &
                    .OR.                                                                    &
                    ( (options_water==kWATER_LAKE) .AND.                     &    ! if CLM lake model is selected, handle
                      ( (landmask(i,j)==kLC_WATER)                            &    !   any water-tagged cell (by landmask
                        .OR. (vegtype(i,j)==options_water_cat) ) .AND.        &    !   or by LU index) that isn't a lake
                      (lakemask(i,j).ne.1)                                  & !   (lakes go to water_lake)
                    )                                                                         &
                    .OR.                                                                    &
                    ( (options_water==kWATER_FLAKE) .AND.                    &    ! same split for FLake:
                      ( (landmask(i,j)==kLC_WATER)                            &    !   any water-tagged cell (by landmask
                        .OR. (vegtype(i,j)==options_water_cat) ) .AND.        &    !   or by LU index) that isn't a lake
                      (lakemask(i,j).ne.1)                                  & !   (lakes go to flake_step)
                    )                                                                         &
                )then
                    qv_surf(i,j) = 0.98 * sat_mr(sst(i,j),psfc(i,j)) ! multiply by 0.98 to account for salinity
                    ! Floor only - do NOT clamp qv_surf to qv_air; the qv_surf > qv_air
                    ! gradient is exactly what drives evaporation. The old min() clamp
                    ! silently zeroed latent heat over open water in the normal regime.
                    qv_surf(i,j) = max(qv_surf(i,j), 1.0e-6)

                    call simple_water_flux(coef_heat_exch(i,j), density(i,1,j), &
                                           sst(i,j), temperature(i,1,j), qv_surf(i,j), qv(i,1,j), &
                                           sensible_heat(i,j), evap_flux(i,j), latent_heat(i,j))
                    tskin(i,j)   = sst(i,j)

                endif
            end do
        end do
        !$acc end data

    end subroutine water_simple

end module module_water_simple
