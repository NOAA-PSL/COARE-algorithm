module module_coare36_parameters

  use machine, only :  kind_phys &
    ,kind_rad ! for astronomy (date) calculations

  ! air constants and coefficients from the atmospehric model
  use physcons, only: &
    eps =>  con_eps &
    ,cp_a => con_cp &          !< spec heat air @p    (j/kg/k)
    ,epsm1 => con_epsm1 &
    ,sigma_r => con_sbc  &     !< stefan-boltzmann    (w/m2/k4)
    ,omega => con_omega    &    !< ang vel of earth    (1/s)
    ,rvrdm1 => con_fvirt &
    ,rd => con_rd &
    ,rocp => con_rocp  &        !< r/cp
    ,Rgas => con_rd & ! gas constant air ( J/kg/K )
    ,cpa => con_cp & ! spec heat air at p ( J/kg/K )
    ,T2K => con_t0c & ! temp at 0C (K)
    ,rhow => con_rhw0 & ! sea water reference density ( kg/m3 )
    ,pi => con_pi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1  

  public

  real (kind=kind_phys), parameter :: &
    !%***********  set constants **********************************************
    zref = 10., &
    Beta = 1.2, &
    von  = 0.4, &
    fdg  = 1.00, & !% Turbulent Prandtl number
!bloss(Defined above)    T2K  = 273.16, &
    !%***********  air constants **********************************************
!bloss(Defined above)    Rgas = 287.1, &
!bloss(Defined above)    cpa  = 1004.67, &
    !%***********  cool skin constants  ***************************************
    !%%% includes salinity dependent thermal expansion coeff for water
    !%%%%%%%%%%%%%%%%%%%
    bets = 7.5e-4, & !% salintity expansion coeff; assumes beta independent of
                     !  temperature
    !%%%%  see "Computing the seater expansion coefficients directly from the
    !%%%%  1980 equation of state".  J. Lillibridge, J.Atmos.Oceanic.Tech, 1980.
    cpw  = 4000., &
!bloss(Defined above)    rhow = 1022., &
    visw = 1.e-6, &
    tcw  = 0.6, &
    dT_skin_min = 0.001, & ! minimum value for cool-skin temperature
                           ! depression change
                           ! use it for abs(dT_skin-dT_skin_pre) criterion
                           ! for especaping loop
    lw_emiss = 0.97, &  ! longwave emissivity
    sw_onema = 0.945    ! 1-[shortwave albedo]

end module module_coare36_parameters


