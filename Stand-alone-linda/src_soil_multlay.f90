!+ Source module  "src_soil_multlay"
!------------------------------------------------------------------------------

MODULE src_soil_multlay

!------------------------------------------------------------------------------
!
! Description:
!   The module "src_soil_multlay" performs calculations related to the
!   parameterization of soil processes. It contains the soubroutine
!   terra_multlay.f90 which is the combination of the former parts
!   terra1_multlay.incf and terra2_multlay.incf of the LM.
!   All parametric scalar and array data for this soil model routine are
!   defined in the data module data_soil.f90.
!
!   All global variables of the model that are used by the soil model routine
!   terra_multlay.incf is imported by USE statements below.
!
!   The parameterization package has been provided by B. Ritter in a
!   Plug-compatible Fortran77-Version, which is based on the EM/DM soil model
!   by E. Heise. Some technical modifications have been done for the
!   F90 and the parallel Version:
!   Internal communication by common-blocks is replaced by module parameters,
!   scalars and arrays defined in module data_soil. The plug compatible
!   I/O lists of the subroutines have been replaced by the Module interface
!   defined by Use lists below.
!
! Current Code Owner: DWD, Juergen Helmert
!  phone:  +49  69  8062 2704
!  fax:    +49  69  8062 3721
!  email:  Juergen.Helmert@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 2.17       2002/05/08 Reinhold Schrodin
!  Initial release
! 2.18       2002/07/16 Reinhold Schrodin
!  Corrections and further developments
! 3.6        2003/12/11 Reinhold Schrodin
!  Reformulation of bare soil evaporation, transpiration, snow density, 
!  snow height, snow melting, subsoil runoff calculation
! 3.7        2004/02/18 Reinhold Schrodin / Ulrich Schaettler
!  Adapted use of nold to 2 time level scheme and some corrections
! 3.13       2004/12/03 Reinhold Schrodin, Bodo Ritter, Erdmann Heise
!  Snow albedo:   Snow age indicator calculation (freshsnow)
!  Transpiration: 2m temperature and 10m wind considered instead temperature
!                 and wind at lowest atmospheric main level
!  Runoff:        Reformulation of soil water runoff (gravitational part)
!  Soil water:    Soil water calculations restricted to ke_soil_hy,
!                 determined by soil level nearest 2.5 m, combined with
!                 a reformulation of the tridiagonal equation system
!  Snow melting:  Correction in case of pore volume overshooting of
!                 uppermost soil layer by infiltration of snow water
!  Precautions to avoid
!   - soil water content below air dryness point
!   - soil surface temperature problems at lateral LM boundaries.
!   - excessive turbulent fluxes at the soil surface due to mismatch
!     of atmospheric and soil variables in the initial state
!   - evaporation if soil water content is below air dryness point
! 3.16       2005/07/22 Reinhold Schrodin
!   - Modification of evaporation treatment at plant wilting point
!   - Soil water transport: Soil water diffusion term driven by the gradient
!                        of the total soil water concentration (water + ice)
!   - Combining terra1_multlay and terra2_multlay in terra_multlay
! 3.17       2005/12/12 Reinhold Schrodin
!   LME adaptations to GME soil model modifications (B. Ritter):
!   - Constant snow density (250 kg/m**3) replaced by prognostic snow density
!     formulation (min = 50 kg/m**3, max = 400 kg/m**3)
!   - Adjustment of soil temperatures in presence of snow
!   - Extended formulation of transfer coefficient limitation
! 3.18       2006/03/03 Ulrich Schaettler / Erdmann Heise
!  Adaptations to do some initializations also for restart runs
!  Field h_snow is now time dependent, because of use in the FLake model
!  Increase infiltration and reduce surface runoff (bugfix) (Erdmann Heise)
! 3.19       2006/04/25 Erdmann Heise
!  Remove spurious snow and set consistent t_snow and w_snow at the beginning
! 3.21       2006/12/04 Ulrich Schaettler
!  crsmin, rat_lam put to data_soil; eliminated variables that are not used
!  Use dt for ntstep=0 instead of dt/2
!  Additional use of graupel, if present (Thorsten Reinhardt)
!  Modifications to avoid soil water content below air dryness point
!  Gravitational runoff calculation modified in case of frozen soil layer
!                                                      (Reinhold Schrodin)
! V3_23        2007/03/30 Matthias Raschendorfer
!  Introduction of tfv to consider different laminar resistance for heat
!  and water vapour; thus 'rat_lam' is not used here any longer.
! V4_4         2008/07/16 Ulrich Schaettler
!  Splitting of a loop in Section I.4.3b (m_styp is not defined for sea points 
!  and must not occur together with llandmask in the same IF-clause)
! V4_7         2008/12/12 Ulrich Schaettler
!  There were still some loops left with llandmask and m_styp in one line
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Inserted Call to collapse loops
! V4_10        2009/09/11 Christian Bollmann
!  Added compiler directive to use option _on_adb for NEC
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert, Lucio Torrisi
!  Implementation of multi-layer snow model (EM)
!  Use of an external parameter field for stomata resistance (JH)
!  Implementation of ground water as lower boundary of soil column and
!    soil moisture dependent heat conductivity of the soil (JH)
!  Save additional fluxes and stomata resistance to global memory for output (LT)
! V4_12        2010/05/11 Ulrich Schaettler, Ekaterina Machulskaya
!  Renamed t0 to t0_melt because of conflicting names
!  Renamed prs_min to rsmin2d because of conflicting names
!  Update of the new snow model (EM)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler (from H-J Panitz)
!  Introduced snow_melt and ibot_w_so
! V4_18        2011/05/26 Ulrich Schaettler
!  Run the initial steps also for ndfi=1, when nstart > 0!
!  Changed the code owner
! V4_20        2011/08/31 Juergen Helmert
!  Eliminated use of t_2m and use lowest atmospheric layer now (as is in GME)
!  to remove dependency on diagnostic quantity t_2m
! V4_23        2012/05/10 Oliver Fuhrer, Burkhardt Rockel
!  Removed obsolete Fortran features (OF)
!  Correction for multi layer snow model in case of restart (BR)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!                         Ekaterina Machulskaya
!  Replaced qx-variables by using them from the tracer module
!  Added hail rate (in case of two-moment microphysics) to the 
!   precipitation quantities at the ground where it seems necessary.
!  Further developments in the multi-layer snow model (EM)
! V4_26        2012/12/06 Burkhardt Rockel, Ulrich Schaettler
!  Initialize h_snow in case of restart
!  Correct indices for gravity pre-setting (BR)
!  Adapted variable names of multi-layer snow model to corresponding 
!   short names for I/O (US)
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters, ONLY :   &
    ireals,       & ! KIND-type parameter for real variables
    iintegers       ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke_soil,      & ! number of layers in multi-layer soil model
    ke_snow,      & ! number of layers in multi-layer soil model
! LS, b
    dlat,         & ! grid point distance in meridional direction (in degrees)
    dlon,         & ! grid point distance in zonal direction (in degrees)
    degrad,       & ! factor for transforming degree to rad
! LS, e
    czmls,        & ! depth of the main level soil layers in m
                    ! (see organize_physics.f90)

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.
!    
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt              ! long time-step
!   dt2             ! dt*2.

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. physical constants and related variables
! -------------------------------------------

    t0_melt,      & ! absolute zero for temperature
    r_d,          & ! gas constant for dry air    
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d -1  
    cp_d,         & ! specific heat of dry air at constant pressure
    rdocp,        & ! r_d / cp_d
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion         
    lh_s,         & ! latent heat of sublimation    
    g,            & ! acceleration due to gravity
    sigma,        & ! Boltzmann-constant

! 2. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    rho_w,        & ! density of liquid water  
    r_earth         ! mean radius of the earth

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! base state pressure                           (Pa) 

! 2. external parameter fields                                        (unit)
! ----------------------------
    soiltyp    ,    & ! type of the soil (keys 0-9)                     --
    plcov      ,    & ! fraction of plant cover                         --
    rootdp     ,    & ! depth of the roots                            ( m  )
    sai        ,    & ! surface area index                              --
    tai        ,    & ! transpiration area index                        --
    eai        ,    & ! earth area (evaporative surface area) index     --
    llandmask  ,    & ! landpoint mask                                  --
    rsmin2d    ,    & ! minimum stomata resistance                    ( s/m )


! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_snow    ,     & ! temperature of the snow-surface               (  K  )
    t_snow_mult,    & ! temperature of the snow-surface               (  K  )
    tg_radstep,     & ! temperature of the ground surface
                      ! at the last call of radiation routine         (  K  )
    t_s       ,     & ! temperature of the ground surface             (  K  )
    t_g       ,     & ! weighted surface temperature                  (  K  )
    qv_s      ,     & ! specific humidity at the surface              (kg/kg)
    w_snow    ,     & ! water content of snow                         (m H2O)
    rho_snow  ,     & ! snow density                                  (kg/m**3)
    rho_snow_mult,  & ! snow density                                  (kg/m**3)
    h_snow    ,     & ! snow height                                   (  m  
    w_i       ,     & ! water content of interception water           (m H2O)
    t_so      ,     & ! soil temperature (main level)                 (  K  )
    w_so      ,     & ! total water conent (ice + liquid water)       (m H20)
    swi       ,     & ! soil wetness index                            (  1  )
    w_so_ice  ,     & ! ice content                                   (m H20)
! LS, b
    s_so      ,     & ! saturation                                    (1)
    klw       ,     & ! conductivity                                  (m/s)
    dlw       ,     & ! diffusivity                                   (m/s)
    wt_depth  ,     & ! water table depth                             (m)
    satlev    ,     & ! water table depth                             (m)
    q_roff    ,     & ! water table depth                             (m)
    flx       ,     & ! orography                                     (1)
    s_oro     ,     & ! orography                                     (1)
    hsurf     ,     &
! LS, e
    t_2m      ,     & ! temperature in 2m                             (  K  )
    u_10m     ,     & ! zonal wind in 10m                             ( m/s )
    v_10m     ,     & ! meridional wind in 10m                        ( m/s )
    freshsnow ,     & ! indicator for age of snow in top of snow layer(  -  )
    wliq_snow ,     & ! liquid water content in the snow              (m H2O)
    wtot_snow ,     & ! total (liquid + solid) water content of snow  (m H2O)
    dzh_snow          ! layer thickness between half levels in snow   (  m  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   fields of convective and grid-scale precipitation
    prr_con     ,   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con     ,   & ! precipitation rate of snow, convective        (kg/m2*s)
    prr_gsp     ,   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp     ,   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp     ,   & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
    prh_gsp     ,   & ! precipitation rate of hail, grid-scale        (kg/m2*s)

!   fields of the turbulence scheme
    tch         ,   & ! turbulent transfer coefficient for heat       ( -- )
    tcm         ,   & ! turbulent transfer coefficient for momentum   ( -- )
    tfv         ,   & ! laminar reduction factor for evaporation      ( -- )

!   fields of the radiation
    sobs        ,   & ! solar radiation at the ground                 ( W/m2)
    thbs        ,   & ! thermal radiation at the ground               ( W/m2)
    pabs        ,   & ! photosynthetic active radiation               ( W/m2)

! 7. fields for model output and diagnostics                          (unit )
! ---------------------------------------------------------------
    runoff_s    ,   & ! surface water runoff; sum over forecast      (kg/m2)
    runoff_g    ,   & ! soil water runoff; sum over forecast         (kg/m2)
    snow_melt   ,   & ! snow melt; sum over forecast                 (kg/m2)
    rstom       ,   & ! stomata resistance                           ( s/m )
    lhfl_bs     ,   & ! average latent heat flux from bare soil evap.( W/m2)
    lhfl_pl     ,   & ! average latent heat flux from plants         ( W/m2)
! LS, b
    crlat             ! cosine of transformed latitude
! LS, e
! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold  ,       & ! corresponds to ntstep - 1
    nnow  ,       & ! corresponds to ntstep 
    nnew  ,       & ! corresponds to ntstep + 1
    lmelt ,       & ! soil model with melting process
    lmelt_var,    & ! freezing temperature dependent on water content
    ibot_w_so,    & ! number of hydrological active soil layers
    lmulti_snow,  & ! run the multi-layer snow model

! 3. controlling the physics
! --------------------------
    itype_gscp,   & ! type of grid-scale precipitation physics
    itype_trvg,   & ! type of vegetation transpiration parameterization
    itype_evsl,   & ! type of parameterization of bare soil evaporation
    itype_tran,   & ! type of surface to atmospher transfer
    itype_root,   & ! type of root density distribution
    itype_heatcond,&! type of soil heat conductivity
    itype_hydbound,&! type of hydraulic lower boundary condition
    lstomata,     & ! map of minimum stomata resistance

! 5. additional control variables
! --------------------------
     l2tls,       & ! forecast with 2-TL integration scheme
     nextrad ,    &
     hnextrad,    &
     hincrad ,    &
     nincrad
 
! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_io,            ONLY:  &
    lana_rho_snow   ! if .TRUE., take rho_snow-values from analysis file
                    ! else, it is set in the model

!------------------------------------------------------------------------------

USE data_parallel,      ONLY:  &
    num_compute,    & ! number of compute PEs
    icomm_cart,     & ! communicator for the virtual cartesian topology
                      ! that can be used by MPI_WAIT to identify the send
    imp_reals,      & ! determines the correct REAL type used in the model
                      ! for MPI
    ncomm_type,     & ! type of communication
    my_cart_id,     & ! rank of this subdomain in the global communicator
    my_cart_neigh,  & ! neighbors of this subdomain in the cartesian grid
    nboundlines,    & ! number of boundary lines of the domain for which
    sendbuf,        & ! sending buffer for boundary exchange:
                                ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen       ! length of one column of sendbuf


!------------------------------------------------------------------------------

USE src_tracer,      ONLY : trcr_get, trcr_errorstr

USE data_soil       ! All variables from data module "data_soil are used by
                    ! this module. These variables start with letter "c"

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

! External routines used by this module

USE environment,     ONLY : collapse, model_abort, exchg_boundaries
USE meteo_utilities, ONLY : tgcom

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

IMPLICIT NONE


!     Definition of temperature/water related variables in the soil model
!
!
!     --------------T_snow--------
!                                | 
!           zf_snow              |            1 - zf_snow
!                                |      
!     ______________T_s__________|________________T_s_____________
!     ////////////////////////////////////////////////////////////
!
!     Layer k=1 - - - - -  T_SO(k) , W_SO(k), W_SO_ICE(k) - - - -
!
!     ____________________________________________________________
!                                
!                                .
!                                .
!                                . 
!
!     ------------------------------------------------------------
!
!
!     Layer k=ke_soil+1 - - - (climate layer) - - - - - - - - - -
!
!
!     ____________________________________________________________
!
!
!
!     The surface temperature T_g is the snow fraction (zf_snow)
!     weighted sum of T_snow and T_s!
!==============================================================================

CONTAINS

!==============================================================================
!+ Computation of the first part of the soil parameterization scheme
!------------------------------------------------------------------------------

!option! -pvctl _on_adb
SUBROUTINE terra_multlay (yerror, ierror)

!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------

INTEGER (KIND=iintegers)              ::  &
  ierror                        ! error status variable

CHARACTER (LEN=80)                    ::  &
  yerror

CHARACTER (LEN=80)                    ::  yzroutine

!
! Local parameters:
! ----------------

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals , & ! security constant
    zalfa  = 1.0_ireals    , & ! degree of impliciteness (1: full implicit,
                               !    (0.5: Cranck-Nicholson)
    rho_i  = 910._ireals       ! density of solid ice (soil model)  (kg/m**3)

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &

!   Indices

    nx             , & ! time-level for integration
    kso            , & ! loop index for soil moisture layers           
    ksn            , & ! loop index for snow layers
    k              , & ! loop index for snow layers
    kso_zag        , & ! loop index for snow layers
    ke_soil_hy     , & ! number of active soil moisture layers
    i              , & ! loop index in x-direction              
    j              , & ! loop index in y-direction              
    im1, jm1       , & ! i-1, j-1
    jb             , & ! loop index for soil-type               
    mstyp          , & ! soil type index
    msr_off        , & ! number of layers contributing to surface run off
    istarts        , & ! start index for x-direction      
    iends          , & ! end   index for x-direction     
    jstarts        , & ! start index for y-direction     
    jends          , & ! end   index for y-directin      
    zradstep       , & !
! remove that for climate runs
!   i250           , & ! number of layers above a soil depth of 2.50 m
    k10cm          , & ! index of half level closest to 0.1 m
    k100cm             ! index of half level closest to 1.0 m

  REAL    (KIND=ireals   ) ::  &

!   Timestep parameters

    zdt            , & ! integration time-step [s]
    zdelt_bound    , & ! for surface temperature problem at lateral boundaries
    zdel_t_so      , & ! auxiliary variable
    zdtdrhw        , & ! timestep / density of water
    zrhwddt        , & ! density of water / timestep
    zroffdt        , & ! timestep for runoff calculation
    z1d2dt         , & ! 1./2*timestep

!   Connection to the atmosphere

    zuv            , & ! wind velocity in lowest atmospheric layer
    ztvs           , & ! virtual temperature at the surface
    zplow          , & ! pressure of lowest atmospheric layer
    zqvlow         , & ! specific humidity of lowest atmospheric layer
    zrss           , & ! ice covered part of grid element
    zrww           , & ! water covered part of grid element

!   Snow parameters

    zsnow_rate     , & ! rate of snow fall            [kg/m**2 s]
    zdsn_new       , & ! snow age refresh increment   [-]
    zdsn_old       , & ! snow age decay increment     [-]
    zdz_snow(ie,je)    ! snow depth (not for snow heat content but snow
                       ! density calculation)

  REAL    (KIND=ireals   ) ::  &

!   Multi-snow layer parameters
    zsn_porosity          , & ! snow porosity                    (   -   )
    zp1                   , &
    zfukt                 , &
    zq0                   , &
    zqbase(ie,je)         , &
    zrefr(ie,je)          , & ! rate of liquid water refreezing
    zmelt(ie,je)          , & ! rate of snow melting
    ze_in                 , &
    ze_out(ie,je)         , &
    zflag                 , &
    zadd_dz               , &
    zrho_dry_old(ie,je)   , &
    zeta                  , &
    zdens_old             , &
    zcounter(ie,je)       , &
    ze_rad(ie,je)         , &
    zswitch(ie,je)        , &
    time1                 , &
    time2                 , &
    sum_weight(ie,je)     , &
    zp     (ie,je,ke_snow), &
    t_new  (ie,je,ke_snow), &
    rho_new(ie,je,ke_snow), &
    wl_new (ie,je,ke_snow), &
    z_old  (ie,je,ke_snow), &
    dz_old (ie,je,ke_snow), &
    weight                , &
    thbs_new(ie,je)       , &
    tmp_num               , &
    tmp_the               , &

!   Plant parameters
    zbeta                 , & ! reduction factor for evaporation
    zevap                     ! auxiliary variable

  REAL    (KIND=ireals   ) ::  &

!   Local scalars for BATS-scheme

    zpsi0=.01_ireals,& ! air entry potential at water saturation (m)
    zbf1           , & ! auxiliary variable
    zbf2           , & ! auxiliary variable
    zdmax          , & ! auxiliary variable
    zd             , & ! auxiliary variable
    zck            , & ! auxiliary variable
    z1             , & ! depth of half level closest to 0.1 m
    znull          , & ! depth of half level closest to 1.0 m
    zfqmax         , & ! maximum sustainable flux in uppermost soil layer
    zevapor        , & ! evaporation rate of bare soil in BATS-scheme
    zrla           , & ! atmospheric resistance
    zcatm          , & ! transfer function CA
    zpar           , & ! PAR (interim value)
    zf_wat         , & ! soil water function for stomatal resistance
    zf_tem         , & ! temperature function for stomatal resistance
    zf_sat         , & ! saturation deficit function for stomatal resistance
    zepsat         , & ! saturation vapour pressure at near surface temperature
    zepke          , & ! near surface water vapour pressure
    zrstom         , & ! stomata resistance
    zrveg          , & ! sum of atmospheric and stomata resistance
    zedrstom       , & ! 1./zrstom
    zustar         , & ! friction velocity
    ztrabpf        , & ! area average transpiration rate
    ! variables for non-uniform root-distribution
    ztrfr          , & ! fraction of total transpiration associated with each layer
    zrootfr        , & ! normalized root fraction in each layer (=1 for z_soil=0.0)
    zrootdz        , & ! normalized root fraction * min(layer thickness,rootlength-top_of_layer)

!   Soil parameters

    zice           , & ! indicator of soil type ice
    zwqg           , & ! mean of fcap and pwp
    z4wdpv         , & ! 4*zwqg/porv         
    zroot          , & ! fraction of layer filled by roots (Bucket scheme)
    zropart        , & ! fraction of layer filled by roots (Dickinson scheme)
    zrootfc        , & ! distributes transpiration to soil layers
    zr_root        , & ! zroot/zbwt (Bucket scheme)
    ze_sum             ! sum of all contributions to evapotranspiration

  REAL    (KIND=ireals   ) ::  &

!   Thermal parameters

    ztgt0          , & ! Indicator T_g > T_0
    zgstr          , & ! downward longwave radiation
    zrnet_s(ie,je)        , & ! net radiation
    zshfl_s(ie,je)        , & ! sensible heatflux at soil surface
    zlhfl_s(ie,je)        , & ! latent heatflux at soil surface
    zsprs(ie,je)          , & ! utility variable
    zalas          , & ! heat conductivity of snow
    zrnet_snow     , & ! net radiation at snow surface
    zfak           , & ! utility variable for implicit snow temperature forecast
    ztsnow_im      , & ! utility variable for implicit snow temperature forecast
    ztsnew         , & ! preliminary value of snow surface temperature
    ztsnownew      , & ! preliminary value of snow surface temperature
    zshfl_snow     , & ! sensible heatflux at snow surface
    zlhfl_snow     , & ! latent heatflux at snow surface
    zfor_snow      , & ! total forcing at snow surface
    zfr_melt       , & ! melting snow fraction
!    zdwsnm         , & ! utility variable for snow melt determination
    zdwgme         , & ! utility variable for snow melt determination
    zdelt_s        , & ! utility variable for snow melt determination
    ze_avail       , & ! utility variable for snow melt determination
    ze_total           ! utility variable for snow melt determination

  ! Some of these values have to be fields for the multi-layer snow model
  REAL    (KIND=ireals   ) ::  &
    zalas_mult    (ie,je,ke_snow),    & ! heat conductivity of snow
    ztsnownew_mult(ie,je,0:ke_snow),  & ! preliminary value of snow surface temperature
    zextinct      (ie,je,ke_snow),    & ! solar radiation extinction coefficient in snow (1/m)
    zfor_snow_mult(ie,je),            & ! total forcing at snow surface
    zfor_total    (ie,je)               ! total forcing at surface
  INTEGER (KIND=iintegers   ) ::  &
    zicount1      (ie,je),            &
    zicount2      (ie,je)

  REAL    (KIND=ireals   ) ::  &

!   Hydraulic parameters

    zwinstr        , & ! preliminary value of interception store
    zinfmx         , & ! maximum infiltration rate
    zwimax         , & ! maximum interception store
    zalf           , & ! utility variable
    zvers          , & ! water supply for infiltration
    zro_inf        , & ! surface runoff
    zdwsndtt       , & ! time rate of change of snow store
    zwsnstr        , & ! preliminary value of snow store
    zdwseps        , & ! artificial change of small amount of snow store
    zdwidtt        , & ! time rate of change of interception store
    zdwieps        , & ! artificial change of small amount of interception store
    zro_wi         , & ! surface runoff due to limited infiltration capacity
    zro_sfak       , & ! utility variable 
    zro_gfak       , & ! utility variable
    zfmb_fak       , & ! utility variable
    zdwg           , & ! preliminary change of soil water content
    zwgn           , & ! preliminary change of soil water content
    zredfu         , & ! utility variable for runoff determination
    zro            , & ! utility variable for runoff determination
    zro2           , & ! utility variable for runoff determination
    zkorr          , & ! utility variable for runoff determination
    zfr_ice        , & ! reduction factor for water transport
    zfr_ice_free   , & ! reduction factor for water transport
    zwso_new       , & ! preliminary value of soil water content
    w_p            , & ! preliminary value of soil water content
!    zwsnew         , & ! preliminary value of snow water equivalent
    zw_ovpv        , & ! utility variable

!   Implicit solution of thermal and hydraulic equations

    zakb           , & ! utility variable
    zakb_m1        , & ! utility variable
    zakb_p1        , & ! utility variable
    zzz            , & ! utility variable
    z1dgam1        , & ! utility variable
    zredm          , & ! utility variable
    zredm05        , & ! utility variable
    zredp05        , & ! utility variable
    zred           , & ! utility variable
    zgam2m05       , & ! utility variable
    zgam2p05           ! utility variable

  REAL    (KIND=ireals   ) ::  &

!   Statement functions

    zsf_heav       , & ! Statement function: Heaviside function
    zsf_psat_iw    , & ! Saturation water vapour pressure over ice or water
                       ! depending on temperature "zstx"
    zsf_qsat       , & ! Specific humidity at saturation pressure
                       ! (depending on the saturation water vapour pressure
                       !  "zspsatx" and the air pressure "zspx")
    zsf_dqvdt_iw   , & ! Statement function: First derivative of specific
                       ! saturation humidity
                       ! with respect to temperature (depending on temperature
                       ! "zstx" and saturation specific humidity pressure
                       ! "zsqsatx")
    zstx           , & ! dummy argument for Stmt. function
    zspx           , & ! dummy argument for Stmt. function
    zspsatx        , & ! dummy argument for Stmt. function
    zsqsatx        , & ! dummy argument for Stmt. function
    z2iw           , & ! dummy argument for Stmt. function
    z4iw           , & ! dummy argument for Stmt. function
    z234iw         , & ! dummy argument for Stmt. function
    zqs            , & ! saturation humidty at T_s and ps
    zdqs           , & ! derivative of zqs with respect to T_s
    zqsnow         , & ! saturation humidty at T_snow and ps
    zdqsnow        , & ! derivative of zqs with respect to T_snow

!   Local scalars for transfer coefficient limitation calculations and for
!   treatment of possible problems at lateral boundaries

    zch_soil = 1.400E06_ireals ,& ! approximate average heat capacity of soil (J/m**3 K)
    zlim_dtdt = 2.5_ireals   ,&! maximum allowed temperature increment (K) in
                               ! one time step in  uppermost soil layer

   zdT_snow                            ,& ! maximum possible increment in snow layer before
                                          ! melting sets in
   zg1       (ie         ,je         ) ,& ! heat flux between two uppermost soil layers (W/m**2)
   zlhfl     (ie         ,je         ) ,& ! estimated       latent    heat flux surface (W/m**2)
   zshfl     (ie         ,je         ) ,& ! estimated       sensible  heat flux surface (W/m**2)
   zthfl     (ie         ,je         ) ,& ! estimated total turbulent heat flux surface (W/m**2)
   zradfl    (ie         ,je         ) ,& ! total radiative flux at surface             (W/m**2)
   ze_melt   (ie         ,je         ) ,& ! energy required for melting of snow         (J/m**2)
   zch_snow  (ie         ,je         ) ,& ! heat capacity of snow * w_snow * Rho_w      (J/m**2 K)
   zeb1      (ie         ,je         ) ,& ! estimated energy budget of first soil layer (W/m**2)
   ztchv     (ie         ,je         ) ,& ! transfer coefficient multiplied by wind velocity
   ztchv_max (ie         ,je         ) ,& ! dito, but upper limit
   zrho_atm  (ie         ,je         ) ,& ! air density of lowest atmospheric layer     (kg/m**3)
   zdt_atm   (ie         ,je         ) ,& ! temperature difference between lowest
   zdq_atm   (ie         ,je         ) ,& ! dito for moisture (kg/kg)
   !additional variables for root distribution
   zroota    (ie         ,je         ) ,& ! root density profile parameter (1/m)
   zwrootdz  (ie         ,je, ke_soil) ,& ! mean water content over root depth weighted by root density
   zrootdz_int (ie       ,je         ) ,& ! parameter needed to initialize the root density profile integral
   zwrootdz_int(ie       ,je         )    ! parameter needed to initialize the root water content integral


   INTEGER m_limit                          ! counter for application of limitation
   CHARACTER (LEN=7) yhc                    ! heating or cooling indicator


  REAL    (KIND=ireals   ) ::  &

!   Local scalars for water content dependent freezing/melting

    zliquid        , & ! utility variable
    zaa            , & ! utility variable
    znen           , & ! utility variable
    zargu          , & ! utility variable

!   Water transport

    ! LINDA, b
!    zice_fr_ksom1  , & ! fractional ice content of actual layer - 1
!    zice_fr_kso    , & ! fractional ice content of actual layer
!    zice_fr_ksop1  , & ! fractional ice content of actual layer + 1
!    zlw_fr_ksom1   , & ! fractional liquid water content of actual layer - 1
!    zlw_fr_ksom1_new,& ! fractional liquid water content of actual layer -1
!    zlw_fr_kso     , & ! fractional liquid water content of actual layer
!    zlw_fr_kso_new , & ! fractional liquid water content of actual layer
!    zlw_fr_ksop1   , & ! fractional liquid water content of actual layer + 1
    zlw_fr_ksom05  , & ! fractional liquid water content of actual layer - 1/2
!    zdlw_fr_ksom05 , & ! hydraulic diffusivity coefficient at half level above
!    zklw_fr_ksom05 , & ! hydraulic conductivity coefficient at half level above
!    zklw_fr_kso_new, & ! hydraulic conductivity coefficient at main level
                       !    for actual runoff_g
    zlw_fr_ksop05  , & ! fractional liquid water content of actual layer + 1/2
!    zdlw_fr_ksop05 , & ! hydraulic diffusivity coefficient at half level below
!    zklw_fr_ksop05 , & ! hydraulic conductivity coefficient at half level below

    zinf           , & ! infiltration 
    zlw_fr_kso(ie,je,1:ke_soil+1)     , &  ! fractional liquid water content of actual layer
    zice_fr_kso(ie,je,1:ke_soil+1)    , &  ! fractional liquid water content of actual layer
    zdlw_fr_ksop05(ie,je,1:ke_soil+1) , &  ! hydraulic diffusivity coefficient at half level below
    zklw_fr_ksop05(ie,je,1:ke_soil+1) , &  ! hydraulic conductivity coefficient at half level below
    zdlw_fr_ksop05_fg(ie,je,1:ke_soil+1),& ! first guess hydraulic diffusivity coefficient at half level below
    zklw_fr_ksop05_fg(ie,je,1:ke_soil+1),& ! first guess hydraulic conductivity coefficient at half level below
    zklw_fr_ksop05_sg(ie,je,1:ke_soil+1),& ! second guess hydraulic conductivity coefficient at half level below
    q_kso(ie,je,1:ke_soil+1)            ,& ! ground water flow

    ! LINDA, e
!   Snow density
    ztau_snow      , & ! 'ageing constant' for snow density          (-)
    zrho_snowe     , & ! updated density of existing snow ('ageing') kg/m**3
    zrho_snowf     , & ! density of fresh snow                       kg/m**3
    znorm          , & ! normalisation factor for weighted snow density mH2O


!   Freezing/melting of soil water/ice:

    zenergy        , & ! available melting/freezing energy per time step
    zdelwice       , & ! amount of melted soil ice/frozen soil water
    zdwi_max       , & ! maximum amount of freezing/melting soil water
    ztx                ! water content dependent freezing/melting temperature

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND = ireals) :: &

! Model geometry

  zmls     (ke_soil+1)  , & ! depth of soil main level
  zzhls    (ke_soil+1)  , & ! depth of the half level soil layers in m
  zdzhs    (ke_soil+1)  , & ! layer thickness between half levels
  zdzms    (ke_soil+1)  , & ! distance between main levels
  zdz_snow_fl(ie,je)    , & ! snow depth for snow temperature gradient

  ! Multi-layer snow model
  zhh_snow (ie,je,ke_snow)    , & ! depth of the half level snow layers
  zhm_snow (ie,je,ke_snow)    , & ! depth of snow main levels
  zdzh_snow(ie,je,ke_snow)    , & ! layer thickness between half levels
  zdzm_snow(ie,je,ke_snow)    , & ! distance between main levels

! External parameters

  zbwt     (ie,je)      , & ! root depth (with artificial minimum value)
  zrock    (ie,je)      , & ! ice/rock-indicator: 0 for ice and rock
  zsandf   (ie,je)      , & ! mean fraction of sand (weight percent)
  zclayf   (ie,je)      , & ! mean fraction of clay (weight percent)
  zb_por   (ie,je)      , & ! pore size distribution index
  zpsis    (ie,je)      , & ! air entry potential (m)
  zw_m     (ie,je)          ! maximum of liquid water content  (m)

  INTEGER  (KIND=iintegers ) ::  &
  m_styp   (ie,je)          ! soil type

  REAL    (KIND=ireals   ) ::  &

! Connection to the atmosphere

  zrr      (ie,je)      , & ! total rain rate including formation of dew
  zrs      (ie,je)      , & ! total snow rate including formation of rime
  zesoil   (ie,je)      , & ! evaporation from bare soil
  zrhoch   (ie,je)      , & ! transfer coefficient*rho*g
  zf_snow  (ie,je)      , & ! surface fraction covered by snow
  zth_low  (ie,je)      , & ! potential temperature of lowest layer
  zf_wi    (ie,je)      , & ! surface fraction covered by interception water
  ztmch    (ie,je)      , & ! heat transfer coefficient*density*velocity
  zep_s    (ie,je)      , & ! potential evaporation for t_s
  zep_snow (ie,je)      , & ! potential evaporation for t_snow
  zverbo   (ie,je)      , & ! total evapotranspiration
  zversn   (ie,je)      , & ! total evaporation from snow surface
  zthsoi   (ie,je)      , & ! thermal flux at soil surface
  zthsnw   (ie,je)      , & ! thermal flux at snow surface
  zfor_s   (ie,je)      , & ! total forcing at soil surface
  zgsb     (ie,je)      , & ! heat-flux through snow

! Tendencies

  zdwidt   (ie,je)      , & ! interception store tendency
  zdwsndt  (ie,je)      , & ! snow water content tendency
  zdtsdt   (ie,je)      , & ! tendency of zts
  zdtsnowdt(ie,je)      , & ! tendency of ztsnow
  zdtsnowdt_mult(ie,je,0:ke_snow)      , & ! tendency of ztsnow
  zdwgdt   (ie,je,ke_soil)  ! tendency of water content [kg/(m**3 s)]

  REAL    (KIND=ireals   ) ::  &

! Soil and plant parameters

  zroc(ie,je,ke_soil+1) , & ! heat capacity
  zfcap    (ie,je)      , & ! field capacity of soil
  zadp     (ie,je)      , & ! air dryness point
  zporv    (ie,je)      , & ! pore volume (fraction of volume)
  zdlam    (ie,je)      , & ! heat conductivity parameter
  zdw      (ie,je)      , & ! hydrological diff.parameter
  zdw1     (ie,je)      , & ! hydrological diff.parameter
  ! LINDA, b
  !zkw      (ie,je)      , & ! hydrological cond.parameter
  zkw      (ie,je,ke_soil+1), & ! hydrological cond.parameter
  ! LINDA, e
  zkw1     (ie,je)      , & ! hydrological cond.parameter
  zik2     (ie,je)      , & ! minimum infiltration rate
  zpwp     (ie,je)      , & ! plant wilting point  (fraction of volume)
  ztlpmwp  (ie,je)      , & ! turgor-loss-point minus plant wilting point
  zedb     (ie,je)      , & ! utility variable

! Hydraulic variables

  ztrang(ie,je,ke_soil) , & ! transpiration contribution by the different layers
  ztrangs(ie,je)        , & ! total transpiration (transpiration from all
                            !    soil layers)
  zwin     (ie,je)      , & ! water cont. of interception store   (m H20)
  zwsnow   (ie,je)      , & ! snow water equivalent               (m H20)
  zwsnew   (ie,je)      , & ! snow water equivalent               (m H20)
  zdwsnm   (ie,je)      , & ! utility variable for snow melt determination
  zw_fr    (ie,je,ke_soil+1),&!fractional total water content of soil layers
  zinfil   (ie,je)      , & ! infiltration rate
  zlw_fr   (ie,je,ke_soil+1), & ! fractional liqu. water content of soil layer
  ziw_fr   (ie,je,ke_soil+1), & ! fractional ice content of soil layer
  zwsnn    (ie,je)          , & ! new value of zwsnow
  zflmg    (ie,je,ke_soil+1), & ! flux of water at soil layer interfaces
  zrunoff_grav(ie,je,ke_soil+1),&! main level water gravitation
                                !    flux (for runoff calc.)
  ! LINDA, b
  zsoro    (ie,je)          , & ! orographie parameter for groundwater runoff
  zgamma   (ie,je)
  ! LINDA, e

  REAL    (KIND=ireals   ) ::  &

! Local arrays for BATS-scheme

  zk0di    (ie,je)      , & ! surface type dependent parameter
  zbedi    (ie,je)      , & ! surface type dependent parameter
  zsnull   (ie,je)      , & ! mean relative(zporv) water fraction of the
                            !    so-called total active soil layer (above 1m)
  zs1      (ie,je)      , & ! mean relative (zporv) water fraction  of
                            !    layers above 0.1m
  zf_rad   (ie,je)      , & ! radiation function for stomata resistance
  ztraleav (ie,je)      , & ! transpiration rate of dry leaves

! Root functions

  zwroot   (ie,je)      , & ! mean water content over root depth
  zropartw(ie,je,ke_soil)   ! fraction of layer filled by roots * w_fr

  REAL    (KIND=ireals   ) ::  &

! Thermal variables

  zts      (ie,je)      , & ! soil surface temperature
  ztsnow   (ie,je)      , & ! snow surface temperaure
  ztsnow_mult   (ie,je,0:ke_snow)      , & ! snow surface temperaure
! HEATCOND (soil moisture dependent heat conductivity of the soil)
  zalamtmp (ie,je)      , & ! heat conductivity
  zalam    (ie,je,ke_soil),&! heat conductivity
  zrocg    (ie,je)      , & ! volumetric heat capacity of bare soil
  zrocs    (ie,je)      , & ! heat capacity of snow
  ztsn     (ie,je)      , & ! new value of zts
  ztsnown  (ie,je)      , & ! new value of ztsnow
  ztsnown_mult  (ie,je,0:ke_snow)      , & ! new value of ztsnow
  znlgw1f  (ke_soil)    , & ! utility variable
  zaga(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zagb(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zagc(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zagd(ie,je,0:ke_soil+ke_snow+1),& ! utility variable
  zage(ie,je,0:ke_soil+ke_snow+1),& ! utility variable

! Auxiliary variables

  zw_snow_old(ie,je)    , & !
  zrho_snow_old(ie,je)    , & !
  zdqvtsnow(ie,je)      , & ! first derivative of saturation specific humidity
                            !    with respect to t_snow
  zts_pm   (ie,je)      , & ! indicator zts > < T_melt
  ztsnow_pm(ie,je)          ! indicator ztsnow > < T_melt

  LOGICAL     :: &

  ldebug                , & !
  !LINDA,b
  ldecharme             , & ! Decharme formulation
  limit                 , & !
  !LINDA,e
  limit_tch (ie,je)         ! indicator for flux limitation problem

  ! ground water as lower boundary of soil column
  REAL    (KIND=ireals) ::  &
    zdelta_sm, zdhydcond_dlwfr, zklw_fr_kso, zklw_fr_ksom1, zdlw_fr_kso 
  
  ! HEATCOND
  REAL    (KIND=ireals) ::          &
    hzalam   (ie,je,ke_soil+1),     & ! heat conductivity
    zthetas, zlamli, zlamsat, zlams, rsandf, zlamq, zlam0, zrhod, zlamdry,  &
    zsri, zKe, zthliq, zlamic

! LINDA, first saturated level
   integer  zsatlev(ie,je)
   real  z_wt(ie,je), h_wt(ie,je)
   real (KIND=ireals   )     zklw_fr_ksom05_limit
   INTEGER(KIND=iintegers) :: kzdims(24)
 
! Tracer pointers
!----------------
  REAL (KIND=ireals), POINTER :: &
    qv  (:,:,:) => NULL()                    ! QV at tlev=nx

! LS, b
! Additional fields for soil-moisture relaxation                  (unit )
! ---------------------------------------------------------------

  real (KIND=ireals   )     flux_ksom05,flux_ksom05_pos, flux_ksop05, flux_ksop05_neg, q_kso_pos, alpha_lim
  real (KIND=ireals   )     dx, dy, len, oolen, oolenx, ooleny, grad(9)
! LS, e

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine terra_multlay              
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the diagnostic part I of the soil parameterization scheme
!  In this part, evaporation from the surface and transpiration is calculated.
!  A multi-layer soil water distribution is calculated by simple bulk
!  parameterisation or a Penman-Monteith-type version of evaporation
!  and transpiration which can be used alternatively.
!------------------------------------------------------------------------------

! Declaration of STATEMENT-FUNCTIONS

  zsf_heav     (zstx                    ) = 0.5_ireals+SIGN( 0.5_ireals, zstx )
  zsf_psat_iw  (zstx,z2iw   ,z4iw       )                                     &
                   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))
  zsf_qsat     (zspsatx, zspx           )                                     &
                   = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)
  zsf_dqvdt_iw (zstx,zsqsatx,z4iw,z234iw)                                     &
                   = z234iw*(1._ireals+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2

!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

  ierror    = 0
  yerror    = '        '
  yzroutine = 'terra_multlay'

  ldecharme = .true.
  
! select timelevel and timestep for calculations
! In the new formulation a 2 timelevel scheme is used for the soil model.
! Therefore it is not necessary any more to distinguish between Leapfrog
! and Runge-Kutta scheme.
! But the first timestep must not be done with dt/2, as in the Leapfrog
! scheme
  nx  = nnow
  IF ( (ntstep == 0) .AND. (.NOT. l2tls) ) THEN
    ! use the original dt and not dt/2
    zdt = 2 * dt
  ELSE
    zdt = dt
  ENDIF

  ! time step for run-off computation
  zroffdt = zdt

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar
  jstarts = jstartpar
  jends   = jendpar

#ifdef NECSX
  CALL collapse(.TRUE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif

  ! retrieve of the required microphysics tracers (at timelevel nx)
  CALL trcr_get(ierror, 'QV', ptr_tlev = nx, ptr = qv)
  IF (ierror /= 0) THEN
    yerror = trcr_errorstr(ierror)
    CALL model_abort(my_cart_id, ierror, yerror, yzroutine)
  ENDIF

! Computation of derived constants

  zrhwddt = rho_w/zdt     ! density of liquid water/timestep
  zdtdrhw = zdt/rho_w     ! timestep/density of liquid water

! time constant for infiltration of water from interception store
  ctau_i        = MAX(ctau_i,zdt)

! grids for temperature and water content

  zzhls(1) = 2._ireals*czmls(1)   !depth of first half level
  zdzhs(1) = zzhls(1)      !layer thickness betw. half levels of uppermost layer
  zmls(1)  = czmls(1)      !depth of 1st main level
  zdzms(1) = czmls(1)      !layer thickness between soil surface and main level
                           ! of uppermost layer
 
  DO kso = 2,ke_soil+1
    zzhls(kso)  = zzhls(kso-1) + 2._ireals*(czmls(kso) -zzhls(kso-1))
    zdzhs(kso) = zzhls(kso) - zzhls(kso-1) ! layer thickness betw. half levels
    zmls(kso)  = czmls(kso)                ! depth of main levels
    zdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
  ENDDO


! Prepare basic surface properties (for land-points only)
 
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF(llandmask(i,j)) THEN        ! for land-points only
        mstyp       = NINT(soiltyp(i,j))        ! soil type
        m_styp(i,j) = mstyp                     ! array for soil type
        zdw   (i,j)  = cdw0  (mstyp)
        zdw1  (i,j)  = cdw1  (mstyp)
        ! LINDA, b
        !zkw   (i,j)  = ckw0  (mstyp)
        zkw   (i,j,ke_soil+1) = ckw0  (mstyp)   ! other levels are set below with 3D variables
        ! LINDA, e
        zkw1  (i,j)  = ckw1  (mstyp)
        zik2  (i,j)  = cik2  (mstyp)
        zporv(i,j)  = cporv(mstyp)              ! pore volume
        zpwp (i,j)  = cpwp (mstyp)              ! plant wilting point
        zadp (i,j)  = cadp (mstyp)              ! air dryness point
        zfcap(i,j)  = cfcap(mstyp)              ! field capacity
        zrock(i,j)  = crock(mstyp)              ! rock or ice indicator
        zrocg(i,j)  = crhoc(mstyp)              ! heat capacity
        zdlam(i,j)  = cala1(mstyp)-cala0(mstyp) ! heat conductivity parameter
        zbwt(i,j)   = MAX(0.001_ireals,rootdp(i,j))! Artificial minimum value
                                                ! for root depth
        zroota(i,j) = 3._ireals/zbwt(i,j)       ! root density profile parameter (1/m)
                                                ! zroota=0. creates the original TERRA_LM
                                                ! version with constant root density
                                                ! for root depth
        ! New arrays for BATS-scheme
        zk0di(i,j)  = ck0di(mstyp)              !
        zbedi(i,j)  = cbedi(mstyp)              !
        ! Arrays for soil water freezing/melting
        zsandf(i,j)   = csandf(mstyp)
        zclayf(i,j)   = cclayf(mstyp)
        zpsis(i,j)    = -zpsi0 * 10**(1.88_ireals-0.013_ireals*zsandf(i,j))
        zb_por(i,j)   = 2.91_ireals + .159_ireals*zclayf(i,j)
        zedb(i,j)     = 1._ireals/zb_por(i,j)
        !zsoro(i,j)    = 0.05_ireals
        zgamma(i,j)   = 10.0_ireals
      ENDIF
    ENDDO
  ENDDO

  ! Set three-dimensional variables
  DO kso = 1, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF(llandmask(i,j)) THEN        ! for land-points only
          mstyp           = NINT(soiltyp(i,j))        ! soil type
          zalam(i,j,kso)  = cala0(mstyp)              ! heat conductivity parameter
          IF (ldecharme) then
!<JH
    !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
             zkw   (i,j,kso) = ckw0 (mstyp)*EXP(-2._ireals*(zmls(kso)-rootdp(i,j)))
!>JH
          else
             zkw   (i,j,kso) = ckw0 (mstyp)
          end if
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  DO kso = 1, ke_soil+1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF(llandmask(i,j)) THEN        ! for land-points only
          q_kso(i,j,kso) = 0.0_ireals
          zklw_fr_ksop05_fg(i,j,kso) = 0.0_ireals
          zklw_fr_ksop05_sg(i,j,kso) = 0.0_ireals
          zklw_fr_ksop05(i,j,kso) = 0.0_ireals
          zdlw_fr_ksop05_fg(i,j,kso) = 0.0_ireals
          zdlw_fr_ksop05(i,j,kso) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDDO

! For ntstep=nstart : Some preparations
! =====================================

  IF (ntstep == nstart) THEN
!   Determine constants clgk0 for BATS-scheme
    DO jb       = 1, 10
      clgk0(jb) = LOG10(MAX(zepsi,ck0di(jb)/ckrdi))
    END DO
  ENDIF

! For ntstep=0 : Some preparations
! ================================

  ! it has to be checked, if it is the first step after a forward
  ! launching digital filtering initialization: then lsoilinit_dfi
  ! is set to TRUE in dfi_initialization
  IF ((ntstep == 0) .OR. lsoilinit_dfi) THEN
    lsoilinit_dfi = .FALSE.

    w_so(:,:,:,nnew) = w_so(:,:,:,nx)

!   Provide for a soil moisture 1 % above air dryness point, reset soil
!   moisture to zero in case of ice and rock
    DO kso   = 1,ke_soil+1
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN             ! for land-points only
            IF (m_styp(i,j).ge.3) THEN
              w_so (i,j,kso,:) = MAX(w_so(i,j,kso,:),                     &
                                           1.01_ireals*zadp(i,j)*zdzhs(kso) )
            ELSE
              w_so(i,j,kso,:) = 0.0_ireals
            ENDIF

          END IF   ! land-points
        END DO
      END DO
    END DO

!   adjust temperature profile in lower soil layers, if temperature of first soil
!   layer was reduced due to the detection of snow (e.g. in the analysis)
!   loop over grid points
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF(llandmask(i,j)) THEN   ! for land-points only
          IF (w_snow(i,j,nx) <=1.0E-6_ireals) THEN
            ! spurious snow is removed
            t_snow(i,j,:)=t_so(i,j,0,nx)
            w_snow(i,j,nx)= 0.0_ireals
            IF(lmulti_snow) THEN
              t_snow_mult(i,j,0,nx) = t_so(i,j,0,nx)
              DO ksn = 1, ke_snow
                t_snow_mult  (i,j,ksn,nx) = t_so(i,j,0,nx)
                wliq_snow    (i,j,ksn,nx) = 0.0_ireals
                wtot_snow    (i,j,ksn,nx) = 0.0_ireals
                rho_snow_mult(i,j,ksn,nx) = 0.0_ireals
                dzh_snow     (i,j,ksn,nx) = 0.0_ireals
              END DO
            END IF
          ELSE
!           adjust soil temperatures in presence of snow
            zdel_t_so = MAX (0._ireals,t_so(i,j,0,nx)-t0_melt)
            t_so(i,j,0,:)=MIN( t0_melt,t_so(i,j,0,nx) )
            DO kso=1,ke_soil
              IF ( t_so(i,j,kso,nx) > t0_melt) THEN
                 t_so(i,j,kso,:) = MAX( t0_melt,                  &
                                   t_so(i,j,kso,nx) -zdel_t_so    &
                                   *(zmls(ke_soil+1)-zmls(kso))   &
                                   /(zmls(ke_soil+1)-zmls( 1 )) )
              ENDIF
              IF ( t_so(i,j,kso,nx) <= t0_melt) EXIT
            ENDDO
          ENDIF

          t_s(i,j,:) = t_so(i,j,0,:)

!         Set level 1 to level 0 for t_so for every landpoint
          t_so(i,j,1,:) = t_so(i,j,0,nx)
        ENDIF    ! llandmask
      END DO
    END DO


!   Initialization of soil ice content
!   ----------------------------------
!   Generally the combination lmelt=.true., lemlt_var=.true. has to be used.
!   lmelt_var=.false. (soil water freezing at 273.15K) should be used only
!   for special investigations because at model start a partial frozen layer
!   is only considered if the soil ice water content is read in from
!   a pre run.
!   ----------------------------------
    IF (lmelt .AND. .NOT. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
        DO   j = jstarts, jends
          DO i = istarts, iends
            IF (llandmask(i,j)) THEN             ! for land-points only
              w_so_ice(i,j,kso,:) = 0.0_ireals
              IF (t_so(i,j,kso,nx) < t0_melt) THEN
                w_so_ice(i,j,kso,:) = w_so(i,j,kso,nx)
              END IF ! t_so(kso) < t0_melt
            END IF   ! land-points
          END DO
        END DO
      END DO
    END IF           ! lmelt .AND. .NOT. lmelt_var
    IF(lmelt .AND. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
        DO   j = jstarts, jends
          DO i = istarts, iends
            IF (llandmask(i,j)) THEN             ! for land-points only
              IF (t_so(i,j,kso,nx) < (t0_melt-zepsi)) THEN
                zaa    = g*zpsis(i,j)/lh_f
                zw_m(i,j)     = zporv(i,j)*zdzhs(kso)
                zw_m(i,j)   = zw_m(i,j)*                                          &
                  ((t_so(i,j,kso,nx) - t0_melt)/(t_so(i,j,kso,nx)*zaa))**(-zedb(i,j))
                w_so_ice(i,j,kso,:) = MAX (0.0_ireals,w_so(i,j,kso,nx) - zw_m(i,j))
              END IF
            END IF   ! land-points
          END DO
        END DO
      END DO
    END IF           ! lmelt .AND. lmelt_var

!   Initialization of snow density, if necessary
!   --------------------------------------------
    IF (.NOT. lana_rho_snow) THEN
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN             ! for land-points only
            IF(lmulti_snow) THEN
              DO ksn = 1, ke_snow
                rho_snow_mult(i,j,ksn,nx) = 250.0_ireals    ! average initial density
                wtot_snow    (i,j,ksn,nx) = w_snow(i,j,nx)/ke_snow
                dzh_snow     (i,j,ksn,nx) = w_snow(i,j,nx)*rho_w/rho_snow_mult(i,j,ksn,nx)/ke_snow
                wliq_snow    (i,j,ksn,nx) = 0.0_ireals
              END DO
            ELSE
              rho_snow(i,j,nx) = 250.0_ireals    ! average initial density
            ENDIF
          ENDIF   ! land-points
        ENDDO
      ENDDO
    ELSE
      IF(lmulti_snow) THEN
!The sum of wtot_snow, i.e. the "old" total water equivalent depth
        zw_snow_old  (:,:) = 0.0_ireals
        zrho_snow_old(:,:) = 0.0_ireals
        DO ksn = 1, ke_snow
          DO   j = jstarts, jends
            DO i = istarts, iends
              zw_snow_old(i,j)   =   zw_snow_old(i,j) +   wtot_snow(i,j,ksn,nx)
              zrho_snow_old(i,j) = zrho_snow_old(i,j) + dzh_snow(i,j,ksn,nx)
            END DO
          END DO
        END DO
!The "old" snow density
        DO   j = jstarts, jends
          DO i = istarts, iends
            zrho_snow_old(i,j) = zw_snow_old(i,j) / MAX(zrho_snow_old(i,j),1.E-09_ireals) *rho_w
          END DO
        END DO
        DO   j = jstarts, jends
          DO i = istarts, iends
            zicount1(i,j) = COUNT(wtot_snow(i,j,:,nx).gt.zepsi)
            zicount2(i,j) = COUNT(dzh_snow(i,j,:,nx).gt.zepsi)
          END DO
        END DO
        DO ksn = 1, ke_snow
          DO   j = jstarts, jends
            DO i = istarts, iends
              IF(zicount1(i,j).EQ.ke_snow .AND. zicount2(i,j).EQ.ke_snow) THEN
                ! all snow layers have non-zero thicknesses
                rho_snow_mult(i,j,ksn,nx) = wtot_snow(i,j,ksn,nx)/dzh_snow(i,j,ksn,nx)*rho_w
                wtot_snow(i,j,ksn,nx) = wtot_snow(i,j,ksn,nx)*w_snow(i,j,nx)/zw_snow_old(i,j)
                dzh_snow(i,j,ksn,nx) = dzh_snow(i,j,ksn,nx)*w_snow(i,j,nx)/zw_snow_old(i,j)
                wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx)*w_snow(i,j,nx)/zw_snow_old(i,j)
              ELSE
                IF(rho_snow(i,j,nx) .EQ. 0._ireals) rho_snow(i,j,nx) = 250._ireals
                rho_snow_mult(i,j,ksn,nx) = rho_snow(i,j,nx)
                  wtot_snow(i,j,ksn,nx) = w_snow(i,j,nx)/ke_snow
                dzh_snow(i,j,ksn,nx) = w_snow(i,j,nx)*rho_w/rho_snow_mult(i,j,ksn,nx)/ke_snow
                wliq_snow    (i,j,ksn,nx) = 0.0_ireals
              END IF
            END DO
          END DO
        END DO
      ELSE   !no multi-layer snow model
        DO   j = jstarts, jends
          DO i = istarts, iends
            IF(rho_snow(i,j,nx) .EQ. 0._ireals) rho_snow(i,j,nx) = 250._ireals
          END DO
        END DO
      END IF
    ENDIF

!   Initialization of the local array of the grid in snow
!   --------------------------------------------
    IF(lmulti_snow) THEN
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF(llandmask(i,j)) THEN   ! for land-points only
            h_snow(i,j,nnew) = 0.0_ireals
            DO ksn = 1,ke_snow
              zdzh_snow(i,j,ksn) = dzh_snow(i,j,ksn,nx)
              h_snow(i,j,nnew) = h_snow(i,j,nnew) + zdzh_snow(i,j,ksn)
            END DO
            h_snow(i,j,nnow) = h_snow(i,j,nnew)

            zhh_snow(i,j,1) = - h_snow(i,j,nnow) + dzh_snow(i,j,1,nx)
            DO ksn = 2,ke_snow
              zhh_snow(i,j,ksn) = zhh_snow(i,j,ksn-1) + dzh_snow(i,j,ksn,nx)
            END DO

            zhm_snow(i,j,1) = (-h_snow(i,j,nnow) + zhh_snow(i,j,1))/2._ireals
            DO ksn = 2,ke_snow
              zhm_snow(i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals
            END DO
          ENDIF    ! llandmask
        END DO
      END DO
    END IF

    ! LINDA, b
    ke_soil_hy = MAX(ibot_w_so,2) ! um_ak_20120329 technically it is possible to
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF(llandmask(i,j)) THEN   ! for land-points only
          IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
            zsatlev(i,j) = ke_soil_hy
            h_wt(i,j)    = 0.0_ireals
            z_wt(i,j)    = zzhls(ke_soil_hy+1)
            do kso=ke_soil_hy,ke_soil+1
               w_so(i,j,kso,:) = zporv(i,j)*zdzhs(kso)
            end do
          ENDIF
        ENDIF
      END DO
    END DO
    IF (lmelt .AND. .NOT. lmelt_var) THEN
      DO kso   = ke_soil_hy,ke_soil+1
        DO   j = jstarts, jends
          DO i = istarts, iends
            IF (llandmask(i,j)) THEN             ! for land-points only
              w_so_ice(i,j,kso,:) = 0.0_ireals
              IF (t_so(i,j,kso,nx) < t0_melt) THEN
                w_so_ice(i,j,kso,:) = w_so(i,j,kso,nx)
              END IF ! t_so(kso) < t0_melt
            END IF   ! land-points
          END DO
        END DO
      END DO
    END IF           ! lmelt .AND. .NOT. lmelt_var
    IF(lmelt .AND. lmelt_var) THEN
      DO kso   = ke_soil_hy,ke_soil+1
        DO   j = jstarts, jends
          DO i = istarts, iends
            IF (llandmask(i,j)) THEN             ! for land-points only
              IF (t_so(i,j,kso,nx) < (t0_melt-zepsi)) THEN
                zaa    = g*zpsis(i,j)/lh_f
                zw_m(i,j)     = zporv(i,j)*zdzhs(kso)
                zw_m(i,j)   = zw_m(i,j)*                                          &
                  ((t_so(i,j,kso,nx) - t0_melt)/(t_so(i,j,kso,nx)*zaa))**(-zedb(i,j))
                w_so_ice(i,j,kso,:) = MAX (0.0_ireals,w_so(i,j,kso,nx) - zw_m(i,j))
              END IF
            END IF   ! land-points
          END DO
        END DO
      END DO
    END IF           ! lmelt .AND. lmelt_var

  ENDIF             ! ntstep = 0

  IF (ntstep == nstart) THEN

    ! compute orography parameter

    ! enforce periodic boundary conditions for hsurf if required:
    !kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    !CALL exchg_boundaries                                                  &
    !     ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
    !     kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,      &
    !     lperi_x, lperi_y, l2dim, &
    !     10000, .FALSE., ncomm_type, ierror, yerror,                  &
    !     hsurf )

! s_oro is read from input file
!    dy = r_earth * dlat * degrad
!    ooleny = 1./dy
!    DO   j = 2, je-1
!       dx     = dlon * r_earth * degrad * crlat(j,1)
!       len    = sqrt(dx**2+dy**2)
!       oolen  = 1./len
!       oolenx = 1./dx
!       DO i = 2, ie-1
!          IF(llandmask(i,j)) THEN   ! for land-points only
!             grad(1)      = oolenx * (hsurf(i,j)-hsurf(i-1,j  ))
!             grad(2)      = oolenx * (hsurf(i,j)-hsurf(i+1,j  ))
!             grad(3)      = ooleny * (hsurf(i,j)-hsurf(i  ,j-1))
!             grad(4)      = ooleny * (hsurf(i,j)-hsurf(i  ,j+1))
!             grad(5)      = oolen  * (hsurf(i,j)-hsurf(i-1,j-1))
!             grad(6)      = oolen  * (hsurf(i,j)-hsurf(i-1,j+1))
!             grad(7)      = oolen  * (hsurf(i,j)-hsurf(i+1,j-1))
!             grad(8)      = oolen  * (hsurf(i,j)-hsurf(i+1,j+1))
!             grad(9)      = 0.00_ireals
!             s_oro(i,j)   = maxval(grad)
!!             s_oro(i,j)   = zsoro(i,j)
!             if (my_cart_id==3) then
!!!                !print*,zsoro(istarts,jstarts),zsoro(iends,jends), my_cart_id, grad,hsurf(i,j),hsurf(i  ,j+1)
!!                print*,i,j,hsurf(i,j),hsurf(i  ,j+1),hsurf(i+1  ,j)
!!             end if
!          ENDIF
!       END DO
!    END DO

     ! limit s_oro to minimal value
     DO   j = jstarts, jends
       DO i = istarts, iends
         IF(llandmask(i,j)) THEN   ! for land-points only
             s_oro(i,j)   = max(s_oro(i,j),s_oro_min)
          ENDIF
       END DO
    END DO
!    ! LINDA, e
  ENDIF             ! ntstep = 0

! End of timestep 0 preparations
! ==============================

! To ensure t_s = t_so(1) also at gme2lm-modified lateral boundaries by
! modifying the soil temperature profile:
  DO kso   = 1,ke_soil+1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN             ! for land-points only
        ! Next lines have to be modified/deleted if the soil surface temperature
        ! t_so(i,j,0,nx) predicted by the heat conduction equation is used
          zdelt_bound = t_s(i,j,nx) - t_so(i,j,1,nx)
          t_so(i,j,kso,nx) = t_so(i,j,kso,nx) + zdelt_bound * (1._ireals - &
                           (zmls(kso) - zmls(1))/(zmls(ke_soil+1) - zmls(1)))
          IF(kso.EQ.1) t_so(i,j,0,nx) = t_so(i,j,1,nx)
        ENDIF
      END DO
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
        ztrangs(i,j)      = 0.0_ireals
        zsnull(i,j)       = 0.0_ireals
        zs1(i,j)          = 0.0_ireals
        zwroot(i,j)       = 0.0_ireals
        zdwidt (i,j)      = 0.0_ireals         ! Initialisation of all
        zdwsndt(i,j)      = 0.0_ireals         ! evaporation quantities
        zesoil (i,j)      = 0.0_ireals         ! to be equal to zero
        zrr    (i,j)      = 0.0_ireals         ! in first part formation of dew
        zrs    (i,j)      = 0.0_ireals         ! in first part formation of rime
        zw_fr(i,j,ke_soil+1)  = w_so(i,j,ke_soil+1,nx)/zdzhs(ke_soil+1)
    END DO
  END DO


  DO kso   = 1, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
        zw_fr(i,j,kso)    = w_so(i,j,kso,nx)/zdzhs(kso)
        ztrang (i,j,kso)  = 0._ireals
        zropartw(i,j,kso) = 0._ireals
      END DO
    END DO
  END DO

! Determine the layer indices and total depth for water content averaging
! (in soil evaporation after Dickinson)
  k10cm  = 1
  k100cm = 1
  z1     = zzhls(1)
  znull  = zzhls(1)
  DO kso = 1, ke_soil
    IF (zmls(kso).le.0.1) THEN
      z1      = zzhls(kso)
      k10cm   = kso
    END IF
    IF (zmls(kso).le.1.0_ireals) THEN
      znull   = zzhls(kso)
      k100cm  = kso
    END IF
  END DO

! Determine average soil water content over 10 and 100 cm, respectively.
  DO kso   = 1, k100cm
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN   ! for land-points only
          IF (kso.le.k10cm) zs1(i,j) = zs1(i,j) + w_so(i,j,kso,nx)
          zsnull(i,j)   = zsnull(i,j) + w_so(i,j,kso,nx)
        END IF
      END DO
    END DO
  END DO

IF (itype_heatcond == 2) THEN

! heat conductivity dependent on actual soil water content
! based on Peters-Lidard et al. (1998) and Johansen (1975)
! see also Block, Alexander (2007), Dissertation BTU Cottbus
    
  DO kso = 1, ke_soil+1
   DO j = jstarts, jends
    DO i = istarts, iends
     IF (llandmask(i,j)) THEN          ! land-points only
      zthetas = zporv(i,j)

      zlamli  = 0.57_ireals
      zlamic  = 2.2_ireals
      zthliq  = zthetas - (w_so_ice(i,j,kso,nx)/zdzhs(kso))
      zlamq   = 7.7_ireals

      rsandf=zsandf(i,j)/100._ireals

      if (rsandf >= 0.2_ireals)  zlam0    =  2.0_ireals
      if (rsandf < 0.2_ireals)   zlam0    =  3.0_ireals


      zlams   = zlamq**rsandf* zlam0**(1._ireals-rsandf)

!      zlamsat = (clamso(m_styp(i,j))**(1.0_ireals-zthetas)) * (zlamli**zthliq) * (zlamic**(zthetas-zthliq))
      zlamsat = (zlams**(1.0_ireals-zthetas)) * (zlamli**zthliq) * (zlamic**(zthetas-zthliq))

      zrhod   = 2700.0_ireals*(1.0_ireals-zthetas)
      zlamdry = ( 0.135_ireals*zrhod + 64.7_ireals )                     &
              / ( 2700.0_ireals - 0.947*zrhod )

      zsri    = MIN(1.0_ireals, (w_so(i,j,kso,nx)/zdzhs(kso)) / zthetas)

      IF ( zsri >= 0.1_ireals ) THEN
       IF ( t_so(i,j,kso,nx) < t0_melt) THEN
        zKe = zsri
       ELSE
        zKe = LOG10(zsri) + 1.0_ireals
       ENDIF
       zKe = MAX(0.0_ireals, (MIN(1.0_ireals, zKe)) )
       hzalam(i,j,kso) = zKe*(zlamsat - zlamdry) + zlamdry
      ELSE
       hzalam(i,j,kso) = zlamdry
      ENDIF
     END IF  ! land-points
    ENDDO
   ENDDO
  ENDDO

  DO kso = 1, ke_soil
   DO j = jstarts, jends
    DO i = istarts, iends
     IF (llandmask(i,j)) THEN          ! land-points only
      ! mean heat conductivity
      zalam(i,j,kso) = hzalam(i,j,kso)*hzalam(i,j,kso+1)*(zmls(kso+1)-zmls(kso))    &
                     / ( hzalam(i,j,kso)*(zmls(kso+1)-zzhls(kso))                   &
                     +   hzalam(i,j,kso+1)*(zzhls(kso)-zmls(kso)) )
     END IF  ! land-points
    ENDDO
   ENDDO
  ENDDO

ELSE

! heat conductivity based on assumption of a soil water content which is equal to the
! average between wilting point and field capacity
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        zwqg         = 0.5_ireals*(zfcap(i,j) + zpwp(i,j))
        z4wdpv       = 4._ireals*zwqg/zporv(i,j)
        ! heat conductivity
        zalamtmp(i,j) =              zdlam(i,j)                         &
                      * (0.25_ireals + 0.30_ireals*zdlam(i,j)           &
                      / (1._ireals+0.75_ireals*zdlam(i,j)))             &
                      * MIN (z4wdpv, 1.0_ireals + (z4wdpv-1.0_ireals)   &
                      *(1.0_ireals+0.35_ireals*zdlam(i,j))              &
                      /(1.0_ireals+1.95_ireals*zdlam(i,j)))
      END IF  ! land-points
    ENDDO
  ENDDO
  DO kso = 1, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          zalam(i,j,kso) = zalam(i,j,kso) + zalamtmp(i,j)
        END IF  ! land-points
      ENDDO
    ENDDO
  ENDDO

ENDIF

#ifdef NECSX
  CALL collapse(.FALSE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif
! Initialisations and conversion of tch to tmch

  limit_tch(:,:) = .false.  !  preset problem indicator
  ldebug         = .false.
  ! LINDA, b
  limit          = .false.
  ! LINDA, e

  DO   j = jstarts, jends
    jm1  = MAX( 1, j-1 )
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN     ! for land-points only
        im1        = MAX( 1, i-1)
        zuv        = 0.5*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                               +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )
        ztvs       = t_g (i,j,nx)*(1.0_ireals + rvd_m_o*qv_s(i,j,nx))
        !  'potential' temperature of lowest atmospheric layer
        zplow          = p0(i,j,ke) + pp(i,j,ke,nx)
        zth_low (i,j)  =  t(i,j,ke,nx) *( (ps(i,j,nx)/zplow)**rdocp )
        zdt_atm (i,j)  =  zth_low(i,j)-t_g(i,j,nx)
        zdq_atm (i,j)  =  qv(i,j,ke)-qv_s(i,j,nx)

!   introduce an artificical upper boundary on transfer coefficients in cases
!   where an extreme cooling/heating of topmost soil layer may be caused due
!   to excessive sensible/latent heat fluxes (e.g. after creation of unbalanced
!   structures in the data assimilation or following massive changes in radiative
!   forcing due to infrequent radiation computations)

!   estimate current energy budget of topmost soil layer

!   heat flux between layers 1&2 based on current temperature profile
    zg1(i,j)= zalam(i,j,1)*(t_so(i,j,1,nx)-t_so(i,j,2,nx))/zdzms(2)

!   estimates of sensible and latent heat flux
    zrho_atm(i,j)=ps(i,j,nx)/(r_d*ztvs)
    zshfl(i,j) = tch(i,j)*zuv*zrho_atm(i,j)*cp_d*zdt_atm(i,j)
    zlhfl(i,j) = tch(i,j)*zuv*zrho_atm(i,j)*lh_v*zdq_atm(i,j)

    IF (ABS(zshfl(i,j)) <= zepsi) zshfl(i,j)=SIGN(zepsi,zshfl(i,j))

    zthfl(i,j) = zshfl(i,j)+zlhfl(i,j)

!   net radiative fluxes at surface
    zradfl(i,j) = sobs(i,j)+thbs(i,j)

!   unconstrained estimated energy budget of topmost soil layer
    zeb1(i,j) = zshfl(i,j)+zlhfl(i,j)+zradfl(i,j)-zg1(i,j)

!   energy required to melt existing snow
    ze_melt(i,j)=w_snow(i,j,nx)*rho_w*lh_f     ! (J/m**2)
!   heat capacity of snow layer
    zch_snow(i,j)=w_snow(i,j,nx)*rho_w*chc_i   ! (J/(m**2 K))

!   constrain transfer coefficient, if energy budget  of topmost soil layer is:
!   a) negative & surface layer is unstable (i.e   upward directed turbulent heat flux)
!   b) positive & surface layer is stable   (i.e downward directed turbulent heat flux)

  IF (zeb1(i,j)<0.0_ireals .AND. zthfl(i,j)<0.0_ireals) THEN 
    ! cooling of 1st soil layer&upward SHF+LHF

    ztchv_max(i,j) =  &
               (-zlim_dtdt*(zch_soil*zdzhs(1)+zch_snow(i,j))/zdt            &
                +zg1(i,j)-zradfl(i,j) )                                     &
               / zthfl(i,j) * tch(i,j) * zuv
  ELSEIF (zeb1(i,j)>0.0_ireals .AND. zthfl(i,j)>0.0_ireals) THEN 
    ! heating of 1st soil layer & downward SHF+LHF
    !   Note: The heat capacity of snow is only relevant for the difference 
    !         between the actual temperature and the melting point. The mean 
    !         snow temperature is set to the average of t_snow & t_so(1)
    IF(lmulti_snow) THEN
      zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow_mult(i,j,1,nx)+t_so(i,j,1,nx)))
    ELSE
      zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow(i,j,nx)+t_so(i,j,1,nx)))
    ENDIF
    ztchv_max(i,j) =  &
    ( (zlim_dtdt*zch_soil*zdzhs(1)+zdT_snow*zch_snow(i,j)+ze_melt(i,j))/zdt  &
      + zg1(i,j)-zradfl(i,j) )                                               &
                / zthfl(i,j) * tch(i,j) * zuv
  ELSE                                                    
    ! unlimited transfer coefficient
    ztchv_max(i,j) = HUGE(1._ireals)
  ENDIF

                                                  ! required constraint as non-turbulent
                                                  ! energy budget components alone may
  ztchv_max(i,j) =MAX( ztchv_max(i,j) ,zepsi)     ! exceed the upper limit in the energy
                                                  ! budget of the 1st soil layer

  ztchv(i,j)    = tch(i,j)*zuv  ! transfer coefficient * velocity

        LIM: IF (ztchv(i,j) > ztchv_max(i,j)) THEN
          tch(i,j)=ztchv_max(i,j)/MAX(zuv,1.E-06_ireals)
!         IF (ntstep > 10) THEN          ! control print only after initial adaptation
!                                        ! to analysis state
            limit_tch(i,j) = .true.      ! set switch for later use
        END IF LIM

        ztmch(i,j) = tch(i,j)*zuv*g*zrho_atm(i,j)
      ENDIF
    ENDDO
  ENDDO

! counter for limitation of transfer coefficients
  m_limit = COUNT( limit_tch(:,:) ) 
!*DL: Reduce dayfile printout
! IF (m_limit > 0) THEN
!   WRITE(*,'(1X,A,I8,A,I12)')                        &
!   ' time step:',ntstep,                             &
!   ' no. of transfer coefficient limitations:',m_limit
! ENDIF


! In debugging mode and if transfer coefficient occured for at least one grid point
  IF (m_limit > 0 .AND. ldebug) THEN
    WRITE(*,'(1X,A,/,2(1X,A,F10.2,A,/),1X,A,F10.2,/,1X,A,F10.3,/)')                  &
           'terra1: transfer coefficient had to be constrained',                  &
           'model time step                                 :', zdt     ,' seconds', &
           'max. temperature increment allowed per time step:',zlim_dtdt,' K',       &
           'upper soil model layer thickness                :', zdzhs(1)
    DO j = jstarts, jends
    jm1  = MAX( 1, j-1 )
      DO i = istarts, iends
      im1 = MAX(1,i-1)
      IF (limit_tch(i,j)) THEN
        zuv        = 0.5*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                               +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )
        yhc       ='COOLING'
        IF (zeb1(i,j) > 0._ireals) Yhc='HEATING'
        PRINT 7001,ntstep,i,j,my_cart_id,Yhc,ps(i,j,nx),  &
              t(i,j,ke,nx), zrho_atm(i,j),     &
              w_snow(i,j,nx),zch_snow(i,j),ze_melt(i,j), &
              zth_low(i,j),zdt_atm(i,j),       &
              t_g   (i,j,nx),t_so(i,j,2,nx),t_so(i,j,3,nx), &
              qv(i,j,ke)*1000._ireals,qv_s(i,j,nx)*1000._ireals, &
              w_so(i,j,1,nx)/zdzhs(1),w_so(i,j,2,nx)/zdzhs(2), &
              w_so(i,j,3,nx)/zdzhs(3),                               &
              u(i,j,ke,nx),v(i,j,ke,nx),zuv, &
              sobs(i,j),thbs(i,j),   &
              zshfl(i,j),zlhfl(i,j),           &
              zthfl(i,j),zg1(i,j),zeb1(i,j), &
              ztchv(i,j),ztchv_max(i,j)
     7001 FORMAT(1x,'time step:',I5,'   grid point:',3I4,2X,A8,/,&
                 1X,'surface pressure   :',F12.2,' Pa'     ,/, &
                 1X,'t(ke)              :',F12.2,' K '     ,/, &
                 1X,'air density        :',F12.2,' kg/m**3',/, &
                 1X,'snow water content :',F12.4,' mH2O   ',/, &
                 1X,'snow heat capacity :',F12.2,' J/(m**2 K)',/, &
                 1X,'snow melting energy:',F12.2,' J/m**2   ',/, &
                 1X,'th_low             :',F12.2,' K '     ,/, &
                 1X,'th_low-t_g         :',F12.2,' K '     ,/, &
                 1X,'t_g                :',F12.2,' K '     ,/, &
                 1X,'t_so(2)            :',F12.2,' K '     ,/, &
                 1X,'t_so(3)            :',F12.2,' K '     ,/, &
                 1X,'qv(ke)             :',F12.4,' g/kg'   ,/, &
                 1X,'qv_s               :',F12.4,' g/kg'   ,/, &
                 1X,'zwg_fr(1)          :',F12.2,' - '     ,/, &
                 1X,'zwg_fr(2)          :',F12.2,' - '     ,/, &
                 1X,'zwg_fr(3)          :',F12.2,' - '     ,/, &
                 1X,'u(ke)              :',F12.2,' m/s'    ,/, &
                 1X,'v(ke)              :',F12.2,' m/s'    ,/, &
                 1X,'|u|                :',F12.2,' m/s'    ,/, &
                 1X,'sobs               :',F12.2,' W/m**2' ,/, &
                 1X,'thbs               :',F12.2,' W/m**2' ,/, &
                 1X,'zshfl              :',F12.2,' W/m**2' ,/, &
                 1X,'zlhfl              :',F12.2,' W/m**2' ,/, &
                 1X,'zthfl              :',F12.2,' W/m**2' ,/, &
                 1X,'zg1                :',F12.2,' W/m**2' ,/, &
                 1X,'zeb1               :',F12.2,' W/m**2' ,/, &
                 1X,'tch*|u| original   :',F12.6,' m/s'    ,/, &
                 1X,'tch*|u| after limit:',F12.6,' m/s'    )
        END IF
      END DO
    END DO
  ENDIF

#ifdef NECSX
  CALL collapse(.TRUE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif
! Update indicator for age of snow in top of snow layer

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF(llandmask(i,j)) THEN        ! for land-points only
        IF (w_snow(i,j,nx) <=0.0_ireals) THEN
          ! if no snow exists, reinitialize age indicator
          freshsnow(i,j) = 1.0_ireals
        ELSE
          IF ( itype_gscp >= 2000 ) THEN
            ! only possible when running 2-moment microphysics
            zsnow_rate = prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j)+prh_gsp(i,j) ! [kg/m**2 s]
          ELSEIF ( itype_gscp >= 4 ) THEN
            zsnow_rate = prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j) ! [kg/m**2 s]
          ELSE
            zsnow_rate = prs_gsp(i,j)+prs_con(i,j)              ! [kg/m**2 s]
          ENDIF

          ! linear decay equals 1.0 in 28 days
          zdsn_old   = zdt/86400._ireals/28._ireals

          ! linear growth rate equals 1.0 in 1 day for a constant
          ! snow rate of 5 mmH2O (kg/m**2) per day
          zdsn_new   = zdt*zsnow_rate*0.2_ireals   !

          ! reduce decay rate, if new snow is falling and as function of snow
          ! age itself
          zdsn_old   = (zdsn_old - zdsn_new)*freshsnow(i,j)
          zdsn_old   = MAX(zdsn_old,0._ireals)

          freshsnow(i,j) = freshsnow(i,j) + zdsn_new-zdsn_old

          freshsnow(i,j) = MIN(1._ireals,MAX(0._ireals,freshsnow(i,j)))
        END IF
      END IF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section I.2: temperatures, water contents (in mH2O), surface pressure, 
!------------------------------------------------------------------------------

  IF(lmulti_snow) THEN
    ! US: this loop does not vectorize until the ksn-loops are put outside
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN     ! for land-points only

          IF (w_snow(i,j,nx) > 0.0_ireals) THEN
            ! existence of snow
            ! --> no water in interception store and t_snow < t0_melt
            w_snow(i,j,nx) = w_snow(i,j,nx) + w_i(i,j,nx)
            wtot_snow(i,j,1,nx) = wtot_snow(i,j,1,nx) + w_i(i,j,nx)
            dzh_snow(i,j,1,nx) = dzh_snow(i,j,1,nx) + w_i(i,j,nx)/rho_snow_mult(i,j,1,nx)*rho_w
            w_i   (i,j,nx) = 0.0_ireals
            DO ksn = 0,ke_snow
              t_snow_mult(i,j,ksn,nx) = MIN (t0_melt - zepsi, t_snow_mult(i,j,ksn,nx) )
            END DO
          ELSE IF (t_snow_mult(i,j,ke_snow,nx) >= t0_melt) THEN
            ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
            t_s   (i,j,nx) = MAX (t0_melt + zepsi, t_s(i,j,nx) )
            DO ksn = 0,ke_snow
              t_snow_mult(i,j,ksn,nx) = t_s(i,j,nx)
            END DO
          ELSE
            ! no snow and  t_snow < t0_melt
            ! --> t_snow = t_s and no water w_i in interception store
            t_s   (i,j,nx) = MIN (t0_melt - zepsi, t_s(i,j,nx) )
            DO ksn = 0,ke_snow
              t_snow_mult(i,j,ksn,nx) = t_s(i,j,nx)
            END DO
            w_i   (i,j,nx) = 0.0_ireals
          END IF
        END IF
      ENDDO
    ENDDO

  ELSE ! no multi-layer snow

    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN     ! for land-points only
          IF (w_snow(i,j,nx) > 0.0_ireals) THEN
            ! existence of snow
            ! --> no water in interception store and t_snow < t0_melt
            w_snow(i,j,nx) = w_snow(i,j,nx) + w_i(i,j,nx)
            w_i   (i,j,nx) = 0.0_ireals
            t_snow(i,j,nx) = MIN (t0_melt - zepsi, t_snow(i,j,nx) )
          ELSE IF (t_snow(i,j,nx) >= t0_melt) THEN
            ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
            t_s   (i,j,nx) = MAX (t0_melt + zepsi, t_s(i,j,nx) )
            t_snow(i,j,nx) = t_s(i,j,nx)
          ELSE
            ! no snow and  t_snow < t0_melt
            ! --> t_snow = t_s and no water w_i in interception store
            t_s   (i,j,nx) = MIN (t0_melt - zepsi, t_s(i,j,nx) )
            t_snow(i,j,nx) = t_s(i,j,nx)
            w_i   (i,j,nx) = 0.0_ireals
          END IF
        END IF
      ENDDO
    ENDDO
  ENDIF

  !Inizializations for the next sections
  IF (lmulti_snow) THEN
    ! some preparations for ksn==0 and ksn==1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN   ! for land-points only
          ztsnow_mult   (i,j,0) = t_snow_mult(i,j,0,nx)
          ztsnow_mult   (i,j,1) = t_snow_mult(i,j,1,nx)

          zhh_snow(i,j,1) =  -h_snow(i,j,nnow) + dzh_snow(i,j,1,nx)
          zhm_snow(i,j,1) = (-h_snow(i,j,nnow) + zhh_snow(i,j,1))/2._ireals

          zdzh_snow(i,j,1) = dzh_snow(i,j,1,nx)
          zextinct (i,j,1) = 0.13_ireals*rho_snow_mult(i,j,1,nx)+3.4_ireals

          ! set ztsnow to ztsnow_mult(ksn=1)
          ztsnow   (i,j) = t_snow_mult(i,j,1,nx)

          ! initialize zdz_snow (from Section I.3) to zdzh_snow(i,j,1)
          zdz_snow (i,j) = dzh_snow(i,j,1,nx) ! zdzh_snow
          zwsnow   (i,j) = w_snow(i,j,nx)
        END IF
      ENDDO
    ENDDO

    DO  ksn = 2,ke_snow
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN   ! for land-points only
            ztsnow_mult(i,j,ksn) = t_snow_mult(i,j,ksn,nx)

            zhh_snow   (i,j,ksn) = zhh_snow(i,j,ksn-1) + dzh_snow(i,j,ksn,nx)
            zhm_snow   (i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals

            zdzh_snow  (i,j,ksn) = dzh_snow(i,j,ksn,nx)
            zextinct   (i,j,ksn) = 0.13_ireals*rho_snow_mult(i,j,ksn,nx)+3.4_ireals

            ! build sum over all layers for zdzh_snow in zdz_snow
            zdz_snow   (i,j)     = zdz_snow(i,j) + zdzh_snow(i,j,ksn)
          END IF
        ENDDO
      ENDDO
    ENDDO
  ELSE
    ! set ztsnow to t_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN   ! for land-points only
          ztsnow   (i,j) = t_snow(i,j,nx)
          zwsnow   (i,j) = w_snow(i,j,nx)
          zdz_snow (i,j) = zwsnow(i,j)*rho_w/rho_snow(i,j,nx)
        END IF
      ENDDO
    ENDDO
  ENDIF

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN   ! for land-points only
        ! ztsnow is now set according to lmulti_snow (see above);
        ! therefore we need not take care of it in this loop
        ! ztsnow   (i,j) = t_snow(i,j,nx)
        zts      (i,j) = t_s   (i,j,nx)
        zts_pm   (i,j) = zsf_heav(zts   (i,j) - t0_melt)
        ztsnow_pm(i,j) = zsf_heav(ztsnow(i,j) - t0_melt)
        zwin     (i,j) = w_i   (i,j,nx)

        ! moisture and potential temperature of lowest atmospheric layer
        zplow          = p0(i,j,ke) + pp(i,j,ke,nx)
        zqvlow         = qv(i,j,ke)
        zth_low (i,j)  =  t(i,j,ke,nx) *( (ps(i,j,nx)/zplow)**rdocp )

        ! density*transfer coefficient*wind velocity
        zrhoch(i,j)    = ztmch(i,j)*(1._ireals/g) + zepsi

        ! saturation specific humidity for t_s and t_snow and first derivative
        z2iw        = zts_pm(i,j)*b2w + (1._ireals - zts_pm(i,j))*b2i
        z4iw        = zts_pm(i,j)*b4w + (1._ireals - zts_pm(i,j))*b4i
        z234iw      = z2iw*(b3 - z4iw)
        zqs         = zsf_qsat( zsf_psat_iw(zts(i,j), z2iw,z4iw), ps(i,j,nx) )
        zdqs        = zqvlow - zqs
        IF (ABS(zdqs).LT.zepsi) zdqs = 0._ireals
        z2iw        = ztsnow_pm(i,j)*b2w + (1._ireals - ztsnow_pm(i,j))*b2i
        z4iw        = ztsnow_pm(i,j)*b4w + (1._ireals - ztsnow_pm(i,j))*b4i
        z234iw      = z2iw*(b3 - z4iw)
        zqsnow      = zsf_qsat(zsf_psat_iw(ztsnow(i,j),z2iw,z4iw), ps(i,j,nx))
        zdqvtsnow(i,j)= zsf_dqvdt_iw(ztsnow(i,j), zqsnow, z4iw,z234iw)
        zdqsnow     = zqvlow - zqsnow
        IF (ABS(zdqsnow).LT.zepsi) zdqsnow = 0._ireals

        ! potential evaporation at T_snow and Ts
        zep_snow(i,j) = (1._ireals-ztsnow_pm(i,j))*                       &
                                          tfv(i,j)*zrhoch(i,j)*zdqsnow
        zep_s   (i,j) =                   tfv(i,j)*zrhoch(i,j)*zdqs
      END IF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section I.3: heat conductivity, frozen fraction, snow and
!            water covered fraction, snow density and  height,
!            volumetric heat content
!------------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        ! snow and water covered fraction
        zrss = MAX( 0.01_ireals, MIN(1.0_ireals,zwsnow(i,j)/cf_snow) ) 
        zrww = MAX( 0.01_ireals, 1.0_ireals -                           &
                                 EXP(MAX( -5.0_ireals, - zwin(i,j)/cf_w) ) )
        zf_snow(i,j) = zrss*zsf_heav(zwsnow(i,j) - zepsi)
        zf_wi  (i,j) = zrww*zsf_heav(zwin  (i,j) - zepsi)

! BR 7/2005 prognostic  snow density
! US DMironov: for the FLake Model, h_snow has to be a prognostic variable but
!              here it is used only diagnostically. So the values are put
!              to every time level
        ! zdz_snow has been computed in the initializations above depending
        ! on lmulti_snow
        ! zdz_snow(i,j)=zwsnow(i,j)*rho_w/rho_snow(i,j,nx)
        h_snow(i,j,nnew) = zdz_snow(i,j)
        h_snow(i,j,nnow) = zdz_snow(i,j)
        IF (.NOT. l2tls) h_snow(i,j,nold) = zdz_snow(i,j)

!       constrain snow depth and consider this constraint for the computation
!       of average snow density of snow layer
        zdz_snow (i,j) =  zdz_snow (i,j)/MAX(0.01_ireals,zf_snow(i,j))
        zdz_snow (i,j) =  MAX(cdsmin,zdz_snow(i,j))

!       limitation of snow depth to 1.5m for snow cover heat transfer
        zdz_snow_fl(i,j) = MIN(1.5_iREALs, zdz_snow(i,j))
        IF(.NOT. lmulti_snow) &
          zrocs(i,j) = chc_i*zdz_snow_fl(i,j)*rho_snow(i,j,nx)
!BR 7/2005 End
      END IF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section I.4: Hydrology, 1.Section
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section I.4.1: Evaporation from interception store and from snow cover,
  
  !----------------------------------------------------------------------------
  ! Evaporation and transpiration are negative, dew and rime
  ! positive quantities, since positive sign indicates a flux
  ! directed towards the earth's surface!

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN             ! land points only
        ! Evaporation from interception store if it contains water (wi>0) and
        ! if zep_s<0 indicates potential evaporation for temperature Ts
        ! amount of water evaporated is limited to total content of store
        zdwidt(i,j) = zsf_heav(-zep_s(i,j))      &
                      * MAX(-zrhwddt*zwin(i,j), zf_wi(i,j)*zep_s(i,j))
        ! Evaporation of snow, if snow exists (wsnow>0) and if zep_snow<0
        ! indicates potential evaporation for temperature t_snow
        zdwsndt(i,j) = zsf_heav(-zep_snow(i,j))  &
                       * MAX(-zrhwddt*zwsnow(i,j), zf_snow(i,j)*zep_snow(i,j))
        ! Formation of dew or rime, if zep_s > 0 . distinction between
        ! dew or rime is only controlled by sign of surface temperature
        ! and not effected by presence of snow !
        zrr(i,j)=zsf_heav(zep_s   (i,j))*zep_s   (i,j)*            zts_pm(i,j)
        zrs(i,j)=zsf_heav(zep_snow(i,j))*zep_snow(i,j)*(1.0_ireals-zts_pm(i,j))
      END IF
    ENDDO
  ENDDO
  

  !----------------------------------------------------------------------------
  ! Section I.4.2a: Bare soil evaporation, bucket version
  !----------------------------------------------------------------------------
  
  IF (itype_evsl.EQ.1) THEN   ! Bucket version
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN             ! land-points only
          IF (zep_s(i,j) < 0.0_ireals) THEN  ! upwards directed potential
                                             ! evaporation
            ! reduction factor for evaporation based on water content of
            ! first soil layer, air dryness point and field capacity
            zbeta = MAX( 0.0_ireals, MIN( 1.0_ireals,                        &
                    (zw_fr(i,j,1)-zadp(i,j)) / (zfcap(i,j)-zadp(i,j))))
            zbeta = zbeta**2

            ! if first soil layer is frozen, allow evaporation at potential
            ! rate; if soil type is rock, evaporation is not allowed at all
            zice        = zsf_heav(1.5_ireals - m_styp(i,j)) ! 1 only for ice
            zevap       = zrock(i,j) + zice          ! 1 for all, but rock
            zbeta  = zbeta + (1._ireals - zbeta)*zice
            zesoil(i,j) = zevap*zbeta*zep_s(i,j)      & ! evaporation
                          *(1._ireals - zf_wi  (i,j)) & ! not water covered
                          *(1._ireals - zf_snow(i,j)) & ! not snow covered
                          * eai(i,j)/sai(i,j)        ! relative source surface
                                                     !  of the bare soil
            lhfl_bs(i,j) = lh_v * zesoil(i,j)
          END IF ! upwards directed potential evaporation
        END IF   ! land points
      END DO
    END DO
  END IF         ! Bucket version
  
  !----------------------------------------------------------------------------
  ! Section I.4.2b: Bare soil evaporation, BATS version
  !----------------------------------------------------------------------------
  IF (itype_evsl.EQ.2) THEN
    ! Calculation of bare soil evaporation after Dickinson (1984)
    ! Determination of mean water content relative to volume of voids
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN       ! land points only
          IF (zep_s(i,j) < 0.0_ireals) THEN   ! upwards directed potential
                                              ! evaporation
            zsnull(i,j) = zsnull(i,j)/(znull*zporv(i,j))
            ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
            zice   = zsf_heav(1.5_ireals - m_styp(i,j)) ! 1 only for ice
            zevap  = zrock(i,j) + zice                  ! 1 for all soil types
                                                        ! but rock and ice (=0)
            zbeta  = 0.0_ireals
            IF (m_styp(i,j).ge.3) THEN ! Computations not for ice and rocks
              ! auxiliary quantities
              zbf1   = 5.5_ireals - 0.8_ireals* zbedi(i,j)*                &
                      (1.0_ireals + 0.1_ireals*(zbedi(i,j) - 4.0_ireals)*  &
                       clgk0(m_styp(i,j)) )
              zbf2   = (zbedi(i,j) - 3.7_ireals + 5.0_ireals/zbedi(i,j))/  &
                      (5.0_ireals + zbedi(i,j))
              zdmax  = zbedi(i,j)*cfinull*zk0di(i,j)/crhowm 
              zs1(i,j)  = zs1(i,j)/(z1*zporv(i,j))
              zd     = 1.02_ireals*zdmax*zs1(i,j)**(zbedi(i,j) + 2) *      &
                                                (zsnull(i,j)/zs1(i,j))**zbf1
              zck    = (1.0_ireals + 1550.0_ireals*cdmin/zdmax)*zbf2
              ! maximum sustainable moisture flux in the uppermost surface
              ! layer in kg/(s*m**2)
              zfqmax = - rho_w*zck*zd*zsnull(i,j)/SQRT(znull*z1)
              zevapor= MAX(zep_s(i,j),zfqmax)
              IF(zw_fr(i,j,1)+zevapor*(1.0_ireals - zf_wi(i,j))   &
                                     *(1.0_ireals - zf_snow(i,j)) &
                                     *eai(i,j)/sai(i,j)* zdtdrhw/zdzhs(1) &
                                     .LE.zadp(i,j)) zevapor = 0._ireals
              zbeta  = zevapor/MIN(zep_s(i,j),-zepsi)
            END IF ! Computations not for ice and rocks
            zbeta  = zbeta + (1.0_ireals - zbeta)*zice
            ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
            ! soil types
            ! consideration of plant or snow/water cover
            zesoil(i,j) = zevap*zbeta*zep_s(i,j)       & ! evaporation
                          *(1.0_ireals - zf_wi  (i,j)) & ! not water covered
                          *(1.0_ireals - zf_snow(i,j)) & ! not snow covered
                          * eai(i,j)/sai(i,j) ! relative source surface
                                              ! of the bare soil
            lhfl_bs(i,j) = lh_v * zesoil(i,j)
          END IF  ! upwards directed potential evaporation
        END IF    ! land points
      END DO
    END DO
  END IF ! BATS version

  !----------------------------------------------------------------------------
  ! Section I.4.3a: transpiration by plants, bucket version
  !----------------------------------------------------------------------------

  IF (itype_trvg.EQ.1) THEN   ! Bucket version
  ! for lower layers, even for water saturated soil and and maximum
  ! root extension, transpiration is limited to the potential evaporation rate
  ! the consideration of a root depth (zroot) allows the maximum
  ! transpiration, if the active soil layer is saturated and completely
  ! penetrated by roots

    DO kso = 1,ke_soil  ! loop over soil layers
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN     ! land-points only
            IF (zep_s(i,j) < 0.0_ireals  &    ! potential evaporation and
               .AND. m_styp(i,j) >= 3) THEN   ! neither ice nor rock
              ! turgor-loss-point minus plant wilting point
              ztlpmwp(i,j) = (zfcap(i,j) - zpwp(i,j))*(0.81_ireals +       &
                    0.121_ireals*ATAN(-86400._ireals*zep_s(i,j) - 4.75_ireals))
              ! determine whether neither surface nor lower boundary of
              ! second layer are below freezing point
              ztgt0 = zsf_heav(t_so(i,j,kso,nx) - t0_melt)
              ! determine reduction factor beta for this layer
              zbeta   = ztgt0 * MAX(0.0_ireals, MIN(1.0_ireals,            &
                        (zw_fr(i,j,kso) - w_so_ice(i,j,kso,nx)/zdzhs(kso)  &
                        -  zpwp(i,j)) / ztlpmwp(i,j)))
              zbeta   = zbeta**2

              ! consider the effect of root depth
              zroot = MIN ( zdzhs(kso), MAX(0.0_ireals,                    &
                         zbwt(i,j) - (zmls(kso) - 0.5_ireals*zdzhs(kso))))
              zr_root = zroot/MAX(zepsi,zbwt(i,j))
              ztrang(i,j,kso) =  zr_root*zbeta             & ! reduction
                                   * zep_s(i,j)            & ! transpiration
                                   * (1. - zf_wi(i,j))     & ! non-water
                                   * (1. - zf_snow(i,j))   & ! non-snow
                                   * tai(i,j)/sai(i,j)       ! transp. surface
                                     ! relative source surface of the plants
              lhfl_pl(i,j,kso)= lh_v * ztrang(i,j,kso)
              ztrangs(i,j)    = ztrangs(i,j) + ztrang(i,j,kso)
            END IF  ! upwards directed potential evaporation .AND. m_styp > 2
          END IF    ! land-points only
        END DO
      END DO
    END DO          ! loop over soil layers
  END IF            ! Bucket version


  !----------------------------------------------------------------------------
  ! Section I.4.3b: transpiration by plants, BATS version
  !----------------------------------------------------------------------------

  IF (itype_trvg.EQ.2) THEN   ! BATS version
    ! This version is based on Dickinson's (1984) BATS scheme, simplified by
    ! neglecting the water and energy transports between the soil and the plant
    ! canopy. This leads to a Monteith combination formula for the computation
    ! of plant transpiration.

! Root distribution
    IF (itype_root == 2) THEN
      DO   j = jstarts, jends
        DO i = istarts, iends
          zrootdz_int (i,j)= 0.0_ireals   !initialize the root density profile integral
          zwrootdz_int(i,j)= 0.0_ireals   !initialize the root water content integral
        END DO
      END DO
      DO kso   = 1,ke_soil
        DO   j = jstarts, jends
          DO i = istarts, iends
            IF (llandmask(i,j)) THEN ! land points only,
              IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i,j) < 0.0) THEN  ! upwards directed potential evaporation
                  ! consider the effect of root depth & root density
                  zrootfr = EXP (-zroota(i,j)*zmls(kso)) ! root density
                  zrootdz = zrootfr*MIN(zdzhs(kso),MAX(0.0_ireals, zbwt(i,j)-(zmls(kso) -0.5_ireals*zdzhs(kso)) ) )
                  zrootdz_int(i,j)=zrootdz_int(i,j) + zrootdz
                  zwrootdz(i,j,kso)=zrootdz*(zw_fr(i,j,kso)-w_so_ice(i,j,kso,nx)/zdzhs(kso))
                  zwrootdz_int(i,j)=zwrootdz_int(i,j) + zwrootdz(i,j,kso)
                END IF  ! negative potential evaporation only
              END IF  ! neither ice or rocks
            END IF    ! land-points only
          END DO
        END DO
      END DO

       ! Compute root zone integrated average of liquid water content
       DO   j = jstarts, jends
          DO i = istarts, iends
             zwrootdz_int(i,j)=zwrootdz_int(i,j)/MAX(zrootdz_int(i,j),zepsi)
          END DO
       END DO

    ELSE

      DO kso   = 1,ke_soil
        DO   j = jstarts, jends
          DO i = istarts, iends
            IF (llandmask(i,j)) THEN ! land points only,
              IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i,j) < 0.0_ireals) THEN  ! upwards directed potential
                                                   ! evaporation
                  zropart  = MIN ( zdzhs(kso), MAX(0.0_ireals,                     &
                             zbwt(i,j) - (zmls(kso) - 0.5_ireals*zdzhs(kso))))
                  zropartw(i,j,kso) = zropart*(zw_fr(i,j,kso)-w_so_ice(i,j,kso,nx)/zdzhs(kso))
                  zwroot(i,j) = zwroot(i,j) + zropartw(i,j,kso)/MAX(zepsi,zbwt(i,j))
                END IF  ! negative potential evaporation only
              END IF  ! neither ice or rocks
            END IF    ! land-points only
          END DO
        END DO
      END DO

    ENDIF

    ! Determination of the transfer functions CA, CF, and CV

    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN ! land points only,
          IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
            IF (zep_s(i,j) < 0.0) THEN  ! upwards directed potential evaporation
              im1        = MAX( 1, i-1)
              jm1        = MAX( 1, j-1 )
              zuv        = 0.5*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2 &
                               +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )
              !zuv        = SQRT (u_10m(i,j) **2 + v_10m(i,j)**2 )
              zcatm      = tch(i,j)*zuv           ! Function CA

              !WL2011b
              !This dependency on transfer-scheme is absolutely strange????? Since I am only a PhD-student I do not dare removin it ;-)
              !WL2011e
              IF(itype_tran == 1) THEN
                zustar     = zuv*SQRT(tcm(i,j))
                zrla       = 1.0_ireals/MAX(cdash*SQRT(zustar),zepsi)
              ELSE
                zrla       = 0._ireals
              ENDIF

           
              ! to compute CV, first the stomatal resistance has to be determined
              ! this requires the determination of the F-functions:
              ! Radiation function
              zpar       = pabs(i,j)  !  PAR
              zf_rad(i,j)= MAX(0.0_ireals,MIN(1.0_ireals,zpar/cparcrit))
              ztlpmwp(i,j) = (zfcap(i,j) - zpwp(i,j))*(0.81_ireals +       &
                      0.121_ireals*ATAN(-86400._ireals*zep_s(i,j) - 4.75_ireals))

              ! Soil water function
              IF (itype_root == 2) THEN
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwrootdz_int(i,j) -  &
                                                zpwp(i,j))/ztlpmwp(i,j)))
              ELSE
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwroot(i,j) -  &
                                                zpwp(i,j))/ztlpmwp(i,j)))
              ENDIF

              ! Temperature function
!             if(ntstep.eq.0.and.itype_tran.ne.2) then
!               t_2m(i,j)=t(i,j,ke,nx)
!             endif
!             zf_tem     = MAX(0.0_ireals,MIN(1.0_ireals,4.0_ireals*     &
!                          (t_2m(i,j)-t0_melt)*(ctend-t_2m(i,j))/(ctend-t0_melt)**2))
              zf_tem     = MAX(0.0_ireals,MIN(1.0_ireals,4.0_ireals*     &
                           (t(i,j,ke,nx)-t0_melt)*(ctend-t(i,j,ke,nx))/(ctend-t0_melt)**2))


              ! Saturation deficit function (function not used, but computations
              ! necessary for determination of  slope of the saturation curve)
              z2iw       = zts_pm(i,j)*b2w + (1._ireals - zts_pm(i,j))*b2i
              z4iw       = zts_pm(i,j)*b4w + (1._ireals - zts_pm(i,j))*b4i
              zepsat     = zsf_psat_iw(t(i,j,ke,nx),z2iw,z4iw)
              zepke      = qv(i,j,ke)*ps(i,j,nx)/                   &
                                        (rdv + o_m_rdv*qv(i,j,ke))
              zf_sat     = MAX(0.0_ireals,MIN(1.0_ireals,1.0_ireals -  &
                                              (zepsat - zepke)/csatdef))
              ! zf_sat paralysed:
              zf_sat     = 1.0_ireals

              IF (lstomata) THEN
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/rsmin2d(i,j) -     &
                             1.0_ireals/crsmax)*zf_rad(i,j)*zf_wat*zf_tem*zf_sat
              ELSE
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/crsmin -     &
                             1.0_ireals/crsmax)*zf_rad(i,j)*zf_wat*zf_tem*zf_sat
              END IF

              zrstom     = 1.0_ireals/zedrstom              ! stomatal resistance
              rstom(i,j) = zrstom
              zrveg      = zrla + zrstom
              ! Transpiration rate of dry leaves:
              ztraleav(i,j)=zep_s(i,j)*tai(i,j)/(sai(i,j)+zrveg*zcatm)
            END IF  ! upwards directed potential evaporation only
          END IF    ! m_styp > 2
        END IF      ! land points
      END DO
    END DO

    ! Consideration of water and snow coverage, distribution to the different
    ! soil layers
    DO     kso       = 1,ke_soil
      DO   j         = jstarts, jends
        DO i         = istarts, iends
          IF (llandmask(i,j)) THEN ! land points only,
            IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
              IF (zep_s(i,j) < 0.0_ireals) THEN    ! upwards potential evaporation
                ztrabpf  = ztraleav(i,j)*                   & ! plant covered part
                           (1.0_ireals - zf_wi(i,j))*       & ! not water covered
                           (1.0_ireals - zf_snow(i,j))        ! not snow covered

                ! for root distribution
                IF (itype_root == 2) THEN
                  ztrfr    = zwrootdz(i,j,kso)/(zrootdz_int(i,j)*zwrootdz_int(i,j))
                  ztrang(i,j,kso) = ztrabpf*ztrfr
                ELSE
                  zrootfc = zropartw(i,j,kso)/(zwroot(i,j) + zepsi)
                  ztrang(i,j,kso) = ztrabpf*zrootfc/MAX(zepsi,zbwt(i,j))
                  IF(zw_fr(i,j,kso)+ztrang(i,j,kso)*zdtdrhw/zdzhs(kso) &
                                    .LT.zpwp(i,j)) ztrang(i,j,kso) = 0._ireals
                ENDIF
                lhfl_pl(i,j,kso)= lh_v * ztrang(i,j,kso)
                ztrangs(i,j)    = ztrangs(i,j) + ztrang(i,j,kso)
              END IF  ! upwards directed potential evaporation only
            END IF    ! m_styp > 2
          END IF      ! land points
        END DO
      END DO
    END DO          ! loop over soil layers

  END IF ! BATS version

  !----------------------------------------------------------------------------
  ! Section I.4.4: total evapotranspiration and 
  !              associated ficticious soil humidity qv_s
  !----------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN   ! land-points only
        ze_sum = zdwsndt(i,j  )  & ! evaporation of snow
               + zdwidt (i,j  )  & ! evaporation from interception store
               + zesoil (i,j  )  & ! evaporation from bare soil
               + ztrangs(i,j  )  & ! transpiration from all soil layers
               + zrr    (i,j  )  & ! formation of dew
               + zrs    (i,j  )    ! formation of rime
        qv_s(i,j,nx  ) = qv (i,j,ke) - ze_sum /(zrhoch(i,j) + zepsi)
        qv_s(i,j,nnew) = qv_s(i,j,nx)
      END IF     ! land points
    END DO
  END DO

!------------------------------------------------------------------------------
! End of former module procedure terra1_multlay
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
! Former SUBROUTINE terra2_multlay (yerror, ierror)
!------------------------------------------------------------------------------

!   In the prognostic part II the equation of heat conduction and water
!   transport is solved for a multi-layer soil using the same vertical grid
!   Freezing/melting of soil water/ice is accounted for (optionally). A
!   simple one-layer snow model provides the snow surface temperature and
!   the snow water equivalent.

!------------------------------------------------------------------------------
! Section II.1: Initializations
!------------------------------------------------------------------------------

  ! Computation of derived constants
  z1d2dt = 1._ireals/zdt      ! 1./2*timestep

  ! Number of soil layers contributing to surface run-off
  msr_off  = 0
!
!HJP 2010-07-05 Begin
!! Note:
!!
!! - Introduction of i250, i.e. index of last soil layer completely above
!!   2.5m soil depth
!! - Usage of i250 to compute subsoil runoff (runoff_g) as drainage flux
!!   through bottom of this layer (as suggested by the Rhone-Aggregation
!!   Experiment)and to control switching off of soil moisture gradient
!!   related flux below (i.e. only sedimentation flux
!!   allowed between i250 and i250+1)

! Note:
!
! - Introduction of ibot_w_so, i.e. index of last soil layer completely
!   above the defined soil depth in namelist variable czbot_w_so.
!   ibot_w_so has already been calculated in organize_physics.f90
! - Usage of ibot_w_so to compute subsoil runoff (runoff_g) as drainage flux
!   through bottom of this layer (as suggested by the Rhone-Aggregation
!   Experiment)and to control switching off of soil moisture gradient
!   related flux below (i.e. only sedimentation flux
!   allowed between ibot_w_so  and ibot_w_so+1

  ! no. of layers above 2.5 m soil depth
! i250 = 0
! do kso=1,ke_soil+1
! if (zzhls(kso)<=2.5) i250=kso
! ke_soil_hy = i250
! end do
!  ke_soil_hy = ibot_w_so
  ke_soil_hy = MAX(ibot_w_so,2) ! um_ak_20120329 technically it is possible to
  ! run the multilayer scheme with ke_soil = 1 => index ke_soil_hy-1 = 0
  ! => index out of bounds below in ziw_fr
  ! Does this physically make sense ?
!HJP 2010-07-05 End

  zdtdrhw = zdt/rho_w   ! timestep/density of liquid water

  ! time constant for infiltration of water from interception store
  ! must not be less than 2*time step
  ctau_i  = MAX( ctau_i, zdt )

  ! Utility variable to avoid IF-constructs
  DO kso     = 1,ke_soil
    znlgw1f(kso)  = 0.0_ireals
  END DO
  znlgw1f(1)         = 1.0_ireals

  ! LINDA, b
  ! diagnose water-table depth
  DO   kso = 1,ke_soil+1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN     ! land-points only
          ziw_fr(i,j,kso) = w_so_ice(i,j,kso,nx)/zdzhs(kso)   ! ice frac.
          zlw_fr(i,j,kso) = zw_fr(i,j,kso) - ziw_fr(i,j,kso)  ! liquid water frac.
        END IF
      END DO
    END DO
  END DO      !soil layers

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zsatlev(i,j) = ke_soil_hy
          do while (((zw_fr(i,j,zsatlev(i,j))+100.*zepsi).ge.zporv(i,j)).and.(zsatlev(i,j).gt.1))
            zsatlev(i,j) = zsatlev(i,j)-1
          end do
          zsatlev(i,j) = zsatlev(i,j)+1
           satlev(i,j) = real(zsatlev(i,j),ireals)
        END IF
      END IF
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          DO   kso =2,ke_soil+1
            if(kso==zsatlev(i,j)-1) then! last partly saturated level above ground water
              !zlw_fr_ksom05 = 0.5*(zdzhs(kso-1)*zlw_fr_kso(i,j,kso)+ &
              !     zdzhs(kso)*zlw_fr_kso(i,j,kso-1))/zdzms(kso)
              h_wt(i,j) = max(0.0_ireals,zdzhs(kso)*(zw_fr(i,j,kso)-zw_fr(i,j,kso-1))/(zporv(i,j)-zw_fr(i,j,kso-1))) ! dz in layer
              z_wt(i,j) = zzhls(kso) - h_wt(i,j) ! depth of water table
              wt_depth(i,j) = z_wt(i,j)
            end if
          END DO
        END IF
      END IF
    END DO
  END DO
  ! LINDA, e

! Initialisations
  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN                 ! land-points
            zdtsnowdt_mult(i,j,ksn)  = 0.0_ireals
            zdtsdt(i,j)     = 0.0
          END IF
        END DO
      END DO
    END DO
  ELSE
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN                 ! land-points
          zdtsnowdt(i,j)  = 0.0
          zdtsdt(i,j)     = 0.0
        END IF
      END DO
    END DO
  END IF


!------------------------------------------------------------------------------
! Section II.2: Prepare basic surface properties and create some local
!               arrays of surface related quantities (for land-points only)
!               Initialise some fields
!------------------------------------------------------------------------------

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN                 ! land-points
      mstyp        = NINT(soiltyp(i,j))
      m_styp(i,j)  = mstyp
      zsandf(i,j)  = csandf(mstyp)
      zclayf(i,j)  = cclayf(mstyp)
      zgsb   (i,j) = 0.0_ireals
      ztrangs(i,j) = 0.0_ireals
    END IF
  END DO
END DO


DO   kso = 1,ke_soil+1
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN     ! land-points only
        ziw_fr(i,j,kso) = w_so_ice(i,j,kso,nx)/zdzhs(kso)   ! ice frac.
        zlw_fr(i,j,kso) = zw_fr(i,j,kso) - ziw_fr(i,j,kso)  ! liquid water frac.
        zroc(i,j,kso)   = zrocg(i,j) + rho_w*zlw_fr(i,j,kso)*chc_w +          &
                                       rho_w*ziw_fr(i,j,kso)*chc_i
        IF (kso<=ke_soil) THEN
          ztrangs(i,j) = ztrangs(i,j) + ztrang(i,j,kso)
          zdwgdt(i,j,kso) = 0.0_ireals
        END IF
        zflmg(i,j,kso)  = 0.0_ireals
        zrunoff_grav(i,j,kso)  = 0.0
      END IF
    END DO
  END DO
END DO      !soil layers

!------------------------------------------------------------------------------
! Section II.3: Estimate thermal surface fluxes
!------------------------------------------------------------------------------
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF(llandmask(i,j))THEN     ! land-points only

        ! store forcing terms due to evapotranspiration, formation of dew
        ! and rime for later use
        zverbo(i,j) = zdwidt(i,j) + zesoil(i,j) + ztrangs(i,j) +              &
                        (1._ireals-zf_snow(i,j))*(zrr(i,j) + zrs(i,j))

        zversn(i,j) = zdwsndt(i,j) + zrs(i,j)                                 &
                                  + zsf_heav (zwsnow(i,j) - zepsi) * zrr(i,j)

        ! add grid scale and convective precipitation (and graupel, if present)
        ! to dew and rime
        zrr(i,j) = zrr(i,j) + prr_con(i,j) + prr_gsp(i,j)
        IF ( itype_gscp >= 2000 ) THEN
          ! only possible when running the 2-moment microphysics
          zrs(i,j) = zrs(i,j) + prs_con(i,j) + prs_gsp(i,j) + prg_gsp(i,j) + prh_gsp(i,j)
        ELSEIF ( itype_gscp >= 4 ) THEN
          zrs(i,j) = zrs(i,j) + prs_con(i,j) + prs_gsp(i,j) + prg_gsp(i,j)
        ELSE
          zrs(i,j) = zrs(i,j) + prs_con(i,j) + prs_gsp(i,j)
        ENDIF

        ! infiltration and surface run-off

        ! ice free fraction of first soil layer scaled by pore volume
        ! is used as reduction factor for infiltration rate
        zfr_ice_free     = 1._ireals-ziw_fr(i,j,1)/zporv(i,j)

        ! subtract evaporation from interception store to avoid negative
        ! values due to sum of evaporation+infiltration
        zwinstr = zwin(i,j) + zdwidt(i,j)*zdtdrhw
        zwinstr = MAX(0.0_ireals,zwinstr)

        ! maximum infiltration rate of the soil (rock/ice/water-exclusion
        ! LINDA, b
        IF (LDECHARME) THEN
           zinfmx = zrock(i,j)*zkw(i,j,1)*rho_w
        ELSE
           zinfmx = zrock(i,j)*zfr_ice_free*csvoro &
                    *( cik1*MAX(0.5_ireals,plcov(i,j))*MAX(0.0_ireals,           &
                    zporv(i,j)-zw_fr(i,j,1))/zporv(i,j) + zik2(i,j) )
        END IF
        ! LINDA, e
        ! to avoid pore volume water excess of the uppermost layer by 
        ! infiltration
        zinfmx = MIN(zinfmx, (zporv(i,j) - zw_fr(i,j,1))*zdzhs(1)*zrhwddt)

        ! to avoid infiltration at snow covered parts of soil surface
        zinfmx = zinfmx*(1._ireals - zf_snow(i,j))

        zwimax = cwimax_ml*(1._ireals + plcov(i,j)*5._ireals)
        zalf   = SQRT(MAX(0.0_ireals,1.0_ireals - zwinstr/zwimax))

        ! water supply from interception store (if Ts above freezing)
        zinf   = zfr_ice_free*zwinstr*rho_w/ctau_i

        ! possible contribution of rain to infiltration
        IF (zrr(i,j)-zepsi > 0.0_ireals) THEN
          zalf = MAX( zalf,                                                   &
                 (zrhwddt*MAX(0.0_ireals, zwimax-zwinstr) + zinf)/zrr(i,j) )
          zalf = MAX( 0.01_ireals, MIN(1.0_ireals, zalf) )

          ! if rain falls onto snow, all rain is considered for infiltration
          ! as no liquid water store is considered in the snowpack
          IF (zwsnow(i,j) > 0.0_ireals) zalf = 0.0_ireals
        END IF
        ! Increase infiltration and reduce surface runoff (bugfix)
        zalf = 0.0_ireals
        ! rain freezes on the snow surface
        IF (lmulti_snow .AND. zwsnow(i,j) > 0.0_ireals) zalf = 1.0_ireals
        ! add rain contribution to water supply for infiltration
        zvers = zinf + (1._ireals - zalf)*zrr(i,j)
        ! final infiltration rate limited by maximum value
        zinfil(i,j) = MIN(zinfmx,zvers)

        ! surface run-off (residual of potential minus actual infiltration)
        zro_inf       = zvers - zinfil(i,j)
        runoff_s(i,j) = runoff_s(i,j) + zro_inf*zroffdt

        ! change of snow water and interception water store
        ! (negligible residuals are added to the run-off)

        ! snow store
        zdwsndtt = zrs(i,j) + zdwsndt(i,j)
        zwsnstr  = zwsnow(i,j) + zdwsndtt*zdtdrhw
        zwsnstr  = MAX(0.0_ireals, zwsnstr) ! avoid negative values (security)
        zdwseps  = 0.0_ireals
        IF (zwsnstr > 0.0_ireals .AND. zwsnstr < zepsi) THEN
!         IF (ztsnow_pm(i,j) > 0.0_ireals) THEN
            zdwseps    = zwsnstr*zrhwddt
            runoff_s(i,j) = runoff_s(i,j) + zdwseps*zroffdt
            zdwsndtt   = zdwsndtt - zdwseps
!         END IF
        END IF
        zdwsndt(i,j) = zdwsndtt

        ! interception store
        zdwidtt  = zalf*zrr(i,j) + zdwidt(i,j)-zinf 
        zwinstr  = zwin(i,j) + zdwidtt*zdtdrhw
        zwinstr  = MAX(0.0_ireals, zwinstr) !avoid negative values (security)
        zdwieps  = 0.0_ireals
        IF (zwinstr > 0.0_ireals .AND. zwinstr < zepsi) THEN
          zdwieps    = zwinstr*zrhwddt
          runoff_s(i,j)= runoff_s(i,j) + zdwieps*zroffdt
          zdwidtt    = zdwidtt - zdwieps
          zwinstr    = 0.0_ireals
        END IF
        ! add excess over zwimax to runoff
        zro_wi       = zrhwddt*MAX( 0.0_ireals, zwinstr-zwimax )
        zdwidtt      = zdwidtt - zro_wi
        zdwidt(i,j)  = zdwidtt
        runoff_s(i,j)= runoff_s(i,j) + zro_wi*zroffdt
      END IF            ! land-points only
    END DO
  END DO

!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------

! First loop through all layers from top to bottom, calculate diffusivity
! and first-guess values for the conductivity. These will then be limited
! in the second loop that goes from bottom to top.

! uppermost layer, kso = 1
DO   j = jstarts, jends
  DO i = istarts, iends
    ! sedimentation and capillary transport in soil
    ! Note: The fractional liquid water content (concentration)  of each layer
    !       is normalized by the ice free fraction of each layer in order to
    !       obtain a representative concentration of liquid water in the
    !       'active' part of each soil layer
    !       Hydraulic diffusivity and conductivity coefficients are multiplied
    !       by a reduction factor depending on the maximum ice fraction of the
    !       adjacent layers in order to avoid the transport of liquid water
    !       in to the frozen part of the adjacent layer
    IF (llandmask(i,j)) THEN      ! land-points only
      IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
        zice_fr_kso(i,j,1)   = ziw_fr(i,j,1)
        zice_fr_kso(i,j,2)   = ziw_fr(i,j,2)
        zlw_fr_kso(i,j,1)    = zlw_fr(i,j,1)
        zlw_fr_kso(i,j,2)    = zlw_fr(i,j,2)

        ! compute reduction factor for transport coefficients
        zfr_ice          = max (zice_fr_kso(i,j,1),zice_fr_kso(i,j,2))
        zredp05          = 1._ireals-zfr_ice/MAX(zlw_fr(i,j,1)+zice_fr_kso(i,j,1),zlw_fr(i,j,2)+zice_fr_kso(i,j,2))

        ! interpolated scaled liquid water fraction at layer interface
        zlw_fr_ksop05 = 0.5_ireals*(zdzhs(2)*zlw_fr(i,j,1)+zdzhs(1)*zlw_fr(i,j,2)) &
                                             /zdzms(2)
        zdlw_fr_ksop05_fg(i,j,1) = zredp05*zdw(i,j)*EXP(zdw1(i,j)*                       &
                         (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
        zklw_fr_ksop05_fg(i,j,1) = zredp05*zkw(i,j,1)*EXP(zkw1(i,j)*                       &
                         (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )


      ENDIF
    ENDIF
  END DO
END DO

! inner layers 2 <=kso<=ke_soil_hy-1
DO   kso =2,ke_soil_hy-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      ! sedimentation and capillary transport in soil
      IF (llandmask(i,j)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zice_fr_kso(i,j,kso)   = ziw_fr(i,j,kso  )
          zice_fr_kso(i,j,kso+1) = ziw_fr(i,j,kso+1)
          zlw_fr_kso(i,j,kso)    = zlw_fr(i,j,kso  )
          zlw_fr_kso(i,j,kso+1)  = zlw_fr(i,j,kso+1)
          ! interpolated scaled liquid water content at interface to layer
          ! above and below
          !zlw_fr_ksom05 = 0.5*(zdzhs(kso-1)*zlw_fr_kso(i,j,kso)+   &
          !                       zdzhs(kso)*zlw_fr_kso(i,j,kso-1))/zdzms(kso)
          zlw_fr_ksop05 = 0.5*(zdzhs(kso+1)*zlw_fr(i,j,kso)+   &
                                 zdzhs(kso)*zlw_fr(i,j,kso+1))/zdzms(kso+1)

          ! compute reduction factor for coefficients
          zfr_ice          = max (zice_fr_kso(i,j,kso),zice_fr_kso(i,j,kso+1))
          zredp05 = 1._ireals-zfr_ice/max (zlw_fr_kso(i,j,kso)+zice_fr_kso(i,j,kso),zlw_fr(i,j,kso+1)+zice_fr_kso(i,j,kso+1))
          ! first guess K
          if ((kso+1)==zsatlev(i,j)) then ! last level above saturation, partially saturated

             zklw_fr_ksop05_fg(i,j,kso) = zredp05*zkw(i,j,kso)*EXP( zkw1(i,j)  &
                                          *(zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
             zdlw_fr_ksop05_fg(i,j,kso) = zredp05*zdw(i,j)*EXP( zdw1(i,j)   &
                                         *(zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
             q_kso(i,j,kso) = zredp05*s_oro(i,j)*zgamma(i,j)*rho_w*zkw(i,j,kso)*h_wt(i,j) ! saturated ground water flow at the bottom

          else if (kso.ge.zsatlev(i,j)) then ! entirely saturated

             zklw_fr_ksop05_fg(i,j,kso) = zredp05*zkw(i,j,kso)*EXP( zkw1(i,j)   &
                                         *(zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )!&
             q_kso(i,j,kso) = zredp05*s_oro(i,j)*zgamma(i,j)*rho_w*zkw(i,j,kso)*zdzhs(kso)

          else ! default case, unsaturated

             zklw_fr_ksop05_fg(i,j,kso) = zredp05*zkw(i,j,kso)*EXP( zkw1(i,j)   &
                                         *(zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )!&
             zdlw_fr_ksop05_fg(i,j,kso) = zredp05*zdw(i,j)*EXP( zdw1(i,j)   &
                                         *(zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
             q_kso(i,j,kso) = 0.0_ireals
          end if
        ENDIF
      ENDIF
    END DO
  END DO
END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          ! lowest active hydrological layer ke_soil_hy-1
          zice_fr_kso(i,j,ke_soil_hy)   = ziw_fr(i,j,ke_soil_hy  )
          zlw_fr_kso(i,j,ke_soil_hy)    = zlw_fr(i,j,ke_soil_hy  )
          !zlw_fr_ksom05 = 0.5*(zdzhs(ke_soil_hy-1)*zlw_fr_kso(i,j,ke_soil_hy)+ &
          !                     zdzhs(ke_soil_hy)*zlw_fr_kso(i,j,ke_soil_hy-1))/zdzms(ke_soil_hy)

          zfr_ice          = zice_fr_kso(i,j,ke_soil_hy)
          zred = 1._ireals-zfr_ice/(zlw_fr(i,j,ke_soil_hy)+zice_fr_kso(i,j,ke_soil_hy))

          zklw_fr_ksop05_fg(i,j,ke_soil_hy) = 0.0_ireals
          zdlw_fr_ksop05_fg(i,j,ke_soil_hy) = 0.0_ireals

          if (ke_soil_hy.ge.zsatlev(i,j)) then ! fully saturated

             q_kso(i,j,ke_soil_hy) = zred*s_oro(i,j)*zgamma(i,j)*rho_w*zkw(i,j,ke_soil_hy)*zdzhs(ke_soil_hy)

          else ! partially saturated

             q_kso(i,j,ke_soil_hy) = zred*s_oro(i,j)*zgamma(i,j)*rho_w*zkw(i,j,ke_soil_hy)*h_wt(i,j)

          end if
        ENDIF
      ENDIF
    END DO
  END DO

! second loop, go from bottom to top and limit the fluxes by adapting K.
! then calculate the coefficients for the matrix with these limited
! fluxes.

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          ! lowest active hydrological layer ke_soil_hy-1

          ! limit first guess K such that no oversaturation can occur in this layer
          zklw_fr_ksom05_limit =(zporv(i,j)-zlw_fr(i,j,ke_soil_hy)-zice_fr_kso(i,j,ke_soil_hy  ))/zdt*zdzhs(ke_soil_hy) &
                               + zdlw_fr_ksop05_fg(i,j,ke_soil_hy-1) &
                               *(zlw_fr(i,j,ke_soil_hy)+zice_fr_kso(i,j,ke_soil_hy) - zlw_fr(i,j,ke_soil_hy-1)- zice_fr_kso(i,j,ke_soil_hy-1))/zdzms(ke_soil_hy)

          zklw_fr_ksop05_sg(i,j,ke_soil_hy-1) = min(zklw_fr_ksop05_fg(i,j,ke_soil_hy-1),zklw_fr_ksom05_limit)

          if(zklw_fr_ksop05_sg(i,j,ke_soil_hy-1)==zklw_fr_ksom05_limit) then
             limit=.true.
             !print*,'limit'
          end if

          ! limit this second guess fluxes such that no depletion can occur
          flux_ksom05     = (zdlw_fr_ksop05_fg(i,j,ke_soil_hy-1)*(zlw_fr(i,j,ke_soil_hy)+zice_fr_kso(i,j,ke_soil_hy)-zlw_fr(i,j,ke_soil_hy-1)-zice_fr_kso(i,j,ke_soil_hy))/  &
                             zdzms(ke_soil_hy) - zklw_fr_ksop05_sg(i,j,ke_soil_hy-1))/zdzhs(ke_soil_hy)
          flux_ksom05_pos = max(0.0_ireals,flux_ksom05)
          q_kso_pos       = max(0.0_ireals,q_kso(i,j,ke_soil_hy)/rho_w)
!          alpha_lim       = max(0.0_ireals,min(1.0_ireals,&
!                                zlw_fr(i,j,ke_soil_hy)/(zdt*(flux_ksom05_pos+q_kso_pos+zepsi))))
          alpha_lim       = max(0.0_ireals,min(1.0_ireals,&
                                zlw_fr(i,j,ke_soil_hy)/(zdt*(flux_ksom05_pos+zepsi))))
          zklw_fr_ksop05(i,j,ke_soil_hy-1)    = alpha_lim * zklw_fr_ksop05_sg(i,j,ke_soil_hy-1)
          zdlw_fr_ksop05(i,j,ke_soil_hy-1)    = alpha_lim * zdlw_fr_ksop05_fg(i,j,ke_soil_hy-1)
          q_kso(i,j,ke_soil_hy)               = alpha_lim * q_kso(i,j,ke_soil_hy)
              
          if(alpha_lim.lt.1.0_ireals) then
             print*,'limit alpha',ke_soil_hy, zlw_fr(i,j,ke_soil_hy),zklw_fr_ksop05_fg(i,j,ke_soil_hy-1),zklw_fr_ksop05_sg(i,j,ke_soil_hy-1), alpha_lim, zdlw_fr_ksop05_fg(i,j,ke_soil_hy-1),q_kso(i,j,ke_soil_hy)
          end if


          z1dgam1 = zdt/zdzhs(ke_soil_hy)
          zgam2m05  = zdlw_fr_ksop05(i,j,ke_soil_hy-1)/zdzms(ke_soil_hy)

          zaga(i,j,ke_soil_hy) = -zalfa* zgam2m05*z1dgam1
          zagb(i,j,ke_soil_hy) = 1.+ zalfa*zgam2m05*z1dgam1
          zagc(i,j,ke_soil_hy) = 0.0
          zagd(i,j,ke_soil_hy) = zlw_fr(i,j,ke_soil_hy)+z1dgam1*zklw_fr_ksop05(i,j,ke_soil_hy-1) &
                               + (1.-zalfa)*z1dgam1*                              &
                                 zgam2m05*(zlw_fr(i,j,ke_soil_hy-1)  - zlw_fr(i,j,ke_soil_hy)) &
                               + z1dgam1*                                         &
                                 zgam2m05*(zice_fr_kso(i,j,ke_soil_hy-1) - zice_fr_kso(i,j,ke_soil_hy))   
          zflmg(i,j,ke_soil_hy)= 0.0
        ENDIF
      ENDIF
    END DO
  END DO

! inner layers 2 <=kso<=ke_soil_hy-1
DO kso = ke_soil_hy-1,2,-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      ! sedimentation and capillary transport in soil
      IF (llandmask(i,j)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type

           ! limit first guess K such that no oversaturation can occur in this layer
           zklw_fr_ksom05_limit = (zporv(i,j)-zlw_fr(i,j,kso)-zice_fr_kso(i,j,kso))/zdt*zdzhs(kso) &
                                 + zdlw_fr_ksop05_fg(i,j,kso-1)*(zlw_fr(i,j,kso)   +zice_fr_kso(i,j,kso)  - zlw_fr(i,j,kso-1)-zice_fr_kso(i,j,kso-1))/zdzms(kso)&
                                 - zdlw_fr_ksop05_fg(i,j,kso)  *(zlw_fr(i,j,kso+1) +zice_fr_kso(i,j,kso+1)- zlw_fr(i,j,kso)  -zice_fr_kso(i,j,kso-1))/zdzms(kso+1)
           zklw_fr_ksop05_sg(i,j,kso-1) =  min(zklw_fr_ksop05_fg(i,j,kso-1),zklw_fr_ksom05_limit)
           if(zklw_fr_ksop05_sg(i,j,kso-1)==zklw_fr_ksom05_limit) then
              limit=.true.
             !print*,'limit',kso, zklw_fr_ksop05_sg(i,j,kso-1),zklw_fr_ksom05_limit, zklw_fr_ksop05(i,j,kso-1)
           end if

           ! limit this second guess fluxes such that no depletion can occur
           flux_ksom05 = (zdlw_fr_ksop05_fg(i,j,kso-1)* &
                         (zlw_fr(i,j,kso)+zice_fr_kso(i,j,kso)-zlw_fr(i,j,kso-1)-zice_fr_kso(i,j,kso-1))/  &
                          zdzms(kso) - zklw_fr_ksop05_sg(i,j,kso-1))/zdzhs(kso)
           flux_ksop05 = (zdlw_fr_ksop05_fg(i,j,kso)* &
                         (zlw_fr(i,j,kso+1)+zice_fr_kso(i,j,kso+1)-zlw_fr(i,j,kso)-zice_fr_kso(i,j,kso))/  &
                          zdzms(kso+1) - zklw_fr_ksop05_sg(i,j,kso))/zdzhs(kso)
           flux_ksom05_pos = max(0.0_ireals,flux_ksom05)
           flux_ksop05_neg = min(0.0_ireals,flux_ksop05)
           q_kso_pos       = max(0.0_ireals,q_kso(i,j,kso)/rho_w)
           alpha_lim       = max(0.0_ireals,min(1.0_ireals,&
                                 zlw_fr(i,j,kso)/(zdt*(-flux_ksop05_neg+flux_ksom05_pos+q_kso_pos+zepsi))))
           zklw_fr_ksop05(i,j,kso-1)    = alpha_lim * zklw_fr_ksop05_sg(i,j,kso-1)
           zdlw_fr_ksop05(i,j,kso-1)    = alpha_lim * zdlw_fr_ksop05_fg(i,j,kso-1)
           zklw_fr_ksop05(i,j,kso)      = min(alpha_lim * zklw_fr_ksop05_sg(i,j,kso),zklw_fr_ksop05_fg(i,j,kso))
           zdlw_fr_ksop05(i,j,kso)      = min(alpha_lim * zdlw_fr_ksop05_fg(i,j,kso),zdlw_fr_ksop05(i,j,kso))
           q_kso(i,j,kso)               = alpha_lim * q_kso(i,j,kso)
           if(alpha_lim.lt.1.0_ireals) then
              !print*,'limit alpha',kso, zklw_fr_ksop05_fg(i,j,kso-1),zklw_fr_ksop05_sg(i,j,kso-1), alpha_lim
              print*,'limit alpha',kso, i, j, kso, my_cart_id, flux_ksom05_pos, flux_ksop05_neg, q_kso_pos, alpha_lim
           end if

           ! coefficients for implicit flux computation
           z1dgam1 = zdt/zdzhs(kso)
           zgam2m05  = zdlw_fr_ksop05(i,j,kso-1)/zdzms(kso)
           zgam2p05  = zdlw_fr_ksop05(i,j,kso)/zdzms(kso+1)

           zaga (i,j,kso) = -zalfa*zgam2m05*z1dgam1
           zagc (i,j,kso) = -zalfa*zgam2p05*z1dgam1
           zagb (i,j,kso) = 1._ireals +zalfa*(zgam2m05+zgam2p05)*z1dgam1
           zagd (i,j,kso) = zlw_fr(i,j,kso)+                               &
                            z1dgam1*(-zklw_fr_ksop05(i,j,kso)+zklw_fr_ksop05(i,j,kso-1))+ &
                            (1._ireals-zalfa)*z1dgam1*                &
                            (zgam2p05*(zlw_fr(i,j,kso+1)-zlw_fr(i,j,kso)  )     &
                            -zgam2m05*(zlw_fr(i,j,kso)  -zlw_fr(i,j,kso-1)))    &
                            +z1dgam1*                                  &
                            (zgam2p05*(zice_fr_kso(i,j,kso+1)-zice_fr_kso(i,j,kso)  )   &
                            -zgam2m05*(zice_fr_kso(i,j,kso)-zice_fr_kso(i,j,kso-1))   )

           !soil water flux, explicit part, for soil water flux investigations
           ! only)
           zflmg(i,j,kso) = rho_w*(zdlw_fr_ksop05(i,j,kso-1)*  &
                           (zlw_fr(i,j,kso)+zice_fr_kso(i,j,kso)-zlw_fr(i,j,kso-1)-zice_fr_kso(i,j,kso-1))/  &
                            zdzms(kso) - zklw_fr_ksop05(i,j,kso-1))

        ENDIF
      ENDIF
    END DO
  END DO
END DO

! uppermost layer, kso = 1
DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN      ! land-points only
      IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type


        ! coefficients for implicit flux computation
        z1dgam1 = zdt/zdzhs(1) 
        zgam2p05 = zdlw_fr_ksop05(i,j,1)/zdzms(2)
        zaga(i,j,1) = 0._ireals
        zagb(i,j,1) = 1._ireals+zalfa*zgam2p05*z1dgam1
        zagc(i,j,1) = -zalfa * zgam2p05*z1dgam1
        zagd(i,j,1) = zlw_fr(i,j,1) + zinfil(i,j)*z1dgam1/rho_w  &
                     -zklw_fr_ksop05(i,j,1)*z1dgam1                     &
                     +(1. - zalfa)* zgam2p05*z1dgam1*(zlw_fr(i,j,2) - zlw_fr(i,j,1))  &
                     +              zgam2p05*z1dgam1*(zice_fr_kso(i,j,2)-zice_fr_kso(i,j,1))

        ! explicit part of soil surface water flux:
        zflmg (i,j,1) = - zinfil(i,j)! boundary value for soil water transport
      ENDIF
    ENDIF
  END DO
END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        ! generalized upper boundary condition
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zagc(i,j,1) = zagc(i,j,1)/zagb(i,j,1)
          zagd(i,j,1) = zagd(i,j,1)/zagb(i,j,1)
        ENDIF
      END IF          ! land-points only
    END DO
  END DO

DO kso=2,ke_soil_hy-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zzz = 1./(zagb(i,j,kso) - zaga(i,j,kso)*zagc(i,j,kso-1))
          zagc(i,j,kso) = zagc(i,j,kso) * zzz
          zagd(i,j,kso) = (zagd(i,j,kso) - zaga(i,j,kso)*zagd(i,j,kso-1)) * zzz
        ENDIF
      END IF          ! land-points only
    END DO
  END DO
END DO                ! soil layers

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
           zage(i,j,ke_soil_hy) = (zagd(i,j,ke_soil_hy)-zaga(i,j,ke_soil_hy)*  &
                             zagd(i,j,ke_soil_hy-1))/                          &
                            (zagb(i,j,ke_soil_hy) - zaga(i,j,ke_soil_hy)*      &
                             zagc(i,j,ke_soil_hy-1))
        ENDIF
      END IF          ! land-points only
   END DO
 END DO

 DO kso = ke_soil_hy-1,1,-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zage(i,j,kso)     = zagd(i,j,kso) - zagc(i,j,kso)*zage(i,j,kso+1)
          ! compute implicit part of new liquid water content and add existing
          ! ice content
          w_so(i,j,kso,nnew) = zage(i,j,kso)*zdzhs(kso) + w_so_ice(i,j,kso,nx)
        END IF
      END IF          ! land-points only
    END DO
  END DO
END DO                ! soil layers

!lowest active hydrological level
 DO   j = jstarts, jends
   DO i = istarts, iends
     IF (llandmask(i,j)) THEN          ! land-points only
       IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
         ! boundary values ensure that the calculation below leaves the climate
         ! layer water contents unchanged compute implicit part of new liquid
         ! water content and add existing ice content
         w_so(i,j,ke_soil_hy,nnew) = zage(i,j,ke_soil_hy)*zdzhs(ke_soil_hy) + &
                         w_so_ice(i,j,ke_soil_hy,nx)
       END IF 
     END IF          ! land-points only
   END DO
 END DO

! to ensure vertical constant water concentration profile beginning at 
! layer ke_soil_hy for energetic treatment only
! soil water climate layer(s)

IF (itype_hydbound == 3) THEN
  ! ground water as lower boundary of soil column
  DO kso = ke_soil_hy+1,ke_soil+1
   DO   j = jstarts, jends
     DO i = istarts, iends
       IF (llandmask(i,j)) THEN          ! land-points only
         IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
           w_so(i,j,kso,nnew) = zporv(i,j)*zdzhs(kso)
         END IF
       END IF          ! land-points only
     END DO
   END DO
  END DO
ELSE
  DO kso = ke_soil_hy+1,ke_soil+1
   DO   j = jstarts, jends
     DO i = istarts, iends
       IF (llandmask(i,j)) THEN          ! land-points only
         IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
           w_so(i,j,kso,nnew) = w_so(i,j,kso-1,nnew)*zdzhs(kso)/zdzhs(kso-1)
         END IF
       END IF          ! land-points only
     END DO
   END DO
  END DO
ENDIF

DO  kso = 1,ke_soil
   DO   j = jstarts, jends
      ! utility variables used to avoid if-constructs in following loops
      zro_sfak = zsf_heav(0.5_ireals + msr_off - kso)  ! 1.0 for 'surface runoff'
      zro_gfak = 1._ireals - zro_sfak                  ! 1.0 for 'ground runoff'
      IF (kso==ibot_w_so) THEN
         zfmb_fak = 1.0_ireals
      ELSE
         zfmb_fak = 0.0_ireals
      END IF
      DO i = istarts, iends
         IF (llandmask(i,j)) THEN      ! land-points only
          ! sedimentation and capillary transport in soil
            IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type

               zdwg =  (w_so(i,j,kso,nnew)/zdzhs(kso)-zw_fr(i,j,kso))*zdzhs(kso) &
                                                                   /zdtdrhw

               ! add evaporation (znlgw1f: first layer only)
               ! and transpiration (for each layer)
               zdwg   = zdwg + znlgw1f(kso) * zesoil(i,j) + ztrang (i,j,kso)
               zwgn   = zw_fr(i,j,kso) + zdtdrhw*zdwg/zdzhs(kso)
               zro2   = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zwgn - zporv(i,j))
               zkorr  = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zadp(i,j) - zwgn )
               zdwgdt(i,j,kso)= zdwg + zkorr - zro2
               zro    = zro2
               runoff_s(i,j) = runoff_s(i,j) + zro*zro_sfak*zroffdt

               klw(i,j,kso) = zklw_fr_ksop05(i,j,kso)
               dlw(i,j,kso) = zdlw_fr_ksop05(i,j,kso)
               flx(i,j,kso) = zflmg(i,j,kso)
            END IF
         END IF  ! land-points only
      END DO
   END DO
END DO        ! soil layers

DO  kso = 1,ke_soil_hy
   DO   j = jstarts, jends
      DO i = istarts, iends
         IF (llandmask(i,j)) THEN      ! land-points only
          ! sedimentation and capillary transport in soil
            IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
               ! remove runoff
               zdwgdt(i,j,kso) = zdwgdt(i,j,kso) - q_kso(i,j,kso)*zdzhs(kso)
               runoff_g(i,j) = runoff_g(i,j)- q_kso(i,j,kso)*zdzhs(kso)*zdt
            END IF
         END IF  ! land-points only
      END DO
   END DO
END DO        ! soil layers
DO  kso = 1,ke_soil+1
   DO   j = jstarts, jends
      DO i = istarts, iends
         IF (llandmask(i,j)) THEN      ! land-points only
          ! sedimentation and capillary transport in soil
            IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
               ! remove runoff
               q_roff(i,j,kso) = q_kso(i,j,kso)
            END IF
         END IF  ! land-points only
      END DO
   END DO
END DO        ! soil layers

!------------------------------------------------------------------------------
! Section II.5: Soil surface heat flux (thermal forcing)
!------------------------------------------------------------------------------

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        ! Estimate thermal surface fluxes:
        ! Estimate thermal surface fluxes over snow covered and snow free
        ! part of surface based on area mean values calculated in radiation
        ! code (positive = downward)

        ! Correction of the long-wave radiation flux due to low frequency
        ! of the radiation routine calls (so far for the multi-layer model only)
        IF ((ntstep >= 1) .OR. (nincrad == 1)) THEN
          IF (hincrad /= 0.0) THEN
            zradstep = hnextrad - hincrad
            zradstep = NINT (zradstep * 3600.0_ireals / dt)
          ELSE
            zradstep = nextrad - nincrad
          ENDIF
        ENDIF

        IF ( (ntstep < 2) .OR. (ntstep == zradstep) ) THEN
          tg_radstep (i,j) = t_g(i,j,nx)
        END IF
        thbs_new(i,j) = thbs(i,j) + sigma*(1._ireals - Ctalb)* &
                    (tg_radstep(i,j)**4 - t_g(i,j,nx)**4  )

        IF(lmulti_snow) THEN
          zgstr = sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i,j))* &
                  zts(i,j) + zf_snow(i,j)*ztsnow(i,j) )**4 + thbs_new(i,j)
        ELSE
          zgstr = sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i,j))* &
                  zts(i,j) + zf_snow(i,j)*ztsnow(i,j) )**4 + thbs(i,j)
        END IF
        zthsnw(i,j) = - sigma*(1._ireals - Ctalb)*ztsnow(i,j)**4 + zgstr
        zthsoi(i,j) = - sigma*(1._ireals - Ctalb)*zts(i,j)**4 + zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
        !
        ! Forcing for snow-free soil:
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux
        zrnet_s(i,j) = (1._ireals - zf_snow(i,j))*(sobs(i,j)+zthsoi(i,j))

        zshfl_s(i,j) = (1._ireals - zf_snow(i,j))*cp_d*zrhoch(i,j)*  &
                                                      (zth_low(i,j) - zts(i,j))
        zlhfl_s(i,j) = (zts_pm(i,j)*lh_v + (1._ireals-zts_pm(i,j))*lh_s)*zverbo(i,j)
      END IF          ! land-points only
    END DO
  END DO

IF (lmulti_snow) THEN
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        zsprs(i,j)   = 0.0_ireals
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i,j)*zrs(i,j) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs(i,j)   = - lh_f*zrs(i,j)
          zdwidt (i,j) = zdwidt (i,j) + zrs(i,j)
          zdwsndt(i,j) = zdwsndt(i,j) - zrs(i,j)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
          zwimax       = cwimax_ml*(1._ireals + plcov(i,j)*5._ireals)
          zwinstr      = zwin(i,j) + zdwidt(i,j)*zdtdrhw
          IF (zwinstr > zwimax) THEN  ! overflow of interception store
            zro        = (zwinstr - zwimax)*zrhwddt
            zdwidt(i,j)= zdwidt(i,j) - zro
            runoff_s(i,j)  = runoff_s(i,j) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF ((1._ireals-ztsnow_pm(i,j))*zrr(i,j) > 0.0_ireals) THEN
          zsprs(i,j)   = lh_f*zrr(i,j)
          zdwidt (i,j) = zdwidt (i,j) - zrr(i,j)
          zdwsndt(i,j) = zdwsndt(i,j) + zrr(i,j)
        END IF

!       Influence of heatflux through snow on total forcing:
        zwsnew(i,j)   = zwsnow(i,j) + zdwsndt(i,j)*zdtdrhw
        IF (zwsnew(i,j).GT.zepsi) THEN

          zrho_snowf = crhosminf+(crhosmaxf-crhosminf)* (zth_low(i,j)-csnow_tmin) &
                                                  /(t0_melt          -csnow_tmin)
          zrho_snowf = MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))

          IF(zdwsndt(i,j)-zrs(i,j)-zrr(i,j).GT.0.0_ireals) THEN

            wtot_snow(i,j,1,nx) = max(wtot_snow(i,j,1,nx) + zdwsndt(i,j)*zdtdrhw,0.0_ireals)

            zhm_snow(i,j,1) = zhm_snow(i,j,1) - (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdt/rho_i/2.-  &
              zrs(i,j)*zdt/zrho_snowf/2.- zrr(i,j)*zdt/rho_i/2.
            zdzh_snow(i,j,1) = zdzh_snow(i,j,1) + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdt/rho_i +  &
              zrs(i,j)*zdt/zrho_snowf + zrr(i,j)*zdt/rho_i

            rho_snow_mult(i,j,1,nx) = max(wtot_snow(i,j,1,nx)*rho_w/zdzh_snow(i,j,1),0.0_ireals)
          ELSE

            wtot_snow(i,j,1,nx) = max(wtot_snow(i,j,1,nx) + (zrs(i,j)+zrr(i,j))*zdtdrhw,0.0_ireals)

            zhm_snow(i,j,1)  = zhm_snow(i,j,1) - zrs(i,j)*zdt/zrho_snowf/2.- zrr(i,j)*zdt/rho_i/2.
            zdzh_snow(i,j,1) = zdzh_snow(i,j,1) + zrs(i,j)*zdt/zrho_snowf + zrr(i,j)*zdt/rho_i

            IF(wtot_snow(i,j,1,nx) .GT. 0._ireals) THEN
              rho_snow_mult(i,j,1,nx) = MAX(wtot_snow(i,j,1,nx)*rho_w/zdzh_snow(i,j,1),0.0_ireals)

                wtot_snow(i,j,1,nx) = MAX(wtot_snow(i,j,1,nx) + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdtdrhw, &
                                            0.0_ireals)

              zhm_snow(i,j,1)  = zhm_snow(i,j,1) - (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdt/rho_snow_mult(i,j,1,nx)/2.
              zdzh_snow(i,j,1) = zdzh_snow(i,j,1) + (zdwsndt(i,j)-zrs(i,j)-zrr(i,j))*zdt/rho_snow_mult(i,j,1,nx)
            ELSE
              rho_snow_mult(i,j,1,nx) = 0.0_ireals
              zdzh_snow(i,j,1) = 0.0_ireals
            END IF

          END IF
        END IF
        h_snow(i,j,nnow) = 0.
        sum_weight(i,j) = 0.0_ireals
      END IF          ! land-points only
    END DO
  END DO
  DO ksn = 1,ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF (zwsnew(i,j).GT.zepsi) THEN
            h_snow(i,j,nnow) = h_snow(i,j,nnow) + zdzh_snow(i,j,ksn)
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO

  DO ksn = ke_snow,1,-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN  ! for landpoints only
          IF(zwsnew(i,j) .GT. zepsi) THEN
            IF(zwsnow(i,j) .GT. zepsi) THEN
              dz_old(i,j,ksn) = zdzh_snow(i,j,ksn)
              z_old(i,j,ksn) = -sum_weight(i,j) - zdzh_snow(i,j,ksn)/2._ireals
              sum_weight(i,j) = sum_weight(i,j) + zdzh_snow(i,j,ksn)
              zhh_snow(i,j,ksn) = -h_snow(i,j,nnow)/ke_snow*(ke_snow-ksn)
            ELSE
              zhh_snow(i,j,ksn) = -h_snow(i,j,nnow)/ke_snow*(ke_snow-ksn)
            END IF
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        IF(zwsnew(i,j) .GT. zepsi) THEN
          IF(zwsnow(i,j) .GT. zepsi) THEN
            zhm_snow (i,j,1) = (-h_snow(i,j,nnow) + zhh_snow(i,j,1))/2._ireals
            !layer thickness betw. half levels of uppermost snow layer
            zdzh_snow(i,j,1) = zhh_snow(i,j,1) + h_snow(i,j,nnow)            
            !layer thickness between snow surface and main level of uppermost layer
            zdzm_snow(i,j,1) = zhm_snow(i,j,1) + h_snow(i,j,nnow)            
            IF(dz_old(i,j,1).ne.0..and.rho_snow_mult(i,j,1,nnow).ne.0.) THEN
              wliq_snow(i,j,1,nnow) = wliq_snow(i,j,1,nnow)/dz_old(i,j,1)
            END IF 
          ELSE
            zhm_snow (i,j,1) = (-h_snow(i,j,nnow) + zhh_snow(i,j,1))/2._ireals
            zdzh_snow(i,j,1) = zhh_snow(i,j,1) + h_snow(i,j,nnow)
            zdzm_snow(i,j,1) = zhm_snow(i,j,1) + h_snow(i,j,nnow)
          END IF
        END IF
      END IF          ! land-points only
    END DO
  END DO
  DO ksn = 2,ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN  ! for landpoints only
          IF(zwsnew(i,j) .GT. zepsi) THEN
            IF(zwsnow(i,j) .GT. zepsi) THEN
              zhm_snow  (i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals
              zdzh_snow (i,j,ksn) = zhh_snow(i,j,ksn) - zhh_snow(i,j,ksn-1) ! layer thickness betw. half levels
              zdzm_snow(i,j,ksn ) = zhm_snow(i,j,ksn) - zhm_snow(i,j,ksn-1) ! layer thickness betw. main levels
              IF(dz_old(i,j,ksn).ne.0..and.rho_snow_mult(i,j,ksn,nnow).ne.0.) THEN
                wliq_snow(i,j,ksn,nnow) = wliq_snow(i,j,ksn,nnow)/dz_old(i,j,ksn)
              END IF
            ELSE 
              zhm_snow (i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals
              zdzh_snow(i,j,ksn) = zhh_snow(i,j,ksn) - zhh_snow(i,j,ksn-1) ! layer thickness betw. half levels
              zdzm_snow(i,j,ksn) = zhm_snow(i,j,ksn) - zhm_snow(i,j,ksn-1) ! layer thickness betw. main levels
            END IF
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO

  DO ksn = ke_snow,1,-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        t_new  (i,j,ksn) = 0.0_ireals
        rho_new(i,j,ksn) = 0.0_ireals
        wl_new (i,j,ksn) = 0.0_ireals
      END DO
    END DO

    DO k = ke_snow,1,-1
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN  ! for landpoints only
            IF(zwsnew(i,j) .GT. zepsi .AND. zwsnow(i,j) .GT. zepsi) THEN

              weight = MAX(MIN(z_old(i,j,k)+dz_old(i,j,k)/2._ireals,zhm_snow(i,j,ksn)+zdzh_snow(i,j,ksn)/2._ireals)-   &
                       MAX(z_old(i,j,k)-dz_old(i,j,k)/2._ireals, &
                       zhm_snow(i,j,ksn)-zdzh_snow(i,j,ksn)/2._ireals),0._ireals)/zdzh_snow(i,j,ksn)

              t_new  (i,j,ksn) = t_new  (i,j,ksn) + ztsnow_mult  (i,j,k   )*weight
              rho_new(i,j,ksn) = rho_new(i,j,ksn) + rho_snow_mult(i,j,k,nx)*weight
              wl_new (i,j,ksn) = wl_new (i,j,ksn) + wliq_snow    (i,j,k,nx)*weight
            END IF
          END IF          ! land-points only
        END DO
      END DO
    END DO
  END DO
  DO ksn = ke_snow,1,-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN  ! for landpoints only
          IF(zwsnew(i,j) .GT. zepsi) THEN
            IF(zwsnow(i,j) .GT. zepsi) THEN
              ztsnow_mult  (i,j,ksn   ) = t_new  (i,j,ksn)
              rho_snow_mult(i,j,ksn,nx) = rho_new(i,j,ksn)
              wtot_snow  (i,j,ksn,nx) = rho_new(i,j,ksn)*zdzh_snow(i,j,ksn)/rho_w
              wliq_snow    (i,j,ksn,nx) = wl_new (i,j,ksn)*zdzh_snow(i,j,ksn)
            ELSE
              ztsnow_mult  (i,j,ksn   ) = t_s(i,j,nx)
              rho_snow_mult(i,j,ksn,nx) = rho_snow_mult(i,j,1,nx)
              wtot_snow  (i,j,ksn,nx) = zwsnew(i,j)/ke_snow
            END IF
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO

! heat conductivity of snow as funtion of water content
  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN  ! for landpoints only
          IF (zwsnew(i,j).GT.zepsi) THEN
            zalas_mult(i,j,ksn) = 2.22_ireals*(rho_snow_mult(i,j,ksn,nx)/rho_i)**1.88_ireals
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        IF (zwsnew(i,j).GT.zepsi) THEN
          zgsb(i,j) = (zalas_mult(i,j,ke_snow)*(-zhm_snow(i,j,ke_snow))+zalam(i,j,1)*zdzms(1))/ &
                      (-zhm_snow(i,j,ke_snow)+zdzms(1)) * &
                      (ztsnow_mult(i,j,ke_snow) - t_so(i,j,1,nx))/(-zhm_snow(i,j,ke_snow)+zdzms(1))
        END IF

        ! total forcing for uppermost soil layer
        zfor_s(i,j) = zrnet_s(i,j) + zshfl_s(i,j) + zlhfl_s(i,j) + zsprs(i,j)*(1._ireals - zf_snow(i,j))       &
                         + zf_snow(i,j) * (1._ireals-ztsnow_pm(i,j)) * zgsb(i,j)

        IF(zwsnew(i,j) .GT. zepsi) THEN
          zrnet_snow = sobs(i,j) * (1.0_ireals - EXP(-zextinct(i,j,1)*zdzm_snow(i,j,1))) + zthsnw(i,j)
        ELSE
          zrnet_snow = sobs(i,j) + zthsnw(i,j)
        END IF
        zshfl_snow = zrhoch(i,j)*cp_d*(zth_low(i,j) - ztsnow_mult(i,j,1))
        zlhfl_snow = lh_s*zversn(i,j) 
        ! zfor_snow_mult should be multiplied by zf_snow!
        ! (in one-layer snow model zf_snow is hidden in zdz_snow_fl)
        zfor_snow_mult(i,j)  = zf_snow(i,j)*(zrnet_snow + zshfl_snow + zlhfl_snow + zsprs(i,j))

      END IF          ! land-points only
    END DO
   END DO

  ELSE

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        zsprs(i,j)      = 0.0_ireals
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i,j)*zrs(i,j) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs(i,j)        = - lh_f*zrs(i,j)
          zdwidt (i,j) = zdwidt (i,j) + zrs(i,j)
          zdwsndt(i,j) = zdwsndt(i,j) - zrs(i,j)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
          zwimax       = cwimax_ml*(1._ireals + plcov(i,j)*5._ireals)
          zwinstr      = zwin(i,j) + zdwidt(i,j)*zdtdrhw
          IF (zwinstr > zwimax) THEN  ! overflow of interception store
            zro        = (zwinstr - zwimax)*zrhwddt
            zdwidt(i,j)= zdwidt(i,j) - zro
            runoff_s(i,j)  = runoff_s(i,j) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF (zwsnow(i,j) == 0.0_ireals .AND.                            &
               (1._ireals-ztsnow_pm(i,j))*zrr(i,j) > 0.0_ireals) THEN
          zsprs(i,j)   = lh_f*zrr(i,j)
          zdwidt (i,j) = zdwidt (i,j) - zrr(i,j)
          zdwsndt(i,j) = zdwsndt(i,j) + zrr(i,j)
        END IF

!       Influence of heatflux through snow on total forcing:
        zwsnew(i,j)   = zwsnow(i,j) + zdwsndt(i,j)*zdtdrhw
        IF (zwsnew(i,j).GT.zepsi) THEN
!         heat conductivity of snow as funtion of water content
! BR      zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnow(i,j)))
!
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
          zalas  = 2.22_ireals*(rho_snow(i,j,nx)/rho_i)**1.88_ireals

! BR 11/2005 Use alternative formulation for heat conductivity by Sun et al., 1999
!            The water vapour transport associated conductivity is not included.

!        zalas   = 0.023_ireals+(2.290_ireals-0.023_ireals)* &
!                               (7.750E-05_ireals*rho_snow(i,j,nx) + &
!                                1.105E-06_ireals*prho_snow(i,j,nx)**2)

          zgsb(i,j) = zalas*(ztsnow(i,j) - zts(i,j))/zdz_snow_fl(i,j)
!<em: should be divided by zf_snow?
!em          zgsb(i,j) = zalas*(ztsnow(i,j) - zts(i,j))/zdz_snow_fl(i,j)/MAX(zf_snow(i,j),zepsi)
        END IF

        ! total forcing for uppermost soil layer
        zfor_s(i,j) = zrnet_s(i,j) + zshfl_s(i,j) + zlhfl_s(i,j) + zsprs(i,j)       &
                         + zf_snow(i,j) * (1._ireals-ztsnow_pm(i,j)) * zgsb(i,j)
      END IF          ! land-points only
    END DO
   END DO

 ENDIF ! lmulti_snow

!------------------------------------------------------------------------------
! Section II.6: Solution of the heat conduction equation, freezing/melting
!               of soil water/ice (optionally)
!------------------------------------------------------------------------------

  DO kso = 2, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          ! for heat conductivity: zalam is now 3D
          zakb = zalam(i,j,kso)/zroc(i,j,kso)
          zaga(i,j,kso) = -zalfa*zdt*zakb/(zdzhs(kso)*zdzms(kso))
          zagc(i,j,kso) = -zalfa*zdt*zakb/(zdzhs(kso)*zdzms(kso+1))
          zagb(i,j,kso) = 1._ireals - zaga(i,j,kso) - zagc(i,j,kso)
          zagd(i,j,kso) = t_so(i,j,kso,nx) +                                     &
                 (1._ireals - zalfa)*( - zaga(i,j,kso)/zalfa*t_so(i,j,kso-1,nx)+ &
                 (zaga(i,j,kso)/zalfa + zagc(i,j,kso)/zalfa)*t_so(i,j,kso,nx) -  &
                  zagc(i,j,kso)/zalfa*t_so(i,j,kso+1,nx)  )
        END IF  ! land-points only
      END DO
    END DO
  END DO        ! soil layers

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        ! for heat conductivity: zalam is now 3D: here we need layer 1
        zakb = zalam(i,j,1)/zroc(i,j,1)
        zaga(i,j,  1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(1))
        zagc(i,j,  1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(2))
        zagb(i,j,  1) = 1._ireals - zaga(i,j,1) - zagc(i,j,1)
        zagd(i,j,  1) = t_so(i,j,1,nx) + (1._ireals - zalfa)* (                 &
                        - zaga(i,j,1)/zalfa * t_s(i,j,nx) +                     &
                        (zaga(i,j,1) + zagc(i,j,1))/zalfa * t_so(i,j,1,nx) -    &
                         zagc(i,j,1)/zalfa * t_so(i,j,2,nx)   )
        zaga(i,j,0)    = 0.0_ireals
        zagb(i,j,0)    = zalfa
        zagc(i,j,0)    = -zalfa
        zagd(i,j,0)    = zdzms(1) * zfor_s(i,j)/zalam(i,j,1)+(1._ireals-zalfa)* &
                        (t_so(i,j,1,nx) - t_s(i,j,nx))
        zaga(i,j,ke_soil+1) = 0.0_ireals
        zagb(i,j,ke_soil+1) = 1.0_ireals
        zagc(i,j,ke_soil+1) = 0.0_ireals
        zagd(i,j,ke_soil+1) = t_so(i,j,ke_soil+1,nx)
      END IF          ! land-points only
    END DO
  END DO


  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        zagc(i,j,0) = zagc(i,j,0)/zagb(i,j,0)
        zagd(i,j,0) = zagd(i,j,0)/zagb(i,j,0)
      END IF          ! land-points only
    END DO
  END DO


  DO kso=1,ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          zzz = 1._ireals/(zagb(i,j,kso) - zaga(i,j,kso)*zagc(i,j,kso-1))
          zagc(i,j,kso) = zagc(i,j,kso) * zzz
          zagd(i,j,kso) = (zagd(i,j,kso) - zaga(i,j,kso)*zagd(i,j,kso-1)) * zzz
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! soil layers

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
      zage(i,j,ke_soil+1) = (zagd(i,j,ke_soil+1) - zaga(i,j,ke_soil+1)*       &
                             zagd(i,j,ke_soil))/                              &
                            (zagb(i,j,ke_soil+1) - zaga(i,j,ke_soil+1)*       &
                             zagc(i,j,ke_soil))
      END IF          ! land-points only
    END DO
  END DO

  DO kso = ke_soil,0,-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          zage(i,j,kso)     = zagd(i,j,kso) - zagc(i,j,kso)*zage(i,j,kso+1)
        ! The surface temperature computed by t_so(i,j,0,nnew)=zage(i,j,0) is
        ! presently unused
          t_so(i,j,kso,nnew) = zage(i,j,kso)
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! soil layers

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        t_so(i,j,ke_soil+1,nnew) = zage(i,j,ke_soil+1) ! climate value, unchanged
      END IF          ! land-points only
    END DO
  END DO

IF (lmulti_snow) THEN

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF(zwsnew(i,j) .GE. zepsi) THEN
          zrocs(i,j) = wliq_snow(i,j,1,nx)/zdzh_snow(i,j,1)*chc_w*rho_w + &
            (wtot_snow(i,j,1,nx)-wliq_snow(i,j,1,nx))/zdzh_snow(i,j,1)*chc_i*rho_i
          zakb      = zalas_mult(i,j,1)/zrocs(i,j)  !(chc_i*rho_snow(1,nx))
          zaga(i,j,1) = 0.0_ireals
          zagc(i,j,1) = -zdt*zakb/(zdzh_snow(i,j,1)*zdzm_snow(i,j,2))
          zagb(i,j,1) = 1._ireals - zagc(i,j,1)
          zagd(i,j,1) = ztsnow_mult(i,j,1) + zdt/zdzh_snow(i,j,1)*zfor_snow_mult(i,j)/zrocs(i,j)

          zrocs(i,j) = wliq_snow(i,j,ke_snow,nx)/zdzh_snow(i,j,ke_snow)*chc_w*rho_w + &
            (wtot_snow(i,j,ke_snow,nx) - &
            wliq_snow(i,j,ke_snow,nx))/zdzh_snow(i,j,ke_snow)*chc_i*rho_i
          zakb = zalas_mult(i,j,ke_snow)/zrocs(i,j)  !(chc_i*rho_snow(1,nx))
          zaga(i,j,ke_snow) = -zdt*zakb/(zdzh_snow(i,j,ke_snow)*zdzm_snow(i,j,ke_snow))
          zagb(i,j,ke_snow) = 1.0_ireals - zaga(i,j,ke_snow)
          zagc(i,j,ke_snow) = 0.0_ireals
          zagd(i,j,ke_snow) = ztsnow_mult(i,j,ke_snow) &
                              - zdt/zdzh_snow(i,j,ke_snow)*zgsb(i,j)/zrocs(i,j)
        END IF
      END IF          ! land-points only
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF(zwsnew(i,j) .GE. zepsi) THEN
          zagc(i,j,1) = zagc(i,j,1)/zagb(i,j,1)
          zagd(i,j,1) = zagd(i,j,1)/zagb(i,j,1)
        END IF
      END IF          ! land-points only
    END DO
  END DO

  DO ksn = 2, ke_snow-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zepsi) THEN
            zrocs(i,j) = wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn)*chc_w*rho_w + &
              (wtot_snow(i,j,ksn,nx) - wliq_snow(i,j,ksn,nx))/zdzh_snow(i,j,ksn)*chc_i*rho_i
            zakb = zalas_mult(i,j,ksn)/zrocs(i,j) !(chc_i*rho_snow(ksn,nx))
            zaga(i,j,ksn) = -zdt*zakb/(zdzh_snow(i,j,ksn)*zdzm_snow(i,j,ksn))
            zagc(i,j,ksn) = -zdt*zakb/(zdzh_snow(i,j,ksn)*zdzm_snow(i,j,ksn+1))
            zagb(i,j,ksn) = 1._ireals - zaga(i,j,ksn) - zagc(i,j,ksn)
            zagd(i,j,ksn) = ztsnow_mult(i,j,ksn)
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! snow layers

  DO ksn=2,ke_snow-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zepsi) THEN
            zzz = 1._ireals/(zagb(i,j,ksn) - zaga(i,j,ksn)*zagc(i,j,ksn-1))
            zagc(i,j,ksn) = zagc(i,j,ksn) * zzz
            zagd(i,j,ksn) = (zagd(i,j,ksn) - zaga(i,j,ksn)*zagd(i,j,ksn-1)) * zzz
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! snow layers

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF(zwsnew(i,j) .GE. zepsi) THEN
          zage(i,j,ke_snow) = (zagd(i,j,ke_snow) - zaga(i,j,ke_snow)* zagd(i,j,ke_snow-1))/ &
                              (zagb(i,j,ke_snow) - zaga(i,j,ke_snow)* zagc(i,j,ke_snow-1))
          ztsnown_mult(i,j,ke_snow) = zage(i,j,ke_snow)
        END IF
      END IF          ! land-points only
    END DO
  END DO

  DO ksn = ke_snow-1,1,-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zepsi) THEN
            zage(i,j,ksn) = zagd(i,j,ksn) - zagc(i,j,ksn)*zage(i,j,ksn+1)
            ztsnown_mult(i,j,ksn) = zage(i,j,ksn)
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO                ! snow layers

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        ztsnow(i,j) = 0._ireals
        ztsn  (i,j) = t_so(i,j,1,nnew)
        IF(zwsnew(i,j) .GE. zepsi) THEN
          zrocs(i,j) = wliq_snow(i,j,1,nx)/zdzh_snow(i,j,1)*chc_w*rho_w + &
            (wtot_snow(i,j,1,nx)-wliq_snow(i,j,1,nx))/zdzh_snow(i,j,1)*chc_i*rho_i
          zswitch(i,j) = MAX(-zfor_snow_mult(i,j)/50./zrocs(i,j)*zdt*ke_snow, &
                              zgsb(i,j)/50./zrocs(i,j)*zdt*ke_snow)
          zswitch(i,j) = MAX(zswitch(i,j),1.E-03_ireals)
        END IF
      END IF          ! land-points only
    END DO
  END DO
  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          ztsnow(i,j) = ztsnow(i,j) + ztsnow_mult(i,j,ksn)*zdzh_snow(i,j,ksn)/h_snow(i,j,nx)
        END IF          ! land-points only
      END DO
    END DO
  END DO
  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zepsi .AND. zwsnew(i,j) .LT. zswitch(i,j)) THEN
            ! in case of very thin snow and large negative heat balance at the snow surface
            ! use one-layer snow model
            tmp_num = ztsnow(i,j) + zdt*2._ireals*(zfor_snow_mult(i,j) - zgsb(i,j)*zf_snow(i,j))  &
                      /zrocs(i,j)/(zswitch(i,j)/rho_snow_mult(i,j,1,nx)*rho_w) - ( ztsn(i,j) - zts(i,j) )
            zalas  = 2.22_ireals*(rho_snow_mult(i,j,1,nx)/rho_i)**1.88_ireals
         
            ztsnow_im = - zrhoch(i,j) * (cp_d + zdqvtsnow(i,j) * lh_s) - zalas/h_snow(i,j,nx)
            zfak  = MAX(zepsi,1.0_ireals - zdt*zalfa*ztsnow_im/zrocs(i,j)/h_snow(i,j,nx))
            tmp_num = ztsnow(i,j) + (tmp_num-ztsnow(i,j))/zfak
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO
  DO ksn = 0, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zepsi .AND. zwsnew(i,j) .LT. zswitch(i,j)) THEN
            ztsnown_mult(i,j,ksn) = tmp_num
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF(zwsnew(i,j) .GT. zepsi .and. zwsnew(i,j) .LT. zswitch(i,j)) THEN
          ! in case of very thin snow and large positive heat balance at the snow surface
          ! melt available snow as a one layer
          IF(zfor_snow_mult(i,j)*zdt > zwsnew(i,j)*rho_w*lh_f) THEN
            zfor_snow_mult(i,j) = 0._ireals
            zdwsndt(i,j) = zdwsndt(i,j) - zwsnew(i,j)*rho_w/zdt
            zwsnew(i,j)  = 0._ireals
          END IF
        END IF
      END IF  ! land-points only
    END DO
  END DO
  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF(zwsnew(i,j) .GE. zepsi .AND. zwsnew(i,j) .LT. zswitch(i,j)) THEN
            IF(zfor_snow_mult(i,j) .GT. 0._ireals) THEN
              ztsnown_mult(i,j,ksn) = ztsnow_mult(i,j,ksn) + &
                zfor_snow_mult(i,j)*zdt/(chc_i*wtot_snow(i,j,ksn,nx))/rho_w/ke_snow
              ztsnown_mult(i,j,0) = ztsnown_mult(i,j,1)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO

END IF    ! lmulti_snow

IF(lmelt) THEN
  IF(.NOT.lmelt_var) THEN
    DO kso = 1,ke_soil
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN ! land points only,
            IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
              ! melting or freezing of soil water
              zenergy   = zroc(i,j,kso)*zdzhs(kso)*(t_so(i,j,kso,nnew) - t0_melt)
              zdwi_max  = - zenergy/(lh_f*rho_w)
              zdelwice  = zdwi_max
              zwso_new  = w_so(i,j,kso,nx) + zdt*zdwgdt(i,j,kso)/rho_w
              IF (zdelwice.LT.0.0_ireals)                                    &
                        zdelwice = - MIN( - zdelwice,w_so_ice(i,j,kso,nx))
              IF (zdelwice.GT.0.0_ireals)                                    &
                        zdelwice =   MIN(   zdelwice,zwso_new - w_so_ice(i,j,kso,nx))
              w_so_ice(i,j,kso,nnew) = w_so_ice(i,j,kso,nx) + zdelwice
!             At this point we have 0.0 LE w_so_ice(i,j,kso,nnew) LE w_so(i,j,kso,nnew)
!             If we have 0.0 LT w_so_ice(i,j,kso,nnew) LT w_so(i,j,kso,nnew)
!             the resulting new temperature has to be equal to t0_melt.
!             If not all energy available can be used to melt/freeze soil water, the
!             following line corrects the temperature. It also applies to cases
!             where all soil water is completely ice or water, respectively.
              t_so(i,j,kso,nnew) = t0_melt + (zdelwice - zdwi_max)*(lh_f*rho_w)/     &
                                            (zroc(i,j,kso)*zdzhs(kso))
            END IF  ! m_styp > 2
          END IF  ! land-points only
        END DO
      END DO
    END DO        ! soil layers
  ELSE IF(lmelt_var) THEN
    DO kso = 1,ke_soil
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN ! land points only,
            IF (m_styp(i,j).ge.3) THEN ! neither ice or rocks
              ztx      = t0_melt
              zw_m(i,j)     = zporv(i,j)*zdzhs(kso)
              IF(t_so(i,j,kso,nnew).LT.(t0_melt-zepsi)) THEN
                zaa    = g*zpsis(i,j)/lh_f
                zw_m(i,j) = zw_m(i,j)*((t_so(i,j,kso,nnew) - t0_melt)/(t_so(i,j,kso,nnew)*  &
                         zaa))**(-zedb(i,j))
                zliquid= MAX(zepsi,w_so(i,j,kso,nx) -  w_so_ice(i,j,kso,nx))
                znen   = 1._ireals-zaa*(zporv(i,j)*zdzhs(kso)/zliquid)**zb_por(i,j)
                ztx    = t0_melt/znen
              ENDIF
              ztx      = MIN(t0_melt,ztx)
              zenergy  = zroc(i,j,kso)*zdzhs(kso)*(t_so(i,j,kso,nnew)-ztx)
              zdwi_max = - zenergy/(lh_f*rho_w)
              zdelwice = zdwi_max
              zwso_new  = w_so(i,j,kso,nx) + zdt*zdwgdt(i,j,kso)/rho_w
              zargu = zwso_new - zw_m(i,j) - w_so_ice(i,j,kso,nx)
              IF (zdelwice.LT.0.0_ireals) zdelwice =                           &
                      - MIN( - zdelwice,MIN(-zargu,w_so_ice(i,j,kso,nx)))
              IF (zdelwice.GT.0.0_ireals) zdelwice =                           &
                        MIN(   zdelwice,MAX( zargu,0.0_ireals))
              w_so_ice(i,j,kso,nnew) = w_so_ice(i,j,kso,nx) + zdelwice
!             At this point we have 0.0 LE w_so_ice(i,j,kso,nnew) LE zwso_new
!             If we have 0.0 LT w_so_ice(i,j,kso,nnew) LT zwso_new
!             the resulting new temperature has to be equal to ztx.
!             If not all energy available can be used to melt/freeze soil water,
!             the following line corrects the temperature. It also applies
!             to cases without any freezing/melting, in these cases the
!             original temperature is reproduced.
              t_so(i,j,kso,nnew) = ztx + (zdelwice - zdwi_max)*       &
                                   (lh_f*rho_w)/(zroc(i,j,kso)*zdzhs(kso))
            END IF                  ! m_stpy > 2
          END IF                    ! land-points only
        END DO
      END DO
    ENDDO
  ENDIF   ! lmelt_var
END IF ! lmelt

!------------------------------------------------------------------------------
! Section II.7: Energy budget and temperature prediction at snow-surface
!------------------------------------------------------------------------------
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        ! next line has to be changed if a soil surface temperature is
        ! predicted by the heat conduction equation
        zdtsdt (i,j) = (t_so(i,j,1,nnew) - zts(i,j))*z1d2dt
        ztsn   (i,j) =  t_so(i,j,1,nnew)
        IF(.NOT. lmulti_snow)  &
          ztsnown(i,j) = ztsn(i,j)       ! default setting
        zwsnew(i,j)       = zwsnow(i,j) + zdwsndt(i,j)*zdtdrhw

        ! forcing contributions for snow formation of dew and rime are
        ! contained in ze_ges, heat fluxes must not be multiplied by
        ! snow covered fraction

        IF(.NOT. lmulti_snow) THEN

          zrnet_snow = sobs(i,j) + zthsnw(i,j)
          zshfl_snow = zrhoch(i,j)*cp_d*(zth_low(i,j) - ztsnow(i,j))
          zlhfl_snow = lh_s*zversn(i,j)
          zfor_snow  = zrnet_snow + zshfl_snow + zlhfl_snow

          ! forecast of snow temperature Tsnow
          IF (ztsnow(i,j) < t0_melt .AND. zwsnew(i,j) > zepsi) THEN
            ztsnown(i,j) = ztsnow(i,j) + zdt*2._ireals*(zfor_snow - zgsb(i,j))  &
                           /zrocs(i,j) - ( ztsn(i,j) - zts(i,j) )
            ! implicit formulation
! BR        zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnew(i,j)))
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
            zalas  = 2.22_ireals*(rho_snow(i,j,nx)/rho_i)**1.88_ireals

            ztsnow_im    = - zrhoch(i,j) * (cp_d + zdqvtsnow(i,j) * lh_s)       &
                                         - zalas/zdz_snow_fl(i,j)
            zfak  = MAX(zepsi,1.0_ireals - zdt*zalfa*ztsnow_im/zrocs(i,j))
            ztsnown(i,j) = ztsnow(i,j) + (ztsnown(i,j)-ztsnow(i,j))/zfak
          END IF

          zdtsnowdt(i,j) = (ztsnown(i,j) - ztsnow(i,j))*z1d2dt
        ENDIF

      END IF          ! land-points only
    END DO
  END DO

  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN          ! land-points only
            IF (zwsnew(i,j) > zepsi) THEN
              zdtsnowdt_mult(i,j,ksn) = (ztsnown_mult(i,j,ksn) - ztsnow_mult(i,j,ksn))*z1d2dt
            END IF
          ENDIF
        END DO
      END DO
    ENDDO
  END IF

!------------------------------------------------------------------------------
! Section II.8: Melting of snow ,infiltration and surface runoff of snow water
!------------------------------------------------------------------------------

! If the soil surface temperature predicted by the equation of heat conduction
! is used instead of using T_s = T_so(1), the following section has to be
! adjusted accordingly.
!
! Basically this snow model uses heat fluxes to either heat the uppermost soil
! layer, if the snow surface temperature exceeds t0_melt, or to heat the snow, if
! the temperature of the uppermost soil layer exceeds t0_melt. Melting is considered
! after this process, if the snow temperature equals t0_melt AND the temperature of
! the uppermost soil layer exceeds t0_melt. The excess heat (t_so(1)-t0_melt) is used for
! melting snow. In cases this melting may be postponed to the next time step.

IF (.NOT. lmulti_snow) THEN

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN          ! land-points only

      zwsnn  (i,j)  = zwsnow(i,j) + zdtdrhw*zdwsndt(i,j)
      zwsnew (i,j)  = zwsnn(i,j)
      ztsnownew     = ztsnown(i,j)
      ze_avail      = 0.0_ireals
      ze_total      = 0.0_ireals
      zfr_melt      = 0.0_ireals

      IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only
        ! first case: T_snow > t0_melt: melting from above
        ! ----------
        IF (ztsnown(i,j) > t0_melt .AND. t_so(i,j,1,nnew) < t0_melt ) THEN
          zdwsnm(i,j)   = zwsnew(i,j)*.5_ireals*(ztsnown(i,j) - (t0_melt - zepsi))/ &
                          (.5_ireals* (zts(i,j) - (t0_melt - zepsi)) - lh_f/chc_i)
          zdwsnm(i,j)   = zdwsnm(i,j)*z1d2dt*rho_w 
!
!HJP 2010-07-02 Begin
          snow_melt(i,j)= snow_melt(i,j) - zdwsnm(i,j)/z1d2dt
!HJP 2010-07-02 End
!
          zdwsndt(i,j)  = zdwsndt (i,j) + zdwsnm(i,j)
          ztsnownew       = t0_melt - zepsi
          zdtsnowdt(i,j)  = zdtsnowdt(i,j) + (ztsnownew - ztsnown(i,j))*z1d2dt
          runoff_s (i,j)= runoff_s(i,j) - zdwsnm(i,j)*zroffdt
        ENDIF ! melting from above

        IF (t_so(i,j,1,nnew) .gt. t0_melt) THEN
          !second case:  temperature of uppermost soil layer > t0_melt. First a
          !-----------   heat redistribution is performed. As a second step,
          !              melting of snow is considered.
          ! a) Heat redistribution
          ztsnew = t0_melt + zepsi
          ztsnownew      = ztsn(i,j) + ztsnown(i,j) - ztsnew +  &
               2._ireals*(ztsn(i,j) - ztsnew)*zroc(i,j,1)*zdzhs(1)/zrocs(i,j)
          zdtsdt(i,j)    = zdtsdt(i,j) + (ztsnew - t_so(i,j,1,nnew))*z1d2dt
          zdtsnowdt(i,j) = zdtsnowdt(i,j) + (ztsnownew - ztsnown(i,j))*z1d2dt
          ! b) Melting of snow (if possible)
          IF (ztsnownew > t0_melt) THEN
            ze_avail     = 0.5_ireals*(ztsnownew - t0_melt)*zrocs(i,j)
            ze_total     = lh_f*zwsnew(i,j)*rho_w
            zfr_melt     = MIN(1.0_ireals,ze_avail/ze_total)
            zdtsnowdt(i,j)= zdtsnowdt(i,j) + (t0_melt - ztsnownew)*z1d2dt
            zdelt_s      = MAX(0.0_ireals,(ze_avail - ze_total)/(zroc(i,j,1)* &
                                                                    zdzhs(1)))
            zdtsdt(i,j)  = zdtsdt(i,j) + zdelt_s*z1d2dt
            ! melted snow is allowed to penetrate the soil (up to field
            ! capacity), if the soil type is neither ice nor rock (zrock = 0);
            ! else it contributes to surface run-off;
            ! fractional water content of the first soil layer determines
            ! a reduction factor which controls additional run-off
            zdwsnm(i,j)   = zfr_melt*zwsnew(i,j)*z1d2dt*rho_w  ! available water
!
!HJP 2010-07-02 Begin
            snow_melt(i,j)= snow_melt(i,j) + zdwsnm(i,j)/z1d2dt
!HJP 2010-07-02 End
!
            zdwsndt(i,j)  = zdwsndt (i,j) - zdwsnm(i,j)
            zdwgme        = zdwsnm(i,j)*zrock(i,j)             ! contribution to w_so
            zro           = (1._ireals - zrock(i,j))*zdwsnm(i,j)      ! surface runoff
            zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,j,1) -  &
                            zfcap(i,j))/MAX(zporv(i,j)-zfcap(i,j), zepsi)))
            zdwgdt(i,j,1) = zdwgdt(i,j,1) + zdwgme*(1._ireals - zredfu)
            zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                   ! for this fraction

            ! zro-, zdw_so_dt-correction in case of pore volume overshooting
            zw_ovpv = MAX(0._ireals, zw_fr(i,j,1)* zdzhs(1) * zrhwddt +  &
                       zdwgdt(i,j,1) - zporv(i,j) * zdzhs(1) * zrhwddt)
            zro = zro + zw_ovpv
            zdwgdt(i,j,1)= zdwgdt(i,j,1) - zw_ovpv


            IF (zfr_melt > 0.9999_ireals) zdwsndt(i,j)= -zwsnow(i,j)*zrhwddt
            runoff_s (i,j)= runoff_s(i,j) + zro*zroffdt
!           IF(runoff_s(i,j).LT.0.) THEN
!             WRITE (yerror,'(A,I5,A,I5,A,F10.5)')                        &
!              'Negative runoff_s for i= ',i,' and j = ',j,': ',runoff_s(i,j)
!             ierror = 1
!           ENDIF
!           IF(runoff_g(i,j).LT.0.) THEN
!             WRITE (yerror,'(A,I5,A,I5,A,F10.5)')                        &
!              'Negative runoff_g for i= ',i,' and j = ',j,': ',runoff_g(i,j)
!             ierror = 1
!           ENDIF
          END IF   ! snow melting
        END IF     ! snow and/or soil temperatures
      END IF       ! points with snow cover only
    END IF         ! land-points only
  END DO
END DO

ELSE       ! new snow scheme

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN          ! land-points only

      zwsnew(i,j) = zwsnow(i,j) + zdtdrhw*zdwsndt(i,j)
      zdwsnm(i,j) = 0.0_ireals

      ze_out  (i,j) = 0.0_ireals
      zqbase  (i,j) = 0.0_ireals
      zcounter(i,j) = 0.0_ireals

      ze_rad(i,j) = 0.0_ireals
      IF(zextinct(i,j,1).gt.0.0_ireals) ze_rad(i,j) = zf_snow(i,j) * sobs(i,j) 

      ztsnownew_mult(i,j,0) = ztsnown_mult(i,j,0)
    END IF         ! land-points only
  END DO
END DO

DO ksn = 1,ke_snow
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only

          zrefr(i,j) = 0.0_ireals
          zmelt(i,j) = 0.0_ireals
          ztsnownew_mult(i,j,ksn) = ztsnown_mult(i,j,ksn)

          IF(zdzh_snow(i,j,ksn) - wliq_snow(i,j,ksn,nx).GT.zepsi .OR. &
            wtot_snow(i,j,ksn,nx) - wliq_snow(i,j,ksn,nx).GT.zepsi) THEN
            zrho_dry_old(i,j) = MAX(wtot_snow(i,j,ksn,nx) - wliq_snow(i,j,ksn,nx), zepsi)*   &
                                rho_w/(zdzh_snow(i,j,ksn) - wliq_snow(i,j,ksn,nx))
          ELSE
            zrho_dry_old(i,j) = rho_w
          END IF


          ztsnownew_mult(i,j,ksn) = (ztsnown_mult(i,j,ksn)*wtot_snow(i,j,ksn,nx) + t0_melt*zqbase(i,j)*zdt)/  &
            (zqbase(i,j)*zdt+wtot_snow(i,j,ksn,nx)) 

          IF(zextinct(i,j,ksn).eq.0.0_ireals) THEN
            ze_in = ze_out(i,j)
          ELSE
            IF(ksn.eq.ke_snow) THEN
              ze_in = ze_out(i,j) + (ze_rad(i,j) - zcounter(i,j))
            ELSEIF(ksn.eq.1) then
              ze_in = ze_rad(i,j) * (EXP (-zextinct(i,j,1)*zdzm_snow(i,j,1)) - EXP (-zextinct(i,j,1)*zdzh_snow(i,j,1)))
              zcounter(i,j) = ze_rad(i,j) * (1._ireals - EXP (-zextinct(i,j,1)*zdzh_snow(i,j,1)))
            ELSE
              ze_in = ze_out(i,j) + (ze_rad(i,j) - zcounter(i,j)) -  &
                (ze_rad(i,j) - zcounter(i,j)) * EXP (-zextinct(i,j,ksn)*zdzh_snow(i,j,ksn))
              zcounter(i,j) = ze_rad(i,j) - (ze_rad(i,j) - zcounter(i,j)) * EXP (-zextinct(i,j,ksn)*zdzh_snow(i,j,ksn))
            END IF
          END IF

          ztsnownew_mult(i,j,ksn) = ztsnownew_mult(i,j,ksn) + ze_in*zdt/(chc_i*wtot_snow(i,j,ksn,nx))/rho_w
          wtot_snow(i,j,ksn,nx) = wtot_snow(i,j,ksn,nx) + zqbase(i,j)*zdt
          wliq_snow  (i,j,ksn,nx) = wliq_snow  (i,j,ksn,nx) + zqbase(i,j)*zdt

          zdzh_snow(i,j,ksn) = zdzh_snow(i,j,ksn) + zqbase(i,j)*zdt

          rho_snow_mult(i,j,ksn,nx) = MAX(wtot_snow(i,j,ksn,nx)*rho_w/zdzh_snow(i,j,ksn),0.0_ireals)

          IF(ztsnownew_mult(i,j,ksn) .GT. t0_melt) THEN

            IF(wtot_snow(i,j,ksn,nx) .LE. wliq_snow(i,j,ksn,nx)) THEN
              ze_out(i,j) = chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt*rho_w
              zmelt(i,j) = 0.0_ireals
            ELSEIF(chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)/lh_f .LE.  &
              wtot_snow(i,j,ksn,nx)-wliq_snow(i,j,ksn,nx)) THEN
              zmelt(i,j) = chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt/lh_f
              ze_out(i,j) = 0.0_ireals
              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) + zmelt(i,j)*zdt
            ELSE
              zmelt(i,j) = (wtot_snow(i,j,ksn,nx)-wliq_snow(i,j,ksn,nx))*z1d2dt
              ze_out(i,j) = chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt*rho_w &
                             - zmelt(i,j)*lh_f*rho_w
              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) + zmelt(i,j)*zdt
            END IF
            ztsnownew_mult(i,j,ksn) = t0_melt

          ELSE
!       T<0
            IF(wliq_snow(i,j,ksn,nx) .GT. -chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)/lh_f) THEN
              zrefr(i,j)            = -chc_i*wtot_snow(i,j,ksn,nx)*(ztsnownew_mult(i,j,ksn) - t0_melt)*z1d2dt/lh_f
              ztsnownew_mult(i,j,ksn)   = t0_melt
              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) - zrefr(i,j)*zdt
            ELSE
              zrefr(i,j)            = wliq_snow(i,j,ksn,nx)*z1d2dt
              wliq_snow(i,j,ksn,nx) = 0.0_ireals
              ztsnownew_mult(i,j,ksn)   = ztsnownew_mult(i,j,ksn) + zrefr(i,j)*zdt*lh_f/(chc_i*wtot_snow(i,j,ksn,nx))
            END IF
            ze_out(i,j) = 0.0_ireals

          END IF

          zdtsnowdt_mult(i,j,ksn) = zdtsnowdt_mult(i,j,ksn) + &
                                    (ztsnownew_mult(i,j,ksn) - ztsnown_mult(i,j,ksn))*z1d2dt
          IF(wtot_snow(i,j,ksn,nx) .LE. wliq_snow(i,j,ksn,nx)) THEN
            zqbase(i,j)           = wliq_snow(i,j,ksn,nx)*z1d2dt
            wliq_snow    (i,j,ksn,nx) = 0.0_ireals
            wtot_snow  (i,j,ksn,nx) = 0.0_ireals
            zdzh_snow    (i,j,ksn)    = 0.0_ireals
            rho_snow_mult(i,j,ksn,nx) = 0.0_ireals
          ELSE
            IF(zrefr(i,j) .GT. 0.0_ireals .OR. zmelt(i,j) .GT. 0.0_ireals) THEN
              zadd_dz = MAX(zrefr(i,j),0._ireals)*(-1.0_ireals + 1.0_ireals/rho_i*rho_w)*zdt
              zadd_dz = MAX(zmelt(i,j),0._ireals)*(-1.0_ireals/zrho_dry_old(i,j)*rho_w + 1.0_ireals)*zdt
              zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn) + zadd_dz
              rho_snow_mult(i,j,ksn,nx) = MAX(wtot_snow(i,j,ksn,nx)*rho_w/zdzh_snow(i,j,ksn),0.0_ireals)
              IF(wtot_snow(i,j,ksn,nx) .LE. 0.0_ireals) zdzh_snow(i,j,ksn) = 0.0_ireals
              IF(rho_snow_mult(i,j,ksn,nx) .GT. rho_w) THEN
                zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn)*rho_snow_mult(i,j,ksn,nx)/rho_w
                rho_snow_mult(i,j,ksn,nx) = rho_w
              END IF
            END IF

            zsn_porosity = 1._ireals - (rho_snow_mult(i,j,ksn,nx)/rho_w -  &
                           wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn))/rho_i*rho_w - &
                           wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn)
            zsn_porosity = MAX(zsn_porosity,cwhc + 0.1_ireals)
            zp1 = zsn_porosity - cwhc

            IF (wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn) .GT. cwhc) THEN
              zfukt             = (wliq_snow(i,j,ksn,nx)/zdzh_snow(i,j,ksn) - cwhc)/zp1
              zq0               = chcond * zfukt**3.0_ireals
              zqbase(i,j)       = MIN(zq0*zdt,wliq_snow(i,j,ksn,nx))
              wliq_snow(i,j,ksn,nx) = wliq_snow(i,j,ksn,nx) - zqbase(i,j)
              wtot_snow(i,j,ksn,nx) = wtot_snow(i,j,ksn,nx) - zqbase(i,j)

              zdzh_snow(i,j,ksn) = zdzh_snow(i,j,ksn) - zqbase(i,j)

              IF(zdzh_snow(i,j,ksn) .LT. zepsi/ke_snow) THEN
                zqbase(i,j) = zqbase(i,j) + wtot_snow(i,j,ksn,nx)
!                zdwsndt(i,j)  = zdwsndt(i,j) - wtot_snow(i,j,ksn,nx)*rho_w/zdt
                wliq_snow(i,j,ksn,nx) = 0.0_ireals
                wtot_snow(i,j,ksn,nx) = 0.0_ireals
                zdzh_snow(i,j,ksn)    = 0.0_ireals
                rho_snow_mult(i,j,ksn,nx)  = 0.0_ireals
              ELSE
                rho_snow_mult(i,j,ksn,nx) = MAX(wtot_snow(i,j,ksn,nx)*rho_w/zdzh_snow(i,j,ksn),0.0_ireals)
                IF(wtot_snow(i,j,ksn,nx) .LE. 0.0_ireals) zdzh_snow(i,j,ksn) = 0.0_ireals
                IF(rho_snow_mult(i,j,ksn,nx) .GT. rho_w) THEN
                  zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn)*rho_snow_mult(i,j,ksn,nx)/rho_w
                  rho_snow_mult(i,j,ksn,nx) = rho_w
                END IF
              END IF
              zqbase(i,j) = zqbase(i,j)*z1d2dt
            ELSE
              zqbase(i,j) = 0.0_ireals
            END IF
          END IF

        END IF       ! points with snow cover only
      END IF         ! land-points only
    END DO
  END DO
END DO        ! snow layers

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN          ! land-points only
      IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only
        zdwsnm(i,j) = zqbase(i,j)*rho_w       ! ksn == ke_snow
!        zdwsnm(i,j) = zqbase(i,j)*rho_w*zf_snow(i,j)       ! ksn == ke_snow
      END IF       ! points with snow cover only
    END IF         ! land-points only
  END DO
END DO

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN          ! land-points only
      IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only
        zdwsndt(i,j)  = zdwsndt(i,j) - zdwsnm(i,j)

          ! melted snow is allowed to penetrate the soil (up to field
          ! capacity), if the soil type is neither ice nor rock (zrock = 0);
          ! else it contributes to surface run-off;
          ! fractional water content of the first soil layer determines
          ! a reduction factor which controls additional run-off

        zdwgme        = zdwsnm(i,j)*zrock(i,j)             ! contribution to w_so
        zro           = (1._ireals - zrock(i,j))*zdwsnm(i,j)      ! surface runoff
        zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,j,1) -  &
                        zfcap(i,j))/MAX(zporv(i,j)-zfcap(i,j), zepsi)))
        zdwgdt(i,j,1) = zdwgdt(i,j,1) + zdwgme*(1._ireals - zredfu)
        zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                   ! for this fraction

        ! zro-, zdw_so_dt-correction in case of pore volume overshooting
        zw_ovpv = MAX(0._ireals, zw_fr(i,j,1)* zdzhs(1) * zrhwddt +  &
                  zdwgdt(i,j,1) - zporv(i,j) * zdzhs(1) * zrhwddt)
        zro = zro + zw_ovpv
        zdwgdt(i,j,1)= zdwgdt(i,j,1) - zw_ovpv

        runoff_s(i,j) = runoff_s(i,j) + zro*zroffdt
      END IF       ! points with snow cover only
    END IF         ! land-points only
  END DO
END DO

! snow densification due to gravity and metamorphism 

  DO ksn = 2, ke_snow
    zp(:,:,ksn) = 0.0_ireals                         ! gravity, Pa
    DO k = ksn,1,-1
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN          ! land-points only
            zp(i,j,ksn) = zp(i,j,ksn) + rho_snow_mult(i,j,k,nx)*g*zdzh_snow(i,j,ksn)
          END IF         ! land-points only
        END DO
      END DO
    END DO
  END DO

  DO ksn = 2, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN          ! land-points only
          IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only
            IF(rho_snow_mult(i,j,ksn,nx) .LT. 600._ireals .AND. rho_snow_mult(i,j,ksn,nx) .NE. 0.0_ireals) THEN
              zdens_old = rho_snow_mult(i,j,ksn,nx)
              zeta =         &! compactive viscosity of snow
                ca2*EXP(19.3_ireals*rho_snow_mult(i,j,ksn,nx)/rho_i)* &
                EXP(67300._ireals/8.31_ireals/ztsnownew_mult(i,j,ksn))
              rho_snow_mult(i,j,ksn,nx) = rho_snow_mult(i,j,ksn,nx)+zdt*rho_snow_mult(i,j,ksn,nx)*(csigma+zp(i,j,ksn))/zeta
              rho_snow_mult(i,j,ksn,nx) = MIN(rho_snow_mult(i,j,ksn,nx),rho_i)
              zdzh_snow(i,j,ksn)   = zdzh_snow(i,j,ksn) * zdens_old/rho_snow_mult(i,j,ksn,nx)
            END IF
          END IF       ! points with snow cover only
        END IF         ! land-points only
      END DO
    END DO
  END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF (zwsnew(i,j) > zepsi) THEN        ! points with snow cover only

          IF(ztsnownew_mult(i,j,0) .GT. t0_melt) THEN
            ztsnownew_mult(i,j,0) = t0_melt
            zdtsnowdt_mult(i,j,0) = zdtsnowdt_mult(i,j,0) +     &
                                    (ztsnownew_mult(i,j,0) - ztsnown_mult(i,j,0))*z1d2dt
          END IF
        END IF       ! points with snow cover only
      END IF         ! land-points only
    END DO
  END DO
END IF

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

IF (lmulti_snow) THEN
  ! First for ksn == 0
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        t_snow_mult  (i,j,0,nnew) = t_snow_mult(i,j,0,nx) + zdt*zdtsnowdt_mult(i,j,0)
        t_snow       (i,j  ,nnew) = t_snow_mult(i,j,0,nnew)
        rho_snow     (i,j  ,nnew) = 0.0_ireals
      ENDIF
    ENDDO
  ENDDO

  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN  ! for landpoints only
          t_snow_mult  (i,j,ksn,nnew) = t_snow_mult(i,j,ksn,nx) + zdt*zdtsnowdt_mult(i,j,ksn)

          dzh_snow(i,j,ksn,nnew) = zdzh_snow(i,j,ksn)
          wtot_snow  (i,j,ksn,nnew) = wtot_snow(i,j,ksn,nx)
          rho_snow_mult(i,j,ksn,nnew) = rho_snow_mult(i,j,ksn,nx)
          wliq_snow    (i,j,ksn,nnew) = wliq_snow(i,j,ksn,nx)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        t_snow(i,j,nnew)  = t_snow(i,j,nx) + zdt*zdtsnowdt(i,j)
      ENDIF
    ENDDO
  ENDDO
ENDIF

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN  ! for landpoints only
      ! t_snow is computed above
      ! t_snow(i,j,nnew)  = t_snow(i,j,nx) + zdt*zdtsnowdt(i,j)
      t_so(i,j,1,nnew)  = t_so(i,j,1,nx) + zdt*zdtsdt   (i,j)
      ! Next line has to be changed, if the soil surface temperature
      ! t_so(i,j,0,nnew) predicted by the heat conduction equation is used
      t_s   (i,j,nnew)    = t_so(i,j,1,nnew)
      t_so  (i,j,0,nnew)  = t_so(i,j,1,nnew)
      w_snow(i,j,nnew)  = w_snow(i,j,nx) + zdt*zdwsndt  (i,j)/rho_w
      w_i   (i,j,nnew)  = w_i   (i,j,nx) + zdt*zdwidt   (i,j)/rho_w

      IF (w_snow(i,j,nnew) <= zepsi) THEN
        w_snow(i,j,nnew) = 0.0_ireals
        t_snow(i,j,nnew) = t_so(i,j,0,nnew)
      ENDIF
      IF (w_i(i,j,nnew) <= zepsi) w_i(i,j,nnew) = 0.0_ireals
    END IF          ! land-points only
  END DO
END DO

! Eliminate snow for multi-layer snow model, if w_snow = 0
IF (lmulti_snow) THEN
  DO ksn = 1, ke_snow
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN  ! for landpoints only
          IF (w_snow(i,j,nnew) == 0.0_ireals) THEN
            t_snow_mult(i,j,ksn,nnew) = t_so(i,j,0,nnew)
            wliq_snow(i,j,ksn,nnew) = 0.0_ireals
            wtot_snow(i,j,ksn,nnew) = 0.0_ireals
            rho_snow_mult (i,j,ksn,nnew) = 0.0_ireals
            dzh_snow (i,j,ksn,nnew) = 0.0_ireals
          ENDIF
        END IF          ! land-points only
      END DO
    END DO
  END DO
ENDIF

IF(.NOT. lmulti_snow) THEN

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN  ! for landpoints only
!BR 7/2005 Update snow density
!
!     a) aging of existing snow
!
!    temperature dependence of relaxation/ageing constant
     ztau_snow = crhosmint+(crhosmaxt-crhosmint)*(t_snow(i,j,nnew)-csnow_tmin) &
                                                /(t0_melt         -csnow_tmin)
     ztau_snow = MAX(crhosmint,MIN(crhosmaxt,ztau_snow))
     zrho_snowe= crhosmax_ml+(rho_snow(i,j,nx)-crhosmax_ml)* &
                           EXP(-ztau_snow*zdt/86400._ireals)
!     b) density of fresh snow
!
     zrho_snowf= crhosminf+(crhosmaxf-crhosminf)* (zth_low(i,j)-csnow_tmin) &
                                                /(t0_melt      -csnow_tmin)
     zrho_snowf= MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))
!
!     c) new snow density is weighted average of existing and new snow
!
     IF ( itype_gscp >= 2000 ) THEN
       ! only possible when running the 2-moment microphysics
!!$ UB: does that really make sense to integrate hail into snow density at the surface?
       znorm=MAX(w_snow(i,j,nx)+(prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j)+prh_gsp(i,j))      &
                 *zdtdrhw,zepsi)
       rho_snow(i,j,nnew)  = ( zrho_snowe*w_snow(i,j,nx) + &
                          zrho_snowf*(prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j)) &
                             *zdtdrhw )    /znorm
     ELSEIF ( itype_gscp >= 4 ) THEN
       znorm=MAX(w_snow(i,j,nx)+(prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j))      &
                 *zdtdrhw,zepsi)
       rho_snow(i,j,nnew)  = ( zrho_snowe*w_snow(i,j,nx) + &
                          zrho_snowf*(prs_gsp(i,j)+prs_con(i,j)+prg_gsp(i,j)) &
                             *zdtdrhw )    /znorm
     ELSE
       znorm=MAX(w_snow(i,j,nx)+(prs_gsp(i,j)+prs_con(i,j) )      &
                 *zdtdrhw,zepsi)
       rho_snow(i,j,nnew)  = ( zrho_snowe*w_snow(i,j,nx) + &
                          zrho_snowf*(prs_gsp(i,j)+prs_con(i,j) ) &
                             *zdtdrhw )    /znorm
     ENDIF
     rho_snow(i,j,nnew) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow(i,j,nnew)))
    END IF          ! land-points only
  END DO
END DO

ELSE   ! new snow scheme

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN  ! for landpoints only
      h_snow(i,j,nnew) = 0.0_ireals
    END IF          ! land-points only
  END DO
END DO
DO ksn = 1,ke_snow
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        IF(w_snow(i,j,nnew) .GT. zepsi) THEN
          h_snow(i,j,nnew) = h_snow(i,j,nnew) + zdzh_snow(i,j,ksn)
        END IF
      END IF          ! land-points only
    END DO
  END DO
END DO

DO   j = jstarts, jends
  DO i = istarts, iends
    sum_weight(i,j) = 0.0_ireals
  END DO
END DO

DO ksn = ke_snow,1,-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        IF(w_snow(i,j,nnew) .GT. zepsi) THEN
          dz_old(i,j,ksn) = dzh_snow(i,j,ksn,nnew)
          z_old(i,j,ksn) = -sum_weight(i,j) - dzh_snow(i,j,ksn,nnew)/2._ireals
          sum_weight(i,j) = sum_weight(i,j) + dzh_snow(i,j,ksn,nnew)
          zhh_snow(i,j,ksn) = -h_snow(i,j,nnew)/ke_snow*(ke_snow-ksn)
        END IF
      END IF          ! land-points only
    END DO
  END DO
END DO

DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN  ! for landpoints only
      IF(w_snow(i,j,nnew) .GT. zepsi) THEN
        zhm_snow(i,j,1) = (-h_snow(i,j,nnew) + zhh_snow(i,j,1))/2._ireals
        !layer thickness betw. half levels of uppermost snow layer
        dzh_snow (i,j,1,nnew) = zhh_snow(i,j,1) + h_snow(i,j,nnew)            
        !layer thickness between snow surface and main level of uppermost layer
        zdzm_snow(i,j,1     ) = zhm_snow(i,j,1) + h_snow(i,j,nnew)            
        IF(dz_old(i,j,1).ne.0..and.rho_snow_mult(i,j,1,nnew).ne.0.) THEN
          wliq_snow(i,j,1,nnew) = wliq_snow(i,j,1,nnew)/dz_old(i,j,1)
        END IF
      END IF
    END IF          ! land-points only
  END DO
END DO
DO ksn = 2,ke_snow
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        IF(w_snow(i,j,nnew) .GT. zepsi) THEN
          zhm_snow(i,j,ksn) = (zhh_snow(i,j,ksn) + zhh_snow(i,j,ksn-1))/2._ireals
          dzh_snow (i,j,ksn,nnew) = zhh_snow(i,j,ksn) - zhh_snow(i,j,ksn-1) ! layer thickness betw. half levels
          zdzm_snow(i,j,ksn     ) = zhm_snow(i,j,ksn) - zhm_snow(i,j,ksn-1) ! layer thickness betw. main levels
          IF(dz_old(i,j,ksn).ne.0..and.rho_snow_mult(i,j,ksn,nnew).ne.0.) THEN
            wliq_snow(i,j,ksn,nnew) = wliq_snow(i,j,ksn,nnew)/dz_old(i,j,ksn)
          END IF
        END IF
      END IF          ! land-points only
    END DO
  END DO
END DO

DO ksn = ke_snow,1,-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      t_new  (i,j,ksn) = 0.0_ireals
      rho_new(i,j,ksn) = 0.0_ireals
      wl_new (i,j,ksn) = 0.0_ireals
    END DO
  END DO

  DO k = ke_snow,1,-1
    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN  ! for landpoints only
          IF(w_snow(i,j,nnew) .GT. zepsi) THEN

            weight = MAX(MIN(z_old(i,j,k)+dz_old(i,j,k)/2._ireals,zhm_snow(i,j,ksn)+dzh_snow(i,j,ksn,nnew)/2._ireals)-   &
                     MAX(z_old(i,j,k)-dz_old(i,j,k)/2._ireals, &
                     zhm_snow(i,j,ksn)-dzh_snow(i,j,ksn,nnew)/2._ireals),0._ireals)/dzh_snow(i,j,ksn,nnew)

            t_new  (i,j,ksn) = t_new  (i,j,ksn) + t_snow_mult  (i,j,k,nnew)*weight
            rho_new(i,j,ksn) = rho_new(i,j,ksn) + rho_snow_mult(i,j,k,nnew)*weight
            wl_new (i,j,ksn) = wl_new (i,j,ksn) + wliq_snow    (i,j,k,nnew)*weight
          END IF
        END IF          ! land-points only
      END DO
    END DO
  END DO
END DO
DO ksn = ke_snow,1,-1
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        IF(w_snow(i,j,nnew) .GT. zepsi) THEN
          t_snow_mult  (i,j,ksn,nnew) = t_new  (i,j,ksn)
          rho_snow_mult(i,j,ksn,nnew) = rho_new(i,j,ksn)
          wtot_snow  (i,j,ksn,nnew) = rho_new(i,j,ksn)*dzh_snow(i,j,ksn,nnew)/rho_w
          wliq_snow    (i,j,ksn,nnew) = wl_new (i,j,ksn)*dzh_snow(i,j,ksn,nnew)
        END IF
      END IF          ! land-points only
    END DO
  END DO
END DO
DO   j = jstarts, jends
  DO i = istarts, iends
    IF (llandmask(i,j)) THEN  ! for landpoints only
      IF(w_snow(i,j,nnew) .GT. zepsi) rho_snow(i,j,nnew) = w_snow(i,j,nnew)/h_snow(i,j,nnew)*rho_w
      t_snow(i,j,nnew) = t_snow_mult(i,j,1,nnew)
      IF(w_snow(i,j,nnew) .GT. zepsi) &
        zrho_snow_old(i,j) = (wtot_snow(i,j,1,nnew)+wtot_snow(i,j,2,nnew)) /  &
                      (dzh_snow(i,j,1,nnew)+dzh_snow(i,j,2,nnew))*rho_w
    END IF          ! land-points only
  END DO
END DO

ENDIF ! lmulti_snow

DO kso = 1,ke_soil
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        w_so(i,j,kso,nnew) = w_so(i,j,kso,nx) + zdt*zdwgdt(i,j,kso)/rho_w
      END IF  ! land-points only
    END DO
  END DO
END DO        ! soil layers

! LS, b, calculation of saturation
DO kso = 1,ke_soil+1
  DO   j = jstartpar, jendpar
    DO i = istartpar, iendpar
      IF (llandmask(i,j)) THEN  ! for landpoints only

        s_so(i,j,kso)=w_so(i,j,kso,nnew)*100.0_ireals/(zporv(i,j)*zdzhs(kso))
        if (s_so(i,j,kso).lt.0.0_ireals) then
           print*,"negative",i,j,kso,my_cart_id,s_so(i,j,kso)
        end if
      END IF  ! land-points only
    END DO
  END DO
END DO        ! soil layers
! LS, e

! computation of the temperature at the boundary soil/snow-atmosphere
CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), t_s(:,:,nnew),        &
             w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow,     &
             istarts, iends, jstarts, jends )

#ifdef NECSX
  CALL collapse(.FALSE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif

!------------------------------------------------------------------------------
! End of module procedure terra_multlay
!------------------------------------------------------------------------------

END SUBROUTINE terra_multlay

!==============================================================================

!------------------------------------------------------------------------------
! End of module src_soil_multlay
!------------------------------------------------------------------------------

END MODULE src_soil_multlay
