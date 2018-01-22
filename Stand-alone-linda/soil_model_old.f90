!------------------------------------------------------------------------------
! Stand-alone implementation of the unsaturated flow in a soil column
! which is a part of the module "src_soil_multlay" in the COSMO model.
!
! Adapted and simplified by Lukas Strebel as his bachelor thesis.
!
! Note that still a lot of parameters that would not be needed (e.g. number of grid points in zonal direction) in a stand-alone 1 column version are left here for comparability and compatibility reasons.
!
! This is the main program and requires the modules 
! - data_namelist    in which the variables and paramters usually provided by other processes are provided.
! - data_soil    in which the soil properties are contained.
! - data_parameters    in which integer and double are redefined
!
! Cautionary note: for small time steps and large running time this program can generate a lot of data fast (which it saves in plain text format).
!
!------------------------------------------------------------------------------
PROGRAM soil_model

USE data_parameters, ONLY :   &
    ireals,       & ! KIND-type parameter for real variables
    iintegers       ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_namelist, ONLY :   &

! horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke_soil,      & ! number of layers in multi-layer soil model
    czmls,        & ! depth of the main level soil layers in m
    
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           &  ! long time-step
    starttime,    &  ! 
    endtime,      &  !
    outtime,      &  !
    nsteps,       &  !
    noutput,      &  ! 
! constants for parametrizations
! ---------------------------------
    rho_w,        &  ! density of liquid water  
!------------------------------------------------------------------------------
! external parameter fields                                        (unit)
! ----------------------------
    soiltyp    ,    & ! type of the soil (keys 0-9)                     --
    rootdp     ,    & ! depth of the roots                            ( m  )
    plcov      ,    & ! fraction of plant cover                         --
    llandmask  ,    & ! landpoint mask  
    czbot_w_so ,    & ! depth of last hydrological active soil layer in m                                --
    czrootdp   ,    & ! LINDA, rootdepth in m                           --
    itype_hydbound, & ! LINDA
! fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    w_so_init      ,     & ! initial total water conent (liquid water)       (m H20)
    prec           ,     & ! precipitation rate 
    evap           ,     & ! evaporation rate 

    init_fields ! init subfunction

USE data_soil       ! All variables from data module "data_soil are used by
                    ! this module. These variables start with letter "c"

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------
USE netcdf

IMPLICIT NONE


!     Definition of temperature/water related variables in the soil model
!!
!     ___________________________|________________________________
!     ////////////////////////////////////////////////////////////
!
!     Layer k=1 - - - - -  W_SO(k) - - - - - - - - - - - - - - - - 
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

! Local parameters:
! ----------------

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals , & ! security constant
    zalfa  = 1.0_ireals     ! degree of impliciteness (1: full implicit,
                               !    (0.5: Cranck-Nicholson)

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &

!   Indices

    nx             , & ! time-level for integration
    nnew           , & ! time-level for next time step
    nout           , & ! output index for this time step
    kso            , & ! loop index for soil moisture layers           
    ke_soil_hy     , & ! number of active soil moisture layers
    i              , & ! loop index in x-direction              
    j              , & ! loop index in y-direction              
    mstyp          , & ! soil type index
    msr_off        , & ! number of layers contributing to surface run off
    ibot_w_so      , & ! index of last hydrological active layer
    irootdp        , & ! LINDA, index of last layer in root zone
    istarts        , & ! start index for x-direction      
    iends          , & ! end   index for x-direction     
    jstarts        , & ! start index for y-direction     
    jends              ! end   index for y-directin      

  REAL    (KIND=ireals   ) ::  &

!   Timestep parameters

    zdt            , & ! integration time-step [s]
    time           , & ! current time [s]
    zdtdrhw        , & ! timestep / density of water
    zrhwddt        , & ! density of water / timestep
    zroffdt         ! timestep for runoff calculation
    

  REAL    (KIND=ireals   ) ::  &

!   Hydraulic parameters

    zinfmx         , & ! maximum infiltration rate
    zro_inf        , & ! surface runoff
    zro_sfak       , & ! utility variable 
    zro_gfak       , & ! utility variable
    zfmb_fak       , & ! utility variable
    zdwg           , & ! preliminary change of soil water content
    zwgn           , & ! preliminary change of soil water content
    zredfu         , & ! utility variable for runoff determination
    zro            , & ! utility variable for runoff determination
    zro2           , & ! utility variable for runoff determination
    zkorr          , & ! utility variable for runoff determination




!   Implicit solution of  hydraulic equations

    zzz            , & ! utility variable
    z1dgam1        , & ! utility variable
    zgam2m05       , & ! utility variable
    zgam2p05           ! utility variable

  REAL    (KIND=ireals   ) ::  &

!   Statement functions

    zsf_heav       , & ! Statement function: Heaviside function
    zstx            ! dummy argument for Stmt. function

  REAL    (KIND=ireals   ) ::  &

!   Water transport

    zlw_fr_ksom1   , & ! fractional liquid water content of actual layer - 1
    zlw_fr_ksom1_new,& ! fractional liquid water content of actual layer -1
    zlw_fr_kso     , & ! fractional liquid water content of actual layer
    zlw_fr_kso_new , & ! fractional liquid water content of actual layer
    zlw_fr_ksop1   , & ! fractional liquid water content of actual layer + 1
    zlw_fr_ksom05  , & ! fractional liquid water content of actual layer - 1/2
    zdlw_fr_ksom05 , & ! hydraulic diffusivity coefficient at half level above
    zklw_fr_ksom05 , & ! hydraulic conductivity coefficient at half level above
    zklw_fr_kso_new, & ! hydraulic conductivity coefficient at main level
                       !    for actual runoff_g
    zlw_fr_ksop05  , & ! fractional liquid water content of actual layer + 1/2
    zdlw_fr_ksop05 , & ! hydraulic diffusivity coefficient at half level below
    zklw_fr_ksop05     ! hydraulic conductivity coefficient at half level below
   

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND = ireals) :: &

! Model geometry

  zmls     (ke_soil+1)  , & ! depth of soil main level
  zzhls    (ke_soil+1)  , & ! depth of the half level soil layers in m
  zdzhs    (ke_soil+1)  , & ! layer thickness between half levels
  zdzms    (ke_soil+1)      ! distance between main levels

  INTEGER  (KIND=iintegers ) ::  &
  m_styp   (ie,je)          ! soil type

  REAL    (KIND=ireals   ) ::  &

! Connection to the atmosphere

  zrr      (ie,je)      , & ! total rain rate including formation of dew

! Tendencies

  zdwgdt   (ie,je,ke_soil)  ! tendency of water content [kg/(m**2 s)]

  REAL    (KIND=ireals   ) ::  &

! Soil and plant parameters

  zfcap    (ie,je)      , & ! field capacity of soil
  zadp     (ie,je)      , & ! air dryness point
  zporv    (ie,je)      , & ! pore volume (fraction of volume)
  zdw      (ie,je)      , & ! hydrological diff.parameter
  zdw1     (ie,je)      , & ! hydrological diff.parameter
  zkw      (ie,je)      , & ! hydrological cond.parameter
  zkw1     (ie,je)      , & ! hydrological cond.parameter
  zik2     (ie,je)      , & ! minimum infiltration rate
  zrock    (ie,je)      , & ! ice/rock-indicator: 0 for ice and rock
! Hydraulic variables
  zw_fr    (ie,je,ke_soil+1),&!fractional total water content of soil layers
  zinfil   (ie,je)      , & ! infiltration rate 
  zflmg    (ie,je,ke_soil+1), & ! flux of water at soil layer interfaces
  zrunoff_grav(ie,je,ke_soil+1), & ! main level water gravitation
                                !    flux (for runoff calc.)
  runoff_s(ie,je)    ,   & ! surface water runoff; sum over forecast      (kg/m2)
  runoff_g(ie,je)          ! soil water runoff; sum over forecast         (kg/m2)

  REAL    (KIND=ireals   ) ::  &

  zaga(ie,je,1:ke_soil+1),& ! utility variable
  zagb(ie,je,1:ke_soil+1),& ! utility variable
  zagc(ie,je,1:ke_soil+1),& ! utility variable
  zagd(ie,je,1:ke_soil+1),& ! utility variable
  zage(ie,je,1:ke_soil+1)   ! utility variable

  ! ground water as lower boundary of soil column
  REAL    (KIND=ireals) ::  &
    zdelta_sm, zdhydcond_dlwfr, zklw_fr_kso, zklw_fr_ksom1, zdlw_fr_kso 

REAL    (KIND=ireals   ) ::  &
  w_so(ie,je, ke_soil+1, nsteps+2) , &! total water content of each layer and each time step
! output variables
  time_out(noutput+1)                  , &
  w_so_out(ie,je,ke_soil+1, noutput+1) , & 
  wl_so(ie,je,ke_soil+1) , & 
  wl_so_out(ie,je,ke_soil+1, noutput+1) , & 
  s_so_out(ie,je,ke_soil+1, noutput+1) , & 
  rain_out(ie,je,noutput+1)             , &
  raino(ie,je)                          , &
  evapo(ie,je)                          , &
  dwdto(ie,je,ke_soil+1)                , &
  evap_out(ie,je,noutput+1)             , &
  dwdt_out(ie,je,ke_soil+1,noutput+1)   , &
  runoff_s_out(ie,je,noutput+1)        , &
  runoff_g_out(ie,je,noutput+1)        , &
  flmg_out(ie,je,ke_soil+1,noutput+1)    ! flux

real (KIND=ireals   ) evap_max

! LINDA, netcdf output
integer  ncid, iret
integer  timedim, zdim, zhdim
integer  timeid, zid, zhid, wsoid, wlsoid, ssoid, dwid, rgid, rsid, fid
integer  rid, eid
!- End of header
!==============================================================================

! Declaration of STATEMENT-FUNCTIONS

  zsf_heav     (zstx                    ) = 0.5_ireals+SIGN( 0.5_ireals, zstx )

!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------
  CALL init_fields
! select timelevel and timestep for calculations
  zdt = dt

  ! time step for run-off computation
  zroffdt = zdt

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar
  jstarts = jstartpar
  jends   = jendpar


! Computation of derived constants

  zrhwddt = rho_w/zdt     ! density of liquid water/timestep
  zdtdrhw = zdt/rho_w     ! timestep/density of liquid water

! grids for temperature and water content

  zzhls(1) = 2._ireals*czmls(1)   !depth of first half level
  zdzhs(1) = zzhls(1)      !layer thickness betw. half levels of uppermost layer
  zmls(1)  = czmls(1)      !depth of 1st main level
  zdzms(1) = czmls(1)      !layer thickness between soil surface and main level
                           ! of uppermost layer
 
  print*,1,zzhls(1),zdzhs(1),zmls(1),zdzms(1), czmls(1)
  DO kso = 2,ke_soil+1
    zzhls(kso)  = zzhls(kso-1) + 2._ireals*(czmls(kso) -zzhls(kso-1))
    zdzhs(kso) = zzhls(kso) - zzhls(kso-1) ! layer thickness betw. half levels
    zmls(kso)  = czmls(kso)                ! depth of main levels
    zdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
    print*,kso,zzhls(kso),zdzhs(kso),zmls(kso),zdzms(kso), czmls(kso)
  END DO


! Prepare basic surface properties (for land-points only)
 
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF(llandmask(i,j)) THEN        ! for land-points only
        mstyp        = NINT(soiltyp(1,1))        ! soil type
        m_styp(i,j)  = mstyp                     ! array for soil type
        zdw   (i,j)  = 0.0_ireals!cdw0(mstyp)
        zdw1  (i,j)  = cdw1(mstyp)
        zkw   (i,j)  = ckw0(mstyp)
        zkw1  (i,j)  = ckw1(mstyp)
        zik2  (i,j)  = cik2(mstyp)
        zporv(i,j)   = cporv(mstyp)              ! pore volume
        zadp (i,j)   = cadp(mstyp)               ! air dryness point
        zfcap(i,j)   = cfcap(mstyp)              ! field capacity
        zrock(i,j)   = crock(mstyp)              ! rock or ice indicator
        rootdp(i,j)  = czrootdp
      END IF
    END DO
  END DO
  
!   Provide for a soil moisture 1 % above air dryness point, reset soil
!   moisture to zero in case of ice and rock
    DO kso   = 1,ke_soil+1
      DO   j = jstarts, jends
        DO i = istarts, iends
          IF (llandmask(i,j)) THEN             ! for land-points only
            IF (m_styp(i,j).ge.3) THEN
!              w_so (i,j,kso,1) = MAX(w_so_init(i,j,kso),                     &
!                                           1.01_ireals*zadp(i,j)*zdzhs(kso) )
              w_so (i,j,kso,1) = 0.5_ireals*zporv(i,j)*zdzhs(kso)
            ELSE
              w_so(i,j,kso,1) = 0.0_ireals
            ENDIF
         END IF   ! land-points
        END DO
      END DO
    END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
        zw_fr(i,j,ke_soil+1)  = w_so(i,j,ke_soil+1,1)/zdzhs(ke_soil+1)
    END DO
  END DO

  DO kso   = 1, ke_soil
    DO   j = jstarts, jends
      DO i = istarts, iends
        zw_fr(i,j,kso)    = w_so(i,j,kso,1)/zdzhs(kso)
      END DO
    END DO
  END DO
!------------------------------------------------------------------------------
! Section II.1: Initializations
!------------------------------------------------------------------------------

  ! Number of soil layers contributing to surface run-off
    msr_off  = 0
    DO   j = jstarts, jends
      DO i = istarts, iends
    runoff_s(i,j) = 0.0_ireals   
    runoff_g(i,j) = 0.0_ireals
      END DO
    END DO
 
! calculate active layers
        czbot_w_so = MIN(czbot_w_so,zzhls(ke_soil))
        ibot_w_so = ke_soil
   DO kso = 1,ke_soil,1
        IF(zzhls(kso) <= czbot_w_so) THEN
               ibot_w_so = kso
        END IF
   END DO
  ke_soil_hy = max(ibot_w_so,2) 
  print*,ke_soil_hy
  DO kso = 1,ke_soil,1
     IF(zzhls(kso) <= czrootdp) THEN
        irootdp = kso
     END IF
  END DO
  time = starttime

! Time stepping:
!-----------------------------------------------------------------------------
    DO nx = 1,nsteps+1,1 ! Beginning time step loop
       ! change next timestep index and current time
        nnew = nx + 1
        time = starttime + (nx-1)*zdt
         
      IF(mod(time,outtime) == 0) THEN
         nout = INT(time/outtime) + 1
            ! Prepare output
         time_out(nout) = time
            ! Precipitation is provided at each output time step therefore change zrr only at time steps which are also output time steps          
         DO   j = jstarts, jends
            DO i = istarts, iends
               zrr(i,j) =  (prec(nout)/rho_w)*(rho_w/zdt) 
! precipitation file in mm but need precipitation rate in kg/m^2*s
            END DO
         END DO
      END IF ! end of output/input preparation


  
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF(llandmask(i,j))THEN     ! land-points only
        ! infiltration and surface run-off
        ! maximum infiltration rate of the soil (rock/ice/water-exclusion
        zinfmx = zrock(i,j)*csvoro &
                 *( cik1*MAX(0.5_ireals,plcov(i,j))*MAX(0.0_ireals,           &
                 zporv(i,j)-zw_fr(i,j,1))/zporv(i,j) + zik2(i,j) )

        ! to avoid pore volume water excess of the uppermost layer by 
        ! infiltration
        zinfmx = MIN(zinfmx, (zporv(i,j) - zw_fr(i,j,1))*zdzhs(1)*zrhwddt)
        ! final infiltration rate limited by maximum value
        zinfil(i,j) = MIN(zinfmx,zrr(i,j))
        ! surface run-off (residual of potential minus actual infiltration)
        zro_inf       = zrr(i,j) - zinfil(i,j)
        runoff_s(i,j) = runoff_s(i,j) + zro_inf*outtime
        raino(i,j) = raino(i,j)+zrr(i,j)*zdt ! in mm
      END IF            ! land-points only
    END DO
  END DO

!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------

! uppermost layer, kso = 1
DO   j = jstarts, jends
  DO i = istarts, iends
    
    IF (llandmask(i,j)) THEN      ! land-points only
      IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
        zlw_fr_kso  = zw_fr(i,j,1)
        zlw_fr_ksop1= zw_fr(i,j,2)

        
        ! interpolated scaled liquid water fraction at layer interface
        zlw_fr_ksop05 = 0.5_ireals*(zdzhs(2)*zlw_fr_kso+zdzhs(1)*zlw_fr_ksop1) &
                                             /zdzms(2)
        wl_so(i,j,1) = zlw_fr_ksop05/zporv(i,j)
        zdlw_fr_ksop05 = zdw(i,j)*EXP(zdw1(i,j)*                       &
                         (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
        zklw_fr_ksop05 = zkw(i,j)*EXP(zkw1(i,j)*                       &
                         (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )


        ! coefficients for implicit flux computation
        z1dgam1 = zdt/zdzhs(1) 
        zgam2p05 = zdlw_fr_ksop05/zdzms(2)
        zaga(i,j,1) = 0._ireals
        zagb(i,j,1) = 1._ireals+zalfa*zgam2p05*z1dgam1
        zagc(i,j,1) = -zalfa * zgam2p05*z1dgam1
        zagd(i,j,1) = zw_fr(i,j,1) + zinfil(i,j)*z1dgam1/rho_w  &
                     -zklw_fr_ksop05*z1dgam1                     &
                     +(1. - zalfa)* zgam2p05*z1dgam1*(zlw_fr_ksop1 - zlw_fr_kso)

        ! explicit part of soil surface water flux:
        zflmg (i,j,1) = - zinfil(i,j)! boundary value for soil water transport
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
          zlw_fr_ksom1  = zw_fr(i,j,kso-1)
          zlw_fr_kso    = zw_fr(i,j,kso  )
          zlw_fr_ksop1  = zw_fr(i,j,kso+1)
          ! interpolated scaled liquid water content at interface to layer
          ! above and below

          zlw_fr_ksom05 = 0.5*(zdzhs(kso-1)*zlw_fr_kso+   &
                                 zdzhs(kso)*zlw_fr_ksom1)/zdzms(kso)
          zlw_fr_ksop05 = 0.5*(zdzhs(kso+1)*zlw_fr_kso+   &
                                 zdzhs(kso)*zlw_fr_ksop1)/zdzms(kso+1)
          wl_so(i,j,kso) = zlw_fr_ksop05/zporv(i,j)

          zdlw_fr_ksom05= zdw(i,j)*EXP( zdw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zdlw_fr_ksop05= zdw(i,j)*EXP( zdw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
          zklw_fr_ksom05= zkw(i,j)*EXP( zkw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zklw_fr_ksop05= zkw(i,j)*EXP( zkw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )

          ! coefficients for implicit flux computation
          z1dgam1 = zdt/zdzhs(kso)
          zgam2m05  = zdlw_fr_ksom05/zdzms(kso)
          zgam2p05  = zdlw_fr_ksop05/zdzms(kso+1)

          !if(zlw_fr_ksop05.lt.zporv(i,j)-0.01) then

             zaga (i,j,kso) = -zalfa*zgam2m05*z1dgam1
             zagc (i,j,kso) = -zalfa*zgam2p05*z1dgam1
             zagb (i,j,kso) = 1._ireals +zalfa*(zgam2m05+zgam2p05)*z1dgam1
             zagd (i,j,kso) = zw_fr(i,j,kso)+                               &
                                z1dgam1*(-zklw_fr_ksop05+zklw_fr_ksom05)+ &
                                (1._ireals-zalfa)*z1dgam1*                &
                                (zgam2p05*(zlw_fr_ksop1-zlw_fr_kso  )     &
                                -zgam2m05*(zlw_fr_kso  -zlw_fr_ksom1)   )
             !soil water flux, explicit part, for soil water flux investigations
             ! only)
             zflmg(i,j,kso) = rho_w*(zdlw_fr_ksom05*(zlw_fr_kso-zlw_fr_ksom1)/  &
                  zdzms(kso) - zklw_fr_ksom05)
          !else
          !   zaga(i,j,kso) = -zalfa* zgam2m05*z1dgam1
          !   zagb(i,j,kso) = 1.+ zalfa*zgam2m05*z1dgam1
          !   zagc(i,j,kso) = 0.0
          !   zagd(i,j,kso) = zw_fr(i,j,kso)+z1dgam1*zklw_fr_ksom05 &
          !        +(1.-zalfa)*z1dgam1*                              &
          !        zgam2m05*(zlw_fr_ksom1  - zlw_fr_kso)
          !   zflmg(i,j,kso)= 0.0
          !end if

          IF(kso==ke_soil_hy-1) THEN
            zflmg(i,j,kso+1)=rho_w*(zdlw_fr_ksop05*(zlw_fr_ksop1-zlw_fr_kso)/ &
                            zdzms(kso+1)  - zklw_fr_ksop05)
          ENDIF
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
          zlw_fr_ksom1  = zw_fr(i,j,ke_soil_hy-1)
          zlw_fr_kso    = zw_fr(i,j,ke_soil_hy  )
          zlw_fr_ksom05 = 0.5*(zdzhs(ke_soil_hy-1)*zlw_fr_kso+ &
                              zdzhs(ke_soil_hy)*zlw_fr_ksom1)/zdzms(ke_soil_hy)

          
          zdlw_fr_ksom05= zdw(i,j)*EXP( zdw1(i,j)* &
                            (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )

          z1dgam1 = zdt/zdzhs(ke_soil_hy)
          zgam2m05  = zdlw_fr_ksom05/zdzms(ke_soil_hy)
          zklw_fr_ksom05= zkw(i,j)*EXP( zkw1(i,j)* &
                            (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zaga(i,j,ke_soil_hy) = -zalfa* zgam2m05*z1dgam1
          zagb(i,j,ke_soil_hy) = 1.+ zalfa*zgam2m05*z1dgam1
          zagc(i,j,ke_soil_hy) = 0.0
          zagd(i,j,ke_soil_hy) = zw_fr(i,j,ke_soil_hy)+z1dgam1*zklw_fr_ksom05 &
                            +(1.-zalfa)*z1dgam1*                              &
                             zgam2m05*(zlw_fr_ksom1  - zlw_fr_kso)
          !zflmg(i,j,ke_soil_hy)= 0.0
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
          ! compute implicit part of new liquid water content  
          w_so(i,j,kso,nnew) = zage(i,j,kso)*zdzhs(kso)
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
         ! water content
         w_so(i,j,ke_soil_hy,nnew) = zage(i,j,ke_soil_hy)*zdzhs(ke_soil_hy)  
       END IF 
     END IF          ! land-points only
   END DO
 END DO
! to ensure vertical constant water concentration profile beginning at 
! layer ke_soil_hy for energetic treatment only
! soil water climate layer(s)


  DO kso = ke_soil_hy+1,ke_soil+1
   DO   j = jstarts, jends
     DO i = istarts, iends
       IF (llandmask(i,j)) THEN          ! land-points only
         IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
            IF (itype_hydbound == 3) THEN
! Constant saturated lower boundary condition: (uncomment line below and comment other boundary condition below)
               w_so(i,j,kso,nnew) = zporv(i,j)*zdzhs(kso)
            ELSE
! No flux Boundary condition i.e. constant water concentration profile at ke_soil_hy: (uncomment line below and comment other boundary condition above)
               w_so(i,j,kso,nnew) = w_so(i,j,kso-1,nnew)*zdzhs(kso)/zdzhs(kso-1) 
            END  IF
         END IF
       END IF          ! land-points only
     END DO
   END DO
  END DO


  DO kso = ke_soil_hy,1,-1
     DO   j = jstarts, jends
        DO i = istarts, iends
           zlw_fr_kso_new  = w_so(i,j,kso  ,nnew)/zdzhs(kso  )
           IF (zlw_fr_kso_new>zporv(i,j)) THEN
              w_so(i,j,kso-1,nnew) = w_so(i,j,kso-1,nnew)+(zlw_fr_kso_new-zporv(i,j))*zdzhs(kso-1)
              w_so(i,j,kso  ,nnew) = zporv(i,j)*zdzhs(kso)
           END IF
        END DO
     END DO
  END DO
! combine implicit part of sedimentation and capillary flux with explicit part
! (for soil water flux investigations only)
DO kso = 2,ke_soil+1
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN          ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
          zlw_fr_ksom1_new= w_so(i,j,kso-1,nnew)/zdzhs(kso-1)
          zlw_fr_kso_new  = w_so(i,j,kso  ,nnew)/zdzhs(kso  )
          zlw_fr_ksom1  = w_so(i,j,kso-1,nx)/zdzhs(kso-1)
          zlw_fr_kso    = w_so(i,j,kso  ,nx)/zdzhs(kso  )
          !... additionally for runoff_g at lower level of lowest active water
          ! layer calculated with (upstream) main level soil water content
          ! compute reduction factor for transport coefficients
          
        ! interpolated liquid water content at interface to layer above
          zlw_fr_ksom05 =0.5*(zdzhs(kso)*zlw_fr_ksom1+zdzhs(kso-1)*zlw_fr_kso) &
                            /zdzms(kso)
          zdlw_fr_ksom05= zdw(i,j)*EXP(zdw1(i,j) *  &
                          (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
          zklw_fr_ksom05= zkw(i,j) * EXP(zkw1(i,j)* &
                          (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )

          IF (kso> ke_soil_hy) zdlw_fr_ksom05=0.0   ! no flux gradient 
                                                    ! contribution below 2.5m
          IF (kso> ke_soil_hy) zklw_fr_ksom05=0.0   ! no gravitation flux below 2.5m
          
          zflmg(i,j,kso) =                              &
                   (1._ireals-zalfa) * zflmg(i,j,kso) + & ! explicit flux component
                              zalfa  * rho_w *          & ! implicit flux component
                   (zdlw_fr_ksom05 * (zlw_fr_kso_new-zlw_fr_ksom1_new)/zdzms(kso) &
                    - zklw_fr_ksom05)

          zklw_fr_kso_new = zkw(i,j) * EXP(zkw1(i,j)* &
                            (zporv(i,j) - zlw_fr_kso_new)/(zporv(i,j) - zadp(i,j)) )
          ! actual gravitation water flux
          IF(w_so (i,j,kso,nnew).LT.1.01_ireals*zadp(i,j)*zdzhs(kso)) zklw_fr_kso_new = 0._ireals
         
          zrunoff_grav(i,j,kso) =  - rho_w * zklw_fr_kso_new

          ! ground water as lower boundary of soil column
          IF ((kso == ke_soil_hy+1).and.(itype_hydbound == 3)) THEN
             zdelta_sm=( zlw_fr_kso_new - zlw_fr_ksom1_new )

             zdlw_fr_kso = zdw(i,j)*EXP(zdw1(i,j) *  &
                  (zporv(i,j)-zlw_fr_kso_new)/(zporv(i,j)-zadp(i,j)) )
             zklw_fr_kso = zkw(i,j) * EXP(zkw1(i,j)* &
                  (zporv(i,j)-zlw_fr_kso_new)/(zporv(i,j)-zadp(i,j)) )
             zklw_fr_ksom1 = zkw(i,j) * EXP(zkw1(i,j)* &
                  (zporv(i,j)-zlw_fr_ksom1_new)/(zporv(i,j)-zadp(i,j)) )

             zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
             zrunoff_grav(i,j,ke_soil_hy)=zrunoff_grav(i,j,ke_soil_hy)+ &
                  zdhydcond_dlwfr / &
                  (1.0-exp(-zdhydcond_dlwfr/zdlw_fr_kso*0.5*zdzms(ke_soil_hy+1)))* &
                  zdelta_sm
          ENDIF
        END IF 
      END IF          ! land-points only
    END DO
  END DO
END DO

  DO  kso = 1,ke_soil
    ! utility variables used to avoid if-constructs in following loops
    zro_sfak = zsf_heav(0.5_ireals + msr_off - kso)  ! 1.0 for 'surface runoff'
    zro_gfak = 1._ireals - zro_sfak                  ! 1.0 for 'ground runoff'

    IF (kso==ibot_w_so) THEN
      zfmb_fak = 1.0_ireals
      !zfmb_fak = 0.0_ireals
    ELSE
      zfmb_fak = 0.0_ireals
    END IF


    DO   j = jstarts, jends
      DO i = istarts, iends
        IF (llandmask(i,j)) THEN      ! land-points only
          ! sedimentation and capillary transport in soil
          IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type

            ! first runoff calculation without consideration of
            ! evapotranspiration
            !zdwg   =  zflmg(i,j,kso+1) - zflmg(i,j,kso)
            !zdwg calculated above by flux divergence has to be aequivalent with 
            zdwg =  (w_so(i,j,kso,nnew)/zdzhs(kso)-zw_fr(i,j,kso))*zdzhs(kso) &
                                                                   /zdtdrhw
            if (kso.le.irootdp) then
               evap_max   = w_so(i,j,kso,nnew)/zdt*irootdp
               evap(nout) = min(evap(nout),evap_max)
               evapo(i,j) = evapo(i,j)-evap(nout)/irootdp*zdt*rho_w ! convert to mm
!                 dwdto(i,j,kso)     =  dwdto(i,j,kso) - zdtdrhw*(evap(nout)/irootdp)*zdzhs(kso)&
!                      -zdtdrhw*q_kso(i,j,kso)
               zdwg = zdwg - rho_w*(evap(nout)/irootdp)
            end if
            zdwg =  zdwg + zrunoff_grav(i,j,kso)*zfmb_fak
            !zredfu =  MAX( 0.0_ireals, MIN( 1.0_ireals,(zw_fr(i,j,kso) -     &
            !           zfcap(i,j))/MAX(zporv(i,j) - zfcap(i,j),zepsi)) )
            zredfu =  MAX( 0.0_ireals, MIN( 1.0_ireals,(zw_fr(i,j,kso) -     &
                       zporv(i,j))/MAX(zporv(i,j) - zporv(i,j),zepsi)) )
            zredfu = zsf_heav(zdwg)*zredfu
            zro    = zdwg*zredfu
            zdwg   = zdwg*(1._ireals - zredfu)
            zwgn   = zw_fr(i,j,kso) + zdtdrhw*zdwg/zdzhs(kso)
            zro2   = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zwgn - zporv(i,j))
            zkorr  = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zadp(i,j) - zwgn )
            zdwgdt(i,j,kso)= zdwg + zkorr - zro2
            zro    = zro      + zro2
            runoff_s(i,j) = runoff_s(i,j) + zro*zro_sfak*zroffdt
            runoff_g(i,j) = runoff_g(i,j) + zro*zro_gfak*zroffdt
            ! runoff_g reformulation:
            runoff_g(i,j) = runoff_g(i,j) - (zrunoff_grav(i,j,kso) * zfmb_fak &
                                          + zkorr) * zroffdt

          END IF          ! ice/rock-exclusion
        END IF   ! land-points only
      END DO
    END DO
  END DO         ! end loop over soil layers

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

DO kso = 1,ke_soil
  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN  ! for landpoints only
        w_so(i,j,kso,nnew) = w_so(i,j,kso,nx) + zdt*zdwgdt(i,j,kso)/rho_w
        zw_fr(i,j,kso)    = w_so(i,j,kso,nnew)/zdzhs(kso)
        dwdto(i,j,kso)    = dwdto(i,j,kso) + (w_so(i,j,kso,nnew) - w_so(i,j,kso,nnew-1))*rho_w
      END IF  ! land-points only
    END DO
  END DO
END DO        ! soil layers


 IF(mod(time,outtime) == 0) THEN   
 ! Prepare output
           DO j = jstarts, jends
              DO i = istarts, iends
                 rain_out(i,j,nout) = raino(i,j)
                 raino(i,j) = 0.0_ireals
                 evap_out(i,j,nout) = evapo(i,j)
                 evapo(i,j) = 0.0_ireals
                 runoff_s_out(i,j,nout) = runoff_s(i,j)
                 runoff_s(i,j) = 0.0_ireals 
                 runoff_g_out(i,j,nout) = runoff_g(i,j)
                 runoff_g(i,j) = 0.0_ireals
                 DO kso=1,ke_soil+1
                    wl_so_out(i,j,kso,nout) = wl_so(i,j,kso)
                    w_so_out(i,j,kso,nout) = w_so(i,j,kso,nnew)
                    s_so_out(i,j,kso,nout) = w_so(i,j,kso,nnew)/(zporv(i,j)*zdzhs(kso))
                    dwdt_out(i,j,kso,nout) = dwdto(i,j,kso)
                    dwdto(i,j,kso) = 0.0_ireals
                 END DO
                 DO kso=1,ke_soil
                    flmg_out(i,j,kso,nout) = zflmg(i,j,kso+1)
                 END DO
                 flmg_out(i,j,ke_soil+1,nout) = 0.0_ireals
              END DO
           END DO

  END IF ! end of output preparation



    END DO ! END of time stepping

!------------------------------------------------------------------------------
! Output into files
!------------------------------------------------------------------------------
! WRITE output into file
 OPEN(UNIT=12, FILE="w_so_out.dat", ACTION="write", STATUS="replace")
 OPEN(UNIT=13, FILE="runoff.dat", ACTION="write", STATUS="replace")
DO j = jstarts, jends
   DO i = istarts, iends
      DO nx = 1,noutput+1,1
   ! Write water content to file output.dat
    WRITE(12,"(F10.0)", advance="no") time_out(nx)
    WRITE(12,"(A1)", advance="no") "  "
         DO kso = 1,ke_soil+1,1
           WRITE(12,"(10F10.8)", advance="no") w_so_out(i,j,kso,nx)
           WRITE(12,"(A2)", advance="no") "  "
         END DO
       WRITE(12,*)
! Write both runoffs to file runoff.dat
       WRITE(13,"(F10.0)", advance="no") time_out(nx)
       WRITE(13,"(A1)", advance="no") " "
       WRITE(13,"(F10.8)", advance="no") runoff_s_out(i,j,nx)
       WRITE(13,"(A1)", advance="no") " "
       WRITE(13,"(F10.8)") runoff_g_out(i,j,nx)

      END DO
    END DO
END DO
  
 CLOSE(UNIT=12)
 CLOSE(UNIT=13)

 !LINDA, write to NetCDF file

 iret = nf90_create ("w_so_out.nc",  NF90_CLOBBER, ncid)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)
        
 !* define dimensions

 iret = nf90_def_dim(ncid, 'z_soil', ke_soil+1, zdim)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 iret = nf90_def_dim(ncid, 'zm_soil', ke_soil+1, zhdim)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 iret = nf90_def_dim(ncid, 'time', NF90_UNLIMITED,timedim)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 iret = nf90_def_var(ncid, 'time', NF90_REAL, timedim, timeid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'z_soil', NF90_REAL, zdim, zid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'zm_soil', NF90_REAL, zhdim, zhid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'W_SO', NF90_REAL, (/zdim, timedim/), wsoid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'WL_SO', NF90_REAL, (/zhdim, timedim/), wlsoid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'S_SO', NF90_REAL, (/zdim, timedim/), ssoid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'DWDT', NF90_REAL, (/zdim, timedim/), dwid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'RUNOFF_S', NF90_REAL, (/timedim/), rsid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'RAIN', NF90_REAL, (/timedim/), rid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'EVAP', NF90_REAL, (/timedim/), eid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'RUNOFF_G', NF90_REAL, (/timedim/), rgid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'FLX', NF90_REAL, (/zhdim, timedim/), fid)
 call check_err(iret)

 iret = nf90_put_att(ncid, timeid,'standard_name',"time")
 call check_err(iret)
 iret = nf90_put_att(ncid, timeid,'long_name',"time")
 call check_err(iret)
 iret = nf90_put_att(ncid, timeid,'units',"seconds since 2006-07-12 00:00:00")
 call check_err(iret)
 iret = nf90_put_att(ncid, timeid,'calendar',"proleptic_gregorian")
 call check_err(iret)

 iret = nf90_put_att(ncid, zid,'standard_name',"depth of soil layer")
 call check_err(iret)
 iret = nf90_put_att(ncid, zid,'long_name',"depth of soil layer")
 call check_err(iret)
 iret = nf90_put_att(ncid, zid,'units',"m")
 call check_err(iret)

 iret = nf90_put_att(ncid, wsoid,'standard_name',"soil water content")
 call check_err(iret)
 iret = nf90_put_att(ncid, wsoid,'long_name',"soil water content")
 call check_err(iret)
 iret = nf90_put_att(ncid, wsoid,'units',"m")
 call check_err(iret)

 iret = nf90_put_att(ncid, ssoid,'standard_name',"soil saturation")
 call check_err(iret)
 iret = nf90_put_att(ncid, ssoid,'long_name',"soil saturation")
 call check_err(iret)
 iret = nf90_put_att(ncid, ssoid,'units',"")
 call check_err(iret)

 iret = nf90_put_att(ncid, rgid,'standard_name',"ground runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rgid,'long_name',"ground water runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rgid,'units',"m/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, rsid,'standard_name',"surface runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rsid,'long_name',"surface water runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rsid,'units',"m/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, fid,'standard_name',"soil-water flux")
 call check_err(iret)
 iret = nf90_put_att(ncid, fid,'long_name',"flux of soil water")
 call check_err(iret)
 iret = nf90_put_att(ncid, fid,'units',"m/s")
 call check_err(iret)

 !* leave define mode
 IRET = NF90_ENDDEF(ncid)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, timeid,time_out)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, zid,zmls)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, zhid,zzhls)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, wsoid, w_so_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, wlsoid, wl_so_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, ssoid, s_so_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, dwid, dwdt_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, rid, rain_out(1,1,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, eid, evap_out(1,1,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, rsid, runoff_s_out(1,1,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, rgid, runoff_g_out(1,1,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, fid, flmg_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_CLOSE(NCID)
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 contains
 subroutine check_err(iret)

   implicit none
   integer iret
   if (iret .ne. NF90_NOERR) then
      print *, nf90_strerror(iret)
      stop
   endif
 end subroutine check_err

END PROGRAM soil_model
