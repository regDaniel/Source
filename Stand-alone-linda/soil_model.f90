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
     !w_so_init      ,     & ! initial total water conent (liquid water)       (m H20)
     prec           ,     & ! precipitation rate 
     evap           ,     & ! evaporation rate 
     ldecharme      ,     & ! Decharme formulation
     s_topo         ,     &
! tuning parameters
!-----------------------------------------------------------------------------
     gamma          ,     & !Scaling parameter for linear runoff-orography scaling
     kexpdec        ,     & !Decharme factor in the exponential function
    
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
     !zredfu         , & ! utility variable for runoff determination
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
     zstx               ! dummy argument for Stmt. function

REAL    (KIND=ireals   ) ::  &

!   Water transport

     zlw_fr_ksom1   , & ! fractional liquid water content of actual layer - 1
     !zlw_fr_ksom1_new,& ! fractional liquid water content of actual layer -1
     !zlw_fr_kso_new , & ! fractional liquid water content of actual layer
     zlw_fr_ksop1   , & ! fractional liquid water content of actual layer + 1
     zlw_fr_ksom05  , & ! fractional liquid water content of actual layer - 1/2
     zlw_fr_ksop05  , & ! fractional liquid water content of actual layer + 1/2
!     zdlw_fr_ksom05 , & ! hydraulic diffusivity coefficient at half level above
!     zklw_fr_ksom05 , & ! hydraulic conductivity coefficient at half level above
!     zklw_fr_kso_new, & ! hydraulic conductivity coefficient at main level
                        !    for actual runoff_g

     zlw_fr_kso(ie,je,1:ke_soil+1)     , &  ! fractional liquid water content of actual layer
     zdlw_fr_ksop05(ie,je,1:ke_soil+1) , &  ! hydraulic diffusivity coefficient at half level below
     zklw_fr_ksop05(ie,je,1:ke_soil+1) , &  ! hydraulic conductivity coefficient at half level below
     zdlw_fr_ksop05_fg(ie,je,1:ke_soil+1),& ! first guess hydraulic diffusivity coefficient at half level below
     zklw_fr_ksop05_fg(ie,je,1:ke_soil+1),& ! first guess hydraulic conductivity coefficient at half level below
     zklw_fr_ksop05_sg(ie,je,1:ke_soil+1),& ! second guess hydraulic conductivity coefficient at half level below
     h_wt_k(ie,je,1:ke_soil+1)           ,& ! portion of the box that is saturated
     q_kso(ie,je,1:ke_soil+1)            ,& ! ground water flow, aquifer
     q_kso_fc(ie,je,1:ke_soil+1)            ! runoff, field capacity

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

       zdwgdt   (ie,je,ke_soil)  ! tendency of water content [kg/(m**3 s)]

   REAL    (KIND=ireals   ) ::  &

! Soil and plant parameters

     zfcap    (ie,je)      , & ! field capacity of soil
     zadp     (ie,je)      , & ! air dryness point
     zporv    (ie,je)      , & ! pore volume (fraction of volume)
     zdw      (ie,je)      , & ! hydrological diff.parameter
     zdw1     (ie,je)      , & ! hydrological diff.parameter
     zkw      (ie,je,ke_soil+1), & ! hydrological cond.parameter
     zkw1     (ie,je)      , & ! hydrological cond.parameter
     zik2     (ie,je)      , & ! minimum infiltration rate
     zrock    (ie,je)      , & ! ice/rock-indicator: 0 for ice and rock
! Hydraulic variables
     zw_fr    (ie,je,ke_soil+1)   ,&!fractional total water content of soil layers
     zinfil   (ie,je)             , & ! infiltration rate 
     zflmg    (ie,je,ke_soil+1)   , & ! flux of water at soil layer interfaces
     !zrunoff_grav(ie,je,ke_soil+1), & ! main level water gravitation
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
!   REAL    (KIND=ireals) ::  &
     !zdelta_sm, zdhydcond_dlwfr, zklw_fr_kso, zklw_fr_ksom1, zdlw_fr_kso 

   REAL    (KIND=ireals   ) ::  &
     w_so(ie,je, ke_soil+1, nsteps+2) , &! total water content of each layer and each time step
! output variables
     time_out(noutput+1)                   , &
     w_so_out(ie,je,ke_soil+1, noutput+1)  , & 
     wl_so_out(ie,je,ke_soil+1, noutput+1) , & 
     s_so_out(ie,je,ke_soil+1, noutput+1)  , & 
     rain_out(ie,je,noutput+1)             , &
     raino(ie,je)                          , &
     evapo(ie,je)                          , &
     dwdto(ie,je,ke_soil+1)                , &
     evap_out(ie,je,noutput+1)             , &
     dwdt_out(ie,je,ke_soil+1,noutput+1)   , &
     runoff_s_out(ie,je,noutput+1)         , &
     runoff_g_out(ie,je,noutput+1)         , &
     flmg_out(ie,je,ke_soil+1,noutput+1)   , & ! flux
     klwp_out(ie,je,ke_soil+1,noutput+1)    , & ! conductivity
     klwm_out(ie,je,ke_soil+1,noutput+1)    , & ! conductivity
     dlwpf_out(ie,je,ke_soil+1,noutput+1)   , & ! conductivity
     dlwp_out(ie,je,ke_soil+1,noutput+1)    , & ! conductivity
     q_kso_out(ie,je,ke_soil+1,noutput+1)   , & ! runoff from aquifer
     q_kso_fc_out(ie,je,ke_soil+1,noutput+1)    ! runoff field capacity exceeded
! LINDA, netcdf output
   integer  ncid, iret
   integer  timedim, zdim, zhdim
   integer  timeid, zid, zhid, wsoid, wlsoid, ssoid, dwid, rgid, rsid, fid
   integer  slid, rslid, kmid, kpid, dpfid, dpid, qid, qfcid, rid, eid
! LINDA, first saturated level
   integer  zsatlev(ie,je), satlev_out(ie,je,noutput+1), count
   real (KIND=ireals   )     z_soil_hy, z_wt, h_wt , zklw_fr_ksom05_limit
   real (KIND=ireals   )     flux_ksom05,flux_ksom05_pos, flux_ksop05, flux_ksop05_neg, q_kso_pos, alpha_lim, evap_max
   real (KIND=ireals   )     hlp!, q_kso_fc_max
   real (KIND=ireals   )     rsatlev_out(ie,je,noutput+1)
   logical limit
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
 
   DO kso = 2,ke_soil+1
      zzhls(kso)  = zzhls(kso-1) + 2._ireals*(czmls(kso) -zzhls(kso-1))
      zdzhs(kso) = zzhls(kso) - zzhls(kso-1) ! layer thickness betw. half levels
      zmls(kso)  = czmls(kso)                ! depth of main levels
      zdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
   END DO


   ! Prepare basic surface properties (for land-points only)
 
   DO   j = jstarts, jends
      DO i = istarts, iends
         IF(llandmask(i,j)) THEN        ! for land-points only
            mstyp        = NINT(soiltyp(1,1))        ! soil type
            m_styp(i,j)  = mstyp                     ! array for soil type
            zdw   (i,j)  = cdw0(mstyp)
            zdw1  (i,j)  = cdw1(mstyp)
            zkw   (i,j,ke_soil+1) = ckw0  (mstyp)   ! other levels are set below with 3D variables
            !zkw   (i,j)  = ckw0(mstyp)
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

   ! Set three-dimensional variables
   DO kso = 1, ke_soil+1
      DO   j = jstarts, jends
         DO i = istarts, iends
            IF(llandmask(i,j)) THEN        ! for land-points only
               mstyp           = NINT(soiltyp(i,j))        ! soil type
               IF (ldecharme) then
!<JH
    !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
                  zkw   (i,j,kso) = ckw0 (mstyp)*EXP(-kexpdec*(zmls(kso)-rootdp(i,j)))
!>JH
               else
                  zkw   (i,j,kso) = ckw0 (mstyp)
               end if
            ENDIF
         ENDDO
      ENDDO
   ENDDO


   !   Provide for a soil moisture 1 % above air dryness point, reset soil
   !   moisture to zero in case of ice and rock
   DO kso   = 1,ke_soil+1
      DO   j = jstarts, jends
         DO i = istarts, iends
            IF (llandmask(i,j)) THEN             ! for land-points only
               IF (m_styp(i,j).ge.3) THEN
                  q_kso(i,j,kso) = 0.0_ireals
                  q_kso_fc(i,j,kso) = 0.0_ireals
!                  w_so (i,j,kso,1) = MAX(w_so_init(i,j,kso),                     &
!                       1.01_ireals*zadp(i,j)*zdzhs(kso) )
                w_so (i,j,kso,1) = 0.5_ireals*zporv(i,j)*zdzhs(kso)
!                w_so (i,j,kso,1) = 0.95_ireals*zporv(i,j)*zdzhs(kso)
!                w_so (i,j,kso,1) = 1.01_ireals*zadp(i,j)*zdzhs(kso)
               ELSE
                  w_so(i,j,kso,1) = 0.0_ireals
               ENDIF
            END IF   ! land-points
         END DO
      END DO
   END DO
!   w_so (1,1,2,1) = 1.02_ireals*zadp(1,1)*zdzhs(2)
!   w_so (1,1,3,1) = 1.02_ireals*zadp(1,1)*zdzhs(3)
!   w_so (1,1,6,1) = zporv(1,1)*zdzhs(6)
!   w_so (1,1,7,1) = zporv(1,1)*zdzhs(7)

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
         DO kso = 1,ke_soil+1
            dwdto(i,j,kso) = 0.0_ireals
         END DO
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
   DO kso = 1,ke_soil,1
      IF(zzhls(kso) <= czrootdp) THEN
         irootdp = kso
      END IF
   END DO

   ke_soil_hy = max(ibot_w_so,2) 
   z_soil_hy  = zmls(ke_soil_hy+1) ! depth of lower layer bounding the active layers 
   print*,ke_soil_hy, z_soil_hy, zzhls(ke_soil_hy)
   z_wt = zzhls(ke_soil+1)
   h_wt = 0.0_ireals
   print*,'start',z_wt

   time = starttime

   count= 0
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
               zrr(i,j) =  (prec(nout)/1000.0_ireals)*(rho_w)/zdt
               ! precipitation file in mm but need precipitation rate in kg/m^2*s
            END DO
         END DO
      END IF ! end of output/input preparation

      ! determine saturation level
      zsatlev(:,:) = ke_soil_hy
      DO   j = jstarts, jends
         DO i = istarts, iends
               do while (((zw_fr(i,j,zsatlev(i,j))+100.*zepsi).ge.zporv(i,j)).and.(zsatlev(i,j).gt.0))
                  zsatlev(i,j) = zsatlev(i,j)-1
               end do
               zsatlev(i,j) = zsatlev(i,j)+1
         END DO
      END DO
  
      ! diagnose water-table depth
      DO   kso =1,ke_soil+1
         DO   j = jstarts, jends
            DO i = istarts, iends
               zlw_fr_ksom1  = zw_fr(i,j,kso-1)
               zlw_fr_kso(i,j,kso)    = zw_fr(i,j,kso  )
               zlw_fr_ksom05 = 0.5*(zdzhs(kso-1)*zlw_fr_kso(i,j,kso)+ &
                    zdzhs(kso)*zlw_fr_ksom1)/zdzms(kso)
               h_wt_k(i,j,kso) = min(zdzhs(kso),max(0.0,zdzhs(kso)* & 
                    (zlw_fr_kso(i,j,kso)-zfcap(i,j))/zporv(i,j))) ! dz in layer
               if(kso==zsatlev(i,j)-1) then! last partly saturated level above ground water
                  h_wt = min(zdzhs(kso),max(0.0,zdzhs(kso)*(zlw_fr_kso(i,j,kso)-zlw_fr_ksom1)/(zporv(i,j)-zlw_fr_ksom1))) ! dz in layer
                  z_wt = zzhls(kso) - h_wt ! depth of water table
               end if
            END DO
         END DO
      END DO

      DO   j = jstarts, jends
         DO i = istarts, iends
            IF(llandmask(i,j))THEN     ! land-points only
               ! infiltration and surface run-off
               ! maximum infiltration rate of the soil (rock/ice/water-exclusion
               IF (LDECHARME) THEN
                  zinfmx = zrock(i,j)*zkw(i,j,1)*rho_w
               ELSE
                  zinfmx = zrock(i,j)*csvoro &
                       *( cik1*MAX(0.5_ireals,plcov(i,j))*MAX(0.0_ireals,           &
                       zporv(i,j)-zw_fr(i,j,1))/zporv(i,j) + zik2(i,j) )
               END IF
               ! to avoid pore volume water excess of the uppermost layer by 
               ! infiltration
               zlw_fr_ksop05 = 0.5_ireals*(zdzhs(2)*zlw_fr_kso(i,j,1)+zdzhs(1)*zlw_fr_kso(i,j,2)) &
                                  /zdzms(2)
               zdlw_fr_ksop05(i,j,1) = zdw(i,j)*EXP(zdw1(i,j)*                       &
                                   (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
               
               zinfmx = MIN(zinfmx, (zporv(i,j) - zw_fr(i,j,1))*zdzhs(1)*zrhwddt &
                      - rho_w*zdlw_fr_ksop05(i,j,1)*(zlw_fr_kso(i,j,2) - zlw_fr_kso(i,j,1)  )/zdzms(2))
               ! final infiltration rate limited by maximum value
               zinfil(i,j) = MIN(zinfmx,zrr(i,j))
               ! surface run-off (residual of potential minus actual infiltration)
               zro_inf       = zrr(i,j) - zinfil(i,j)
               runoff_s(i,j) = runoff_s(i,j) - zro_inf*zdt ! in mm
               raino(i,j) = raino(i,j)+zrr(i,j)*zdt ! in mm
            END IF            ! land-points only
         END DO
      END DO

!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------
      ! First loop through all layers from top to bottom, calculate diffusivity D
      ! and first-guess values for the conductivity K. These will then be limited
      ! in the second loop that goes from bottom to top.

      ! uppermost layer, kso = 1
      DO   j = jstarts, jends
         DO i = istarts, iends
    
            IF (llandmask(i,j)) THEN      ! land-points only
               IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type

                  zlw_fr_kso(i,j,1)  = zw_fr(i,j,1)
                  zlw_fr_ksop1= zw_fr(i,j,2)

        
                  ! interpolated scaled liquid water fraction at layer interface
                  zlw_fr_ksop05 = 0.5_ireals*(zdzhs(2)*zlw_fr_kso(i,j,1)+zdzhs(1)*zlw_fr_ksop1) &
                                  /zdzms(2)
                  zdlw_fr_ksop05(i,j,1) = zdw(i,j)*EXP(zdw1(i,j)*                       &
                                   (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
                  zklw_fr_ksop05_fg(i,j,1) = zkw(i,j,1)*EXP(zkw1(i,j)*                       &
                                   (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )!&
                                   !*max(0.0_ireals,min(1.0_ireals,1.0_ireals-(zlw_fr_ksop05-zfcap(i,j))/(zporv(i,j)-zfcap(i,j))))

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
                     limit=.false.
                     zlw_fr_ksom1        = zw_fr(i,j,kso-1)
                     zlw_fr_kso(i,j,kso) = zw_fr(i,j,kso  )
                     zlw_fr_ksop1        = zw_fr(i,j,kso+1)
                     ! interpolated scaled liquid water content at interface to layer
                     ! above and below

                     zlw_fr_ksom05 = 0.5*(zdzhs(kso-1)*zlw_fr_kso(i,j,kso)+   &
                                     zdzhs(kso)*zlw_fr_ksom1)/zdzms(kso)
                     zlw_fr_ksop05 = 0.5*(zdzhs(kso+1)*zlw_fr_kso(i,j,kso)+   &
                                     zdzhs(kso)*zlw_fr_ksop1)/zdzms(kso+1)
                     !zdlw_fr_ksop05_fg(i,j,kso-1) = zdw(i,j)*EXP( zdw1(i,j)*   &
                     !                (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )

                     if ((kso+1)==zsatlev(i,j)) then ! last level above saturation, partially saturated
                     !   !zdlw_fr_ksop05 = 0.0_ireals
                        ! first guess K
                        !zklw_fr_ksop05_fg(i,j,kso-1) = zkw(i,j,kso)*EXP( zkw1(i,j)*   &
                        !     (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
                        zklw_fr_ksop05_fg(i,j,kso) = zkw(i,j,kso)*EXP( zkw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )!&
                        !zklw_fr_ksop05_fg(i,j,kso) = 0.0_ireals ! No flux at the bottom
                        !q_kso_fc(i,j,kso) = 0.0
                        q_kso(i,j,kso) = s_topo*gamma*rho_w*zkw(i,j,kso)*h_wt ! saturated ground water flow at the bottom
                     else if (kso.ge.zsatlev(i,j)) then ! entirely saturated
                     !   zklw_fr_ksom05 = zkw(i,j)*EXP( -fsat(i,j)* (zmls(kso)-dc))
                     !   zklw_fr_ksop05 = zkw(i,j)*EXP( -fsat(i,j)* (zmls(kso)-dc))
                        !zklw_fr_ksop05_fg(i,j,kso-1) = 0.0_ireals
                        zklw_fr_ksop05_fg(i,j,kso) = zkw(i,j,kso)*EXP( zkw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )!&
                        !zklw_fr_ksop05_fg(i,j,kso)      = 0.0_ireals
                        !q_kso_fc(i,j,kso) = 0.0
                        q_kso(i,j,kso) = s_topo*gamma*rho_w*zkw(i,j,kso)*zdzhs(kso)
                     else ! default case, unsaturated
                        !zklw_fr_ksop05_fg(i,j,kso-1)= zkw(i,j,kso)*EXP( zkw1(i,j)*   &
                        !     (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )!&
                        zklw_fr_ksop05_fg(i,j,kso) = zkw(i,j,kso)*EXP( zkw1(i,j)*   &
                             (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )!&
                        zdlw_fr_ksop05(i,j,kso) = zdw(i,j)*EXP( zdw1(i,j)*   &
                                     (zporv(i,j)-zlw_fr_ksop05)/(zporv(i,j)-zadp(i,j)) )
                        ! runoff from this layer if field capacity is exceeded
                        hlp = zlw_fr_kso(i,j,kso) - zfcap(i,j)
                        q_kso(i,j,kso) = 0.0
                        !q_kso_fc_max = hlp*rho_w/(zdt*zdzhs(kso))
                        !q_kso_fc(i,j,kso) = zsf_heav(hlp)*min( q_kso_fc_max,s_topo*gamma*rho_w*zkw(i,j,kso)*h_wt_k(i,j,kso))
                        !if(q_kso_fc(i,j,kso).gt.0.)  print*,q_kso_fc(i,j,kso),q_kso_fc_max,kso, h_wt_k(i,j,kso)
                     end if
                     
                     !          IF(kso==ke_soil_hy-1) THEN
                     !            zflmg(i,j,kso+1)=rho_w*(zdlw_fr_ksop05*(zlw_fr_ksop1-zlw_fr_kso)/ &
                     !                            zdzms(kso+1)  - zklw_fr_ksop05)
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO

  DO   j = jstarts, jends
    DO i = istarts, iends
      IF (llandmask(i,j)) THEN      ! land-points only
        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
           limit=.false.
          ! lowest active hydrological layer ke_soil_hy-1
          zlw_fr_ksom1        = zw_fr(i,j,ke_soil_hy-1)
          zlw_fr_kso(i,j,ke_soil_hy) = zw_fr(i,j,ke_soil_hy  )
          zlw_fr_ksom05       = 0.5*(zdzhs(ke_soil_hy-1)*zlw_fr_kso(i,j,ke_soil_hy)+ &
                              zdzhs(ke_soil_hy)*zlw_fr_ksom1)/zdzms(ke_soil_hy)

          
          !zdlw_fr_ksop05_fg(i,j,ke_soil_hy-1)= zdw(i,j)*EXP( zdw1(i,j)* &
          !                  (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )

          if (ke_soil_hy.ge.zsatlev(i,j)) then ! fully saturated
             !zklw_fr_ksop05_fg(i,j,ke_soil_hy-1) = 0.0_ireals!zkw(i,j)*EXP( -fsat(i,j)* (zmls(ke_soil_hy)-dc))
             !q_ksop05(i,j,ke_soil_hy-1) = 0.0_ireals!zkw(i,j)*EXP( -fsat(i,j)* (zmls(ke_soil_hy)-dc))
             q_kso_fc(i,j,kso) = 0.0
             q_kso(i,j,ke_soil_hy) = s_topo*gamma*rho_w*zkw(i,j,ke_soil_hy)*zdzhs(ke_soil_hy)
          else ! partially saturated
             !zklw_fr_ksom05= zkw(i,j)*EXP( zkw1(i,j)* &
             !  (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )&
             !  *max(0.0_ireals,min(1.0_ireals,1.0_ireals+(zlw_fr_ksom05-zfcap(i,j))/(zporv(i,j)-zfcap(i,j))))
             !zklw_fr_ksop05_fg(i,j,ke_soil_hy-1) = zkw(i,j,ke_soil_hy)*EXP( zkw1(i,j)* &
             !                   (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
             q_kso_fc(i,j,kso) = 0.0
             q_kso(i,j,ke_soil_hy) = s_topo*gamma*rho_w*zkw(i,j,ke_soil_hy)*h_wt
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
        IF (llandmask(i,j)) THEN          ! land-points only
           IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type

              ! limit first guess fluxes such that no depletion can occur
              flux_ksom05 = (- zklw_fr_ksop05_fg(i,j,ke_soil_hy-1))/zdzhs(ke_soil_hy)
              flux_ksom05_pos = max(0.0_ireals,flux_ksom05)
              q_kso_pos       = max(0.0_ireals,q_kso(i,j,ke_soil_hy)/rho_w)
              !alpha_lim       = max(0.0_ireals,min(1.0_ireals,(flux_ksom05_pos*zdt)/(zlw_fr_kso(i,j,ke_soil_hy)+zepsi)))
              alpha_lim       = max(0.0_ireals,min(1.0_ireals,&
                                   zlw_fr_kso(i,j,ke_soil_hy)/(zdt*(flux_ksom05_pos+q_kso_pos+zepsi))))
              zklw_fr_ksop05_sg(i,j,ke_soil_hy-1) = alpha_lim * zklw_fr_ksop05_fg(i,j,ke_soil_hy-1)
              q_kso(i,j,ke_soil_hy)            = alpha_lim * q_kso(i,j,ke_soil_hy)
              
              if(alpha_lim.lt.1.0_ireals) then
                 print*,'limit alpha',ke_soil_hy, zklw_fr_ksop05_sg(i,j,ke_soil_hy-1),zklw_fr_ksop05_fg(i,j,ke_soil_hy-1), alpha_lim
              end if

              ! limit second guess K such that no oversaturation can occur in this layer
              zklw_fr_ksom05_limit =(zporv(i,j)-zlw_fr_kso(i,j,ke_soil_hy))/zdt*zdzhs(ke_soil_hy) &
!                                   - (q_kso(i,j,ke_soil_hy)+q_kso_fc(i,j,ke_soil_hy))*zdzhs(ke_soil_hy)/rho_w &
                                   - q_kso(i,j,ke_soil_hy)*zdzhs(ke_soil_hy)/rho_w

              zklw_fr_ksop05(i,j,ke_soil_hy-1) = min(zklw_fr_ksop05_sg(i,j,ke_soil_hy-1),zklw_fr_ksom05_limit)

              if(zklw_fr_ksop05(i,j,ke_soil_hy-1)==zklw_fr_ksom05_limit) then
                 limit=.true.
                 !print*,'limit'
              end if

              z1dgam1 = zdt/zdzhs(ke_soil_hy)
              zgam2m05  = zdlw_fr_ksop05(i,j,ke_soil_hy-1)/zdzms(ke_soil_hy)

              zaga(i,j,ke_soil_hy) = -zalfa* zgam2m05*z1dgam1
              zagb(i,j,ke_soil_hy) = 1.+ zalfa*zgam2m05*z1dgam1
              zagc(i,j,ke_soil_hy) = 0.0
              zagd(i,j,ke_soil_hy) = zw_fr(i,j,ke_soil_hy)+z1dgam1*zklw_fr_ksop05(i,j,ke_soil_hy-1) &
                                   + (1.-zalfa)*z1dgam1*                              &
                                     zgam2m05*(zlw_fr_kso(i,j,ke_soil_hy-1)  - zlw_fr_kso(i,j,ke_soil_hy))
              zflmg(i,j,ke_soil_hy)= 0.0
           ENDIF
        ENDIF
     END DO
  END DO

  ! inner loops
  DO kso = ke_soil_hy-1,2,-1
     DO   j = jstarts, jends
        DO i = istarts, iends
           IF (llandmask(i,j)) THEN          ! land-points only
              IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
                 ! limit first guess fluxes such that no depletion can occur
                 flux_ksom05 = (- zklw_fr_ksop05_fg(i,j,kso-1))/zdzhs(kso)
                 flux_ksop05 = (- zklw_fr_ksop05(i,j,kso))/zdzhs(kso)
                 flux_ksom05_pos = max(0.0_ireals,flux_ksom05)
                 flux_ksop05_neg = min(0.0_ireals,flux_ksop05)
                 !q_kso_pos       = max(0.0_ireals,(q_kso(i,j,kso)+q_kso_fc(i,j,kso))/rho_w)
                 q_kso_pos       = max(0.0_ireals,q_kso(i,j,kso)/rho_w)
                 !alpha_lim       = max(0.0_ireals,min(1.0_ireals,&
                 !                 (-flux_ksop05_neg+flux_ksom05_pos)*zdt/(zlw_fr_kso(i,j,kso)+zepsi)))
                 alpha_lim       = max(0.0_ireals,min(1.0_ireals,&
                                   zlw_fr_kso(i,j,kso)/(zdt*(-flux_ksop05_neg+flux_ksom05_pos+q_kso_pos+zepsi))))
                 zklw_fr_ksop05_sg(i,j,kso-1) = alpha_lim * zklw_fr_ksop05_fg(i,j,kso-1)
                 zklw_fr_ksop05(i,j,kso)      = alpha_lim * zklw_fr_ksop05(i,j,kso)
                 q_kso(i,j,kso)               = alpha_lim * q_kso(i,j,kso)
                 !q_kso_fc(i,j,kso)            = alpha_lim * q_kso_fc(i,j,kso)
                 if(alpha_lim.lt.1.0_ireals) then
                    print*,'limit alpha',kso, zklw_fr_ksop05_sg(i,j,kso-1),zklw_fr_ksop05_fg(i,j,kso-1), alpha_lim
                 end if
                 ! limit second guess K such that no oversaturation can occur in this layer
                 zklw_fr_ksom05_limit = (zporv(i,j)-zlw_fr_kso(i,j,kso))/zdt*zdzhs(kso) &
!                                      - (q_kso(i,j,kso)+q_kso_fc(i,j,kso))*zdzhs(kso)/rho_w &
                                      - q_kso(i,j,kso)*zdzhs(kso)/rho_w &
                                      + zklw_fr_ksop05(i,j,kso)
                 zklw_fr_ksop05(i,j,kso-1) = min(zklw_fr_ksop05_sg(i,j,kso-1),zklw_fr_ksom05_limit)
                 !zklw_fr_ksop05(i,j,kso-1) = zklw_fr_ksop05_sg(i,j,kso-1)
                 if(zklw_fr_ksop05(i,j,kso-1)==zklw_fr_ksom05_limit) then
                    limit=.true.
!                    print*,'limit',kso, zklw_fr_ksop05_sg(i,j,kso-1),zklw_fr_ksom05_limit, zklw_fr_ksop05(i,j,kso-1)
                 end if

                 ! coefficients for implicit flux computation
                 z1dgam1 = zdt/zdzhs(kso)
                 zgam2m05  = zdlw_fr_ksop05(i,j,kso-1)/zdzms(kso)
                 zgam2p05  = zdlw_fr_ksop05(i,j,kso)/zdzms(kso+1)

                 zaga (i,j,kso) = -zalfa*zgam2m05*z1dgam1
                 zagc (i,j,kso) = -zalfa*zgam2p05*z1dgam1
                 zagb (i,j,kso) = 1._ireals +zalfa*(zgam2m05+zgam2p05)*z1dgam1
                 zagd (i,j,kso) = zw_fr(i,j,kso)+                               &
                                  z1dgam1*(-zklw_fr_ksop05(i,j,kso)+zklw_fr_ksop05(i,j,kso-1))+ &
                                  (1._ireals-zalfa)*z1dgam1*                &
                                  (zgam2p05*(zlw_fr_kso(i,j,kso+1)-zlw_fr_kso(i,j,kso)  )     &
                                  -zgam2m05*(zlw_fr_kso(i,j,kso)  -zlw_fr_kso(i,j,kso-1)))
                 !soil water flux, explicit part, for soil water flux investigations
                 ! only)
                 zflmg(i,j,kso) = rho_w*(zdlw_fr_ksop05(i,j,kso-1)*(zlw_fr_kso(i,j,kso)-zlw_fr_kso(i,j,kso-1))/  &
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

              ! limit first-guess fluxes such that no depletion can occur
              flux_ksop05 = (- zklw_fr_ksop05(i,j,1))/zdzhs(1)
              flux_ksop05_neg = min(0.0_ireals,flux_ksop05)
              alpha_lim       = max(0.0_ireals,min(1.0_ireals,&
                                   zlw_fr_kso(i,j,1)/(zdt*(-flux_ksop05_neg+zepsi))))
              zklw_fr_ksop05(i,j,1)      = alpha_lim * zklw_fr_ksop05(i,j,1)
              zdlw_fr_ksop05(i,j,1)      = alpha_lim * zdlw_fr_ksop05(i,j,1)
              if(alpha_lim.lt.1.0_ireals) then
                 print*,'limit alpha uppermost layer', zklw_fr_ksop05_fg(i,j,1),zklw_fr_ksop05(i,j,1), alpha_lim
              end if
              ! coefficients for implicit flux computation
              z1dgam1 = zdt/zdzhs(1) 
              zgam2p05 = zdlw_fr_ksop05(i,j,1)/zdzms(2)
              zaga(i,j,1) = 0._ireals
              zagb(i,j,1) = 1._ireals+zalfa*zgam2p05*z1dgam1
              zagc(i,j,1) = -zalfa * zgam2p05*z1dgam1
              zagd(i,j,1) = zw_fr(i,j,1) + zinfil(i,j)*z1dgam1/rho_w  &
                           -zklw_fr_ksop05(i,j,1)*z1dgam1                     &
                           +(1. - zalfa)* zgam2p05*z1dgam1*(zlw_fr_kso(i,j,2) - zlw_fr_kso(i,j,1))
           ENDIF
        ENDIF
     END DO
  END DO
  ! diagnose water-table depth and ground runoff
  !runoff_g(:,:) = 0.0
  DO   kso =1,ke_soil_hy
     DO   j = jstarts, jends
        DO i = istarts, iends
           !if(kso==zsatlev(i,j)-1) then! last partly saturated level above ground water
           !   runoff_g(i,j) = runoff_g(i,j)- q_kso(i,j,kso)*h_wt*zdt/rho_w
           !else if(kso.ge.zsatlev(i,j)) then
              runoff_g(i,j) = runoff_g(i,j)- q_kso(i,j,kso)*zdzhs(kso)*zdt ! in mm
              !runoff_g(i,j) = runoff_g(i,j)- q_kso_fc(i,j,kso)*zdzhs(kso)*zdt ! in mm
           !end if
        END DO
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


!  DO kso = ke_soil_hy,1,-1
!     DO   j = jstarts, jends
!        DO i = istarts, iends
!            dwdto(i,j,kso) =  dwdto(i,j,kso) + zdt*(w_so(i,j,kso,nnew)/zdzhs(kso)-zw_fr(i,j,kso))*zdzhs(kso)
!           zlw_fr_kso_new  = w_so(i,j,kso  ,nnew)/zdzhs(kso  )
!           IF (zlw_fr_kso_new>zporv(i,j)) THEN
!              w_so(i,j,kso-1,nnew) = w_so(i,j,kso-1,nnew)+(zlw_fr_kso_new-zporv(i,j))*zdzhs(kso-1)
!              w_so(i,j,kso  ,nnew) = zporv(i,j)*zdzhs(kso)
!           END IF
!        END DO
!     END DO
!  END DO
! combine implicit part of sedimentation and capillary flux with explicit part
! (for soil water flux investigations only)
!DO kso = 2,ke_soil+1
!  DO   j = jstarts, jends
!    DO i = istarts, iends
!      IF (llandmask(i,j)) THEN          ! land-points only
!        IF (m_styp(i,j) >= 3) THEN   ! neither ice nor rock as soil type
!          zlw_fr_ksom1_new= w_so(i,j,kso-1,nnew)/zdzhs(kso-1)
!          zlw_fr_kso_new  = w_so(i,j,kso  ,nnew)/zdzhs(kso  )
!          zlw_fr_ksom1  = w_so(i,j,kso-1,nx)/zdzhs(kso-1)
!          zlw_fr_kso    = w_so(i,j,kso  ,nx)/zdzhs(kso  )
!          !... additionally for runoff_g at lower level of lowest active water
!          ! layer calculated with (upstream) main level soil water content
!          ! compute reduction factor for transport coefficients
          
!        ! interpolated liquid water content at interface to layer above
!          zlw_fr_ksom05 =0.5*(zdzhs(kso)*zlw_fr_ksom1+zdzhs(kso-1)*zlw_fr_kso) &
!                            /zdzms(kso)
!          zdlw_fr_ksom05= zdw(i,j)*EXP(zdw1(i,j) *  &
!                          (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )
!          zklw_fr_ksom05= zkw(i,j) * EXP(zkw1(i,j)* &
!                          (zporv(i,j)-zlw_fr_ksom05)/(zporv(i,j)-zadp(i,j)) )

!          IF (kso> ke_soil_hy) zdlw_fr_ksom05=0.0   ! no flux gradient 
                                                    ! contribution below 2.5m
!          IF (kso> ke_soil_hy) zklw_fr_ksom05=0.0   ! no gravitation flux below 2.5m
          
!          zflmg(i,j,kso) =                              &
!                   (1._ireals-zalfa) * zflmg(i,j,kso) + & ! explicit flux component
!                              zalfa  * rho_w *          & ! implicit flux component
!                   (zdlw_fr_ksom05 * (zlw_fr_kso_new-zlw_fr_ksom1_new)/zdzms(kso) &
!                    - zklw_fr_ksom05)

!          zklw_fr_kso_new = zkw(i,j) * EXP(zkw1(i,j)* &
!                            (zporv(i,j) - zlw_fr_kso_new)/(zporv(i,j) - zadp(i,j)) )
!          ! actual gravitation water flux
!          IF(w_so (i,j,kso,nnew).LT.1.01_ireals*zadp(i,j)*zdzhs(kso)) zklw_fr_kso_new = 0._ireals
         
!          zrunoff_grav(i,j,kso) =  0.0!- rho_w * zklw_fr_kso_new

!          ! ground water as lower boundary of soil column
!          IF ((kso == ke_soil_hy+1).and.(itype_hydbound == 3)) THEN
!             zdelta_sm=( zlw_fr_kso_new - zlw_fr_ksom1_new )

!             zdlw_fr_kso = zdw(i,j)*EXP(zdw1(i,j) *  &
!                  (zporv(i,j)-zlw_fr_kso_new)/(zporv(i,j)-zadp(i,j)) )
!             zklw_fr_kso = zkw(i,j) * EXP(zkw1(i,j)* &
!                  (zporv(i,j)-zlw_fr_kso_new)/(zporv(i,j)-zadp(i,j)) )
!             zklw_fr_ksom1 = zkw(i,j) * EXP(zkw1(i,j)* &
!                  (zporv(i,j)-zlw_fr_ksom1_new)/(zporv(i,j)-zadp(i,j)) )

!             zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
!             zrunoff_grav(i,j,ke_soil_hy)=zrunoff_grav(i,j,ke_soil_hy)+ &
!                  zdhydcond_dlwfr / &
!                  (1.0-exp(-zdhydcond_dlwfr/zdlw_fr_kso*0.5*zdzms(ke_soil_hy+1)))* &
!                  zdelta_sm
!          ENDIF
!        END IF 
!      END IF          ! land-points only
!    END DO
!  END DO
!END DO

  DO  kso = 1,ke_soil
    ! utility variables used to avoid if-constructs in following loops
    zro_sfak = zsf_heav(0.5_ireals + msr_off - kso)  ! 1.0 for 'surface runoff'
    zro_gfak = 1._ireals - zro_sfak                  ! 1.0 for 'ground runoff'

    IF (kso==ibot_w_so) THEN
      !zfmb_fak = 1.0_ireals
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
!            zdwg =  zdwg + zrunoff_grav(i,j,kso)*zfmb_fak
            ! runoff from soil layer if field capacity is exceeded and negative
            ! flux divergence
            !zredfu =  MAX( 0.0_ireals, MIN( 1.0_ireals,(zw_fr(i,j,kso) -     &
            !           zfcap(i,j))/MAX(zporv(i,j) - zfcap(i,j),zepsi)) )
            !zredfu =  MAX( 0.0_ireals, MIN( 1.0_ireals,(zw_fr(i,j,kso) -     &
            !           zporv(i,j))/MAX(zporv(i,j) - zporv(i,j),zepsi)) )
            !zredfu = zsf_heav(zdwg)*zredfu
            !zro    = -zdwg*zredfu
            !zdwg   = zdwg*(1._ireals - zredfu)

            if (kso.le.irootdp) then
               evap_max   = w_so(i,j,kso,nnew)/zdt*irootdp
               evap(nout) = min(evap(nout),evap_max)
               evapo(i,j) = evapo(i,j)-evap(nout)/irootdp*zdt*rho_w
!                 dwdto(i,j,kso)     =  dwdto(i,j,kso) - zdtdrhw*(evap(nout)/irootdp)*zdzhs(kso)&
!                      -zdtdrhw*q_kso(i,j,kso)
               zdwg = zdwg - rho_w*(evap(nout)/irootdp)
            end if

            zwgn            = zw_fr(i,j,kso) + zdtdrhw*zdwg/zdzhs(kso)
            zro2            = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zwgn - zporv(i,j))
            zkorr           = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zadp(i,j) - zwgn )
            ! first contribution to ground runoff: above field capacity
            zdwgdt(i,j,kso) = zdwg !+ zkorr - zro2
            zro             = 0.0!- zro2
            !runoff_s(i,j) = runoff_s(i,j) + zro*zro_sfak*zroffdt
            !runoff_g(i,j) = runoff_g(i,j) + zro*zro_gfak*zroffdt
!            ! runoff_g reformulation:
!            runoff_g(i,j) = runoff_g(i,j) - (zrunoff_grav(i,j,kso) * zfmb_fak &
!                                          + zkorr) * zroffdt

            w_so(i,j,kso,nnew) = zw_fr(i,j,kso)*zdzhs(kso) +zdtdrhw*zdwgdt(i,j,kso)
          END IF          ! ice/rock-exclusion
        END IF   ! land-points only
      END DO
    END DO
  END DO         ! end loop over soil layers

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

! remove some soil water caused by evaporation

  DO kso = 1,ke_soil
     DO   j = jstarts, jends
        DO i = istarts, iends
           IF (llandmask(i,j)) THEN  ! for landpoints only
              w_so(i,j,kso,nnew) = w_so(i,j,kso,nnew) -zdtdrhw*q_kso(i,j,kso)*zdzhs(kso)
              !w_so(i,j,kso,nnew) = w_so(i,j,kso,nnew) -zdtdrhw*q_kso_fc(i,j,kso)*zdzhs(kso)
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
          satlev_out(i,j,nout) = zsatlev(i,j)
          rsatlev_out(i,j,nout)= z_wt
          DO kso=1,ke_soil+1
             w_so_out(i,j,kso,nout) = w_so(i,j,kso,nnew)
             wl_so_out(i,j,kso,nout) = zlw_fr_kso(i,j,kso)
             s_so_out(i,j,kso,nout) = w_so(i,j,kso,nnew)/(zporv(i,j)*zdzhs(kso))
             dwdt_out(i,j,kso,nout) = dwdto(i,j,kso)
             dwdto(i,j,kso) = 0.0_ireals
             klwm_out(i,j,kso,nout)  = zklw_fr_ksop05_fg(i,j,kso)
             klwp_out(i,j,kso,nout)  = zklw_fr_ksop05(i,j,kso)
             dlwpf_out(i,j,kso,nout) = zdlw_fr_ksop05_fg(i,j,kso)
             dlwp_out(i,j,kso,nout)  = zdlw_fr_ksop05(i,j,kso)
             q_kso_out(i,j,kso,nout) = q_kso(i,j,kso)
             q_kso_fc_out(i,j,kso,nout) = q_kso_fc(i,j,kso)
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

 iret = nf90_def_var(ncid, 'WL_SO', NF90_REAL, (/zdim, timedim/), wlsoid)
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

 iret = nf90_def_var(ncid, 'SATLEV', NF90_INT, (/timedim/), slid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'RSATLEV', NF90_REAL, (/timedim/), rslid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'KP05_FG', NF90_REAL, (/zhdim, timedim/), kmid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'KP05_LIM', NF90_REAL, (/zhdim, timedim/), kpid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'DP05_FG', NF90_REAL, (/zhdim, timedim/), dpfid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'DP05_LIM', NF90_REAL, (/zhdim, timedim/), dpid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'Q', NF90_REAL, (/zhdim, timedim/), qid)
 call check_err(iret)

 iret = nf90_def_var(ncid, 'QFC', NF90_REAL, (/zhdim, timedim/), qfcid)
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

 iret = nf90_put_att(ncid, slid,'standard_name',"first saturated level")
 call check_err(iret)
 iret = nf90_put_att(ncid, slid,'long_name',"first saturated level")
 call check_err(iret)
 iret = nf90_put_att(ncid, slid,'units',"")
 call check_err(iret)

 iret = nf90_put_att(ncid, rslid,'standard_name',"depth of water table")
 call check_err(iret)
 iret = nf90_put_att(ncid, rslid,'long_name',"depth of water table")
 call check_err(iret)
 iret = nf90_put_att(ncid, rslid,'units',"m")
 call check_err(iret)

 iret = nf90_put_att(ncid, kmid,'standard_name',"first guess conductivity")
 call check_err(iret)
 iret = nf90_put_att(ncid, kmid,'long_name',"first guess conductivity at layer +1/2")
 call check_err(iret)
 iret = nf90_put_att(ncid, kmid,'units',"m/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, kpid,'standard_name',"limited conductivity")
 call check_err(iret)
 iret = nf90_put_att(ncid, kpid,'long_name',"limited conductivity at layer +1/2")
 call check_err(iret)
 iret = nf90_put_att(ncid, kpid,'units',"m/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, dpfid,'standard_name',"first guess diffusivity")
 call check_err(iret)
 iret = nf90_put_att(ncid, dpfid,'long_name',"first guess diffusivity at layer +1/2")
 call check_err(iret)
 iret = nf90_put_att(ncid, dpfid,'units',"m/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, dpid,'standard_name',"limited diffusivity")
 call check_err(iret)
 iret = nf90_put_att(ncid, dpid,'long_name',"limited diffusivity at layer +1/2")
 call check_err(iret)
 iret = nf90_put_att(ncid, dpid,'units',"m/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, rgid,'standard_name',"ground water runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rgid,'long_name',"ground water runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rgid,'units',"mm/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, rsid,'standard_name',"surface runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rsid,'long_name',"surface runoff")
 call check_err(iret)
 iret = nf90_put_att(ncid, rsid,'units',"mm/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, qid,'standard_name',"ground water flow, aquifer")
 call check_err(iret)
 iret = nf90_put_att(ncid, qid,'long_name',"ground water flow, aquifer")
 call check_err(iret)
 iret = nf90_put_att(ncid, qid,'units',"mm/s")
 call check_err(iret)

 iret = nf90_put_att(ncid, qfcid,'standard_name',"ground water flow, field capacity exceed")
 call check_err(iret)
 iret = nf90_put_att(ncid, qfcid,'long_name',"ground water flow, field capacity exceed")
 call check_err(iret)
 iret = nf90_put_att(ncid, qfcid,'units',"mm/s")
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

 IRET = NF90_PUT_VAR(ncid, slid, satlev_out(1,1,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, rslid, rsatlev_out(1,1,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, kmid, klwm_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, kpid, klwp_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, dpfid, dlwpf_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, dpid, dlwp_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, qid, q_kso_out(1,1,:,:))
 IF (IRET .NE. NF90_NOERR) PRINT *, NF90_STRERROR(IRET)

 IRET = NF90_PUT_VAR(ncid, qfcid, q_kso_fc_out(1,1,:,:))
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
