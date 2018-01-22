!+ Data module for control parameters, some constants, some input
!------------------------------------------------------------------------------

MODULE data_namelist
!==============================================================================
!
! Declarations:
!
! Modules used:
USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables


!==============================================================================

IMPLICIT NONE

!==============================================================================
INTEGER (KIND=iintegers), PARAMETER ::  &
! horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie = 1_iintegers,           & ! number of grid points in zonal direction
    je = 1_iintegers,           & ! number of grid points in meridional direction
    ke_soil = 9_iintegers,      & ! number of layers in multi-layer soil model
    ie_tot = 1_iintegers,       & ! number of grid points in zonal direction
    je_tot = 1_iintegers,       & ! number of grid points in meridional direction

    istartpar = 1_iintegers,    & ! start index for computations in the parallel program
    iendpar = 1_iintegers,      & ! end index for computations in the parallel program
    jstartpar = 1_iintegers,    & ! start index for computations in the parallel program
    jendpar = 1_iintegers        ! end index for computations in the parallel program
REAL  (KIND=ireals), PARAMETER ::  &
! variables for the time discretization and related variables
! --------------------------------------------------------------
    dt = 20.0_ireals            ,  &  !  time-step
    starttime = 0.0_ireals      ,  & ! start time
    endtime = 80._ireals*604800_ireals   , & ! end time 604800 seconds = 7 days
    outtime = 2.*86400.0_ireals       ! after how much time output is done
INTEGER (KIND=iintegers), PARAMETER ::  &
! integer variables for time discretization
   nsteps =  INT((endtime - starttime)/dt) , & ! number of time steps
   noutput = INT((endtime - starttime)/outtime) ! number of output steps
REAL  (KIND=ireals), PARAMETER ::  &
! constants for parametrizations
! ---------------------------------
    rho_w = 1000.0_ireals      ,    &   ! density of liquid water
    plcov(1,1) = 0.0_ireals             ! fraction of plant cover                         --
INTEGER (KIND=iintegers), PARAMETER ::  &
     itype_hydbound = 1 ! lower boundary
!------------------------------------------------------------------------------
! external parameter fields                                        (unit)
! ----------------------------
REAL (KIND=ireals), PARAMETER ::  &
    soiltyp(1,1) = 5_ireals             ! type of the soil (keys 0-9)
LOGICAL, PARAMETER ::           &
    llandmask(1,1) = .TRUE.             ! landpoint mask

LOGICAL, PARAMETER ::           &
    ldecharme = .TRUE.           ! decharme formulation

LOGICAL, PARAMETER ::           &
    lexpporv  = .TRUE.           ! exponential pore volume profile

REAL (KIND=ireals) ::  &
    rootdp(1,1)             ! root depth

REAL (KIND=ireals   )     s_topo, gamma

REAL  (KIND=ireals) ::  &                         !    --
! fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    w_so_init(1,1,ke_soil+1)      , &  ! total water conent (liquid water)       (m H20)
    prec(noutput+1)    , & ! precipitation for each output timestep
    evap(noutput+1)    , & ! evaporation for each output timestep
    czmls(ke_soil+1)           , &  ! depth of the main level soil layers in m
    czbot_w_so =4.0_ireals,   & ! depth of last hydrological active soil layer in m
    czrootdp = 0.5_ireals      ! depth of root zone in m

! fields for model output and diagnostics                          (unit )
! ---------------------------------------------------------------

REAL  (KIND=ireals) ::  &
    ignore            ! dummy variable for input reading

INTEGER (KIND=iintegers) ::  &
    ind               ! dummy variable for loop index in namelist
!-----------------------------------------------------------------------------
!=======================================================================

CONTAINS

!==============================================================================
!==============================================================================
!+ Subroutine for initialization
!------------------------------------------------------------------------------

SUBROUTINE init_fields
!=============================================================================
! INIT
!=============================================================================

s_topo = 0.05_ireals
gamma  = 0.25_ireals


OPEN(UNIT=13, FILE="w_so_out", ACTION="read", STATUS="old")
read(13,'(F12.0)',advance='no') ignore ! ignore first entry as it contains time
DO ind=1,ke_soil-1,1 ! loop over all layers except last two
       read(13,'(F12.8)', advance='no')  czmls(ind) ! read in all depths except the last
END DO
       read(13,'(F11.8)',advance='no') czmls(ke_soil) ! second last depth
       read(13,'(F13.8)', advance='no') czmls(ke_soil+1) ! last depth
read(13,'(F12.0)') ignore ! ignore first entry as it contains time

read(13,'(F12.0)', advance='no') ignore ! ignore first entry as it contains time
! Read in water content
DO ind=1,ke_soil+1,1 ! loop over all layers
       read(13,'(F10.8, F2.0)', advance='no')  w_so_init(1,1,ind), ignore !
END DO
 CLOSE(13)
!DO ind=1,ke_soil+1,1 ! loop over all layers
!   czmls(ind) = 0.001*2.01**(real(ind))
!       print*,czmls(ind)
!END DO
! Set the precipitation to a constant value or a constant value for part of the runtime e.g. here for half of the total time
        DO ind=1,noutput/2,1
                  prec(ind) = 5.0_ireals/400._ireals*dt*0.0046296_ireals ! Const precipitation
!                  prec(ind) = 20.*0.00463_ireals ! Const precipitation
        END DO
       DO ind=noutput/2+1,noutput+1
                  prec(ind) = 0.0_ireals ! No precipitation
       END DO
!        DO ind=1,noutput
!           prec(ind) = 1./20.*0.00463_ireals ! Const precipitation
!        END DO

! evaporation
        DO ind=1,noutput/2,1
                  evap(ind) = 0.0_ireals ! No evaporation
        END DO
       DO ind=noutput/2+1,noutput+1
                  evap(ind) = 0.5_ireals/400._ireals*1._ireals/rho_w*0.00463_ireals ! Const evaporation
       END DO
!        DO ind=1,noutput
!           prec(ind) = 1./20.*0.00463_ireals ! Const precipitation
!        END DO
!       DO ind=1,noutput+1
!!                  evap(ind) = 1./1000.*0.00463_ireals ! Const evaporation
!                  evap(ind) = 0.00_ireals ! Const evaporation
!       END DO
END SUBROUTINE init_fields

END MODULE data_namelist
