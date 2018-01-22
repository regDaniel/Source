!+ Data module for all parametric data in the soil model "terra"  
!------------------------------------------------------------------------------

MODULE data_soil

!------------------------------------------------------------------------------
!
! Description:
!  This module declares and initializes all parametric scalar and array      
!  data which are used in the soil model (terra1 and terra2)
!
! Current Code Owner: DWD, Juergen Helmert
!  phone:  +49  69  8062 2704
!  fax:    +49  69  8062 3721
!  email:  Juergen.Helmert@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenther Doms
!  Initial release
! 1.30       1999/06/24 Erdmann Heise
! Implementation of variables for simplified BATS-scheme
! 1.33       1999/10/14 Matthias Raschendorfer
!  crsmin is now a namelist-parameter.
! 2.17       2002/05/08 Ulrich Schaettler
!  Additional parameters for soil water content dependent freezing/melting
! 2.18       2002/07/16 Reinhold Schrodin
!  Redefined variable cf_snow (for calculation of fractional snow coverage)
! 3.6        2003/12/11 Reinhold Schrodin
!  Adapted several variables for new multi-layer soil model
! 3.13       2004/12/03 Reinhold Schrodin
!  New variable cwimax_ml (maximum interception water content for multi-layer
!  soil model). Changed values for minimal and maximal density of snow
!  (crhosmin_ml: 100 => 250; crhosmax_ml: 400 => 250) (now consistent with GME)
! 3.17       2005/12/12 Reinhold Schrodin
!  New variables (crhosmin, crhosmaxf, crhosmin, crhosmaxt, csnow_tmin)
!  and changed variables (crhosmin_ml,crhosmax_ml) for ageing of snow density
!  calculation
! 3.18       2006/03/03 Ulrich Schaettler
!  Editorial changes
! 3.21       2006/12/04 Ulrich Schaettler
!  Put declaration of NL parameters crsmin and rat_lam to data_soil
! V3_23        2007/03/30 Matthias Raschendorfer
!  Moving 'rat_lam' into MODULE 'data_turbulence'.
!  Initialisation of 'crsmin' with a new default value of 150.0
! V4_11        2009/11/30 Ekaterina Machulskaya
!  Introduced parameters for the snow model
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced new logical lsoilinit_dfi to initialize soil variables after
!    a DFI forward launching
!  Changed the code owner
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
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. Data arrays for properties of different soil types (array index)     
! -------------------------------------------------------------------
 
  REAL  (KIND=ireals) ::  &
!   a) parameters describing the soil water budget
    cporv (10), &  !  pore volume (fraction of volume)
    cfcap (10), &  !  field capacity (fraction of volume)
    cpwp  (10), &  !  plant wilting point (fraction of volume)
    cadp  (10), &  !  air dryness point (fraction of volume)
    cik2  (10), &  !  minimum infiltration rate (kg/s*m**2)
    ckw0  (10), &  !  parameter for determination of hydr. conductivity (m/s)
    ckw1  (10), &  !  parameter for determination of hydr. conductivity (1)
    cdw0  (10), &  !  parameter for determination of hydr. diffusivity (m**2/s)
    cdw1  (10), &  !  parameter for determination of hydr. diffusivity (1)
    crock (10), &  !  rock/ice/water indicator (hydrological calculations 
                   !  only for crock=1)

!   b) parameters describing the soil heat budget
    cdz1  (10), &  !  top layer thickness (EFR-method)
    crhoc (10), &  !  soil heat capacity  (J/K*m**3)
    cala0 (10), &  !  parameters for the determination of
    cala1 (10), &  !      the soil heat conductivity (W/(K*m))
    csalb (10), &  !  solar albedo for dry soil                            
    csalbw(10), &  !  slope of solar albedo with respect to soil water content     

!   c) additional parameters for the BATS scheme (Dickinson)
    ck0di (10), &  !  (m/s)
    cbedi (10), &  !  (1)
    clgk0 (10), &  !  auxiliary variable

!   d) additional parameters for soil water content dependent freezing/melting
    csandf(10), &  !  mean fraction of sand (weight percent)
    cclayf(10)     !  mean fraction of clay (weight percent)
 

  ! Initialization of soil type parameters except cdz1 
  ! (being calculated during execution)
    
  ! soil type:   ice    rock    sand    sandy   loam   clay      clay    peat    sea     sea  
  ! (by index)                          loam           loam                     water    ice
      
  DATA  cporv / 1.E-10, 1.E-10,  .364 ,  .445 ,  .455 ,  .475 ,  .507 ,  .863 , 1.E-10, 1.E-10 /
  DATA  cfcap / 1.E-10, 1.E-10,  .196 ,  .260 ,  .340 ,  .370 ,  .463 ,  .763 , 1.E-10, 1.E-10 /
  DATA  cpwp  / 0.0   , 0.0   ,  .042 ,  .100 ,  .110 ,  .185 ,  .257 ,  .265 , 0.0   ,  0.0   /
  DATA  cadp  / 0.0   , 0.0   ,  .012 ,  .030 ,  .035 ,  .060 ,  .065 ,  .098 , 0.0   ,  0.0   /
  DATA  crhoc /1.92E6 , 2.10E6, 1.28E6, 1.35E6, 1.42E6, 1.50E6, 1.63E6, 0.58E6, 4.18E6, 1.92E6 /
  DATA  cik2  / 0.0   , 0.0   , 0.0035, 0.0023, 0.0010, 0.0006, 0.0001, 0.0002, 0.0   ,  0.0   /
  DATA  ckw0  / 0.0   , 0.0   , 479E-7, 943E-8, 531E-8, 764E-9,  17E-9,  58E-9, 0.0   ,  0.0   /
  DATA  ckw1  / 0.0   , 0.0   , -19.27, -20.86, -19.66, -18.52, -16.32, -16.48, 0.0   ,  0.0   /
  DATA  cdw0  / 0.0   , 0.0   , 184E-7, 346E-8, 357E-8, 118E-8, 442E-9, 106E-9, 0.0   ,  0.0   /
  DATA  cdw1  / 0.0   , 0.0   , -8.45 ,  -9.47, -7.44 , -7.76 , -6.74 , -5.97 , 0.0   ,  0.0   /
  DATA  crock / 0.0   , 0.0   ,  1.0  ,   1.0 ,  1.0  ,  1.0  ,  1.0  ,  1.0  , 0.0   ,  0.0   /
  DATA  cala0 / 2.26  , 2.41  , 0.30  ,  0.28 ,  0.25 ,  0.21 ,  0.18 ,  0.06 , 1.0   ,  2.26  /
  DATA  cala1 / 2.26  , 2.41  , 2.40  ,  2.40 ,  1.58 ,  1.55 ,  1.50 ,  0.50 , 1.0   ,  2.26  /
  DATA  csalb / 0.70  , 0.30  ,  0.30 ,  0.25 , 0.25  , 0.25  , 0.25  , 0.20  , 0.07  ,  0.70  /
  DATA  csalbw/ 0.00  , 0.00  ,  0.44 ,  0.27 , 0.24  , 0.23  , 0.22  , 0.10  , 0.00  ,  0.00  /
  DATA  ck0di / 1.E-4 , 1.E-4 , 2.E-4 , 2.E-5 , 6.E-6 , 2.E-6 , 1.E-6 , 1.5E-6, 0.00  ,  0.00  /
  DATA  cbedi / 1.00  , 1.00  ,  3.5  ,  4.8  , 6.1   , 8.6   , 10.0  , 9.0   , 0.00  ,  0.00  /
  DATA  csandf/ 0.0   , 0.0   ,  90.  ,  65.  , 40.   , 35.   , 15.  ,  90.   , 0.00  ,  0.00 /
  DATA  cclayf/ 0.0   , 0.0   ,  5.0  ,  10.  , 20.  ,  35.  ,  70.  ,   5.0  , 0.00  ,  0.00 /
 

!==============================================================================

! 2. Additional parameters for the soil model                             
! -------------------------------------------------------------------

  REAL  (KIND=ireals) ::  &
!==============================================================================

    csalb_p    = 0.15_ireals  , & !  solar albedo of ground covered by plants
    csalb_snow = 0.70_ireals  , & !  solar albedo of ground covered by snow
    csalb_snow_min = 0.400_ireals, &
                           ! min. solar albedo of snow for forest free surfaces
    csalb_snow_max = 0.700_ireals, &
                           ! max. solar albedo of snow for forest free surfaces
  ! for possible later use:
    csalb_snow_fe  = 0.200_ireals , &  ! solar albedo of snow for surfaces with evergreen forest
    csalb_snow_fd  = 0.200_ireals , &  ! solar albedo of snow for surfaces with deciduous forest
    ctalb      = 0.004_ireals , & !  thermal albedo ( of all soil types )   
    cf_snow    = 0.0150_ireals, & !  parameter for the calculation of the 
                                  !  fractional snow coverage
  ! for the multi-layer soil model
    cwhc       = 0.04_ireals,   & !  water holding capacity of snow ()
    chcond     = 0.01_ireals,   & !  saturation hydraulic conductivity of snow ()
    ca2        = 6.6E-07_ireals,& !  activation energy (for snow metamorphosis) (J)
    csigma     = 75._ireals,    & !  snow metamorphosis, Pa

  ! cf_w changed from 0.0004 to 0.0010 (in agreement with GME)
    cf_w       = 0.0010_ireals, & !  parameter for the calculation of the
                                  !  fractional water coverage

    csvoro     = 1.0000_ireals, & !  parameter to estimate the subgrid-scale 
                                  !  variation of orography
    cik1       = 0.0020_ireals, & !  parameter for the determination of the 
                                  !  maximum infiltaration
    cwimax     = 0.0005_ireals, & !  parameter for the determination of the 
    cwimax_ml  = 1.E-6_ireals,  & !  maximum interception water content
    ctau_i     = 1000.0_ireals, & !  time constatant for the drainage from the 
                                  !  interception storeage 
    cakw       = 0.8000_ireals, & !  parameter for averaging the water contents
                                  !  of the top and middle soil water layers to 
                                  !  calculate the hydraulic diffusivity and 
                                  !  conductiviy

    ctau1      = 1.0000_ireals, & !  first adjustment time period in EFR-method
    ctau2      = 5.0000_ireals, & !  second adjustment time period in EFR-method
    chc_i      = 2100.0_ireals, & !  heat capacity of ice     
    chc_w      = 4180.0_ireals, & !  heat capacity of water     

    cdzw12     = 0.1000_ireals, & !  thickness of upper soil water layer in 
                                  !  two-layer model         
    cdzw22     = 0.9000_ireals, & !  thickness of lower soil water layer in 
                                  !  two-layer model      
    cdzw13     = 0.0200_ireals, & !  thickness of upper soil water layer in 
                                  !  three-layer model
    cdzw23     = 0.0800_ireals, & !  thickness of middle soil water layer in 
                                  !  three-layer model 
    cdzw33     = 0.9000_ireals    !  thickness of lower soil water layer in 
                                  !  three-layer model

  REAL  (KIND=ireals) ::  &
    cdsmin     = 0.0100_ireals, & !  minimum snow depth
    crhosmin   = 500.00_ireals, & !  minimum density of snow
    crhosmax   = 800.00_ireals, & !  maximum density of snow
    crhosmin_ml=  50.00_ireals, & !  minimum density of snow
    crhosmax_ml= 400.00_ireals, & !  maximum density of snow
    crhosminf  =  50.00_ireals, & !  minimum density of fresh snow
    crhosmaxf  = 150.00_ireals, & !  maximum density of fresh snow
    crhosmint  =   0.20_ireals, & !  minimum value of time constant for ageing 
                                  !  of snow
    crhosmaxt  =   0.40_ireals, & !  maximum value of time constant for ageing 
                                  !  of snow
    csnow_tmin = 258.15_ireals, & !  lower threshold temperature of snow for 
                                  !  ageing and fresh snow density computation 
                                  !  ( = 273.15-15.0)
    crhos_dw   = 300.00_ireals, & !  change of snow density with water content
    calasmin   = 0.2000_ireals, & !  minimum heat conductivity of snow (W/m K)
    calasmax   = 1.5000_ireals, & !  maximum heat conductivity of snow (W/m K)
    calas_dw   = 1.3000_ireals, & !  change of snow heat conductivity with
                                  !  water content                (W/(m**2) K)
   
    crhowm     =    0.8_ireals    , & !  BATS (1)
    cdmin      =    0.25E-9_ireals, & !  BATS (m**2/s)
    cfinull    =    0.2_ireals    , & !  BATS (m)
    ckrdi      =    1.0E-5_ireals , & !  BATS (m/s)
    cdash      =    0.05_ireals   , & !  BATS ((m/s)**1/2)
    clai       =    3.0_ireals    , & !  BATS
    cparcrit   =  100.0_ireals    , & !  BATS (W/m**2)
    ctend      =  313.15_ireals   , & !  BATS (K)
    csatdef    = 4000.0_ireals    , & !  BATS (Pa)

    !Minimum and maximum value of stomatal resistance (s/m)
    !used by the Pen.-Mont. method for vegetation transpiration
    !(itype_trvg=2):
    crsmin     = 150.0_ireals     , & !  BATS (s/m)
    crsmax     = 4000.0_ireals        !  BATS (s/m)

! crsmax increased from 1000 to 4000 s/m (to reduce latent heat flux).

! 3. Additional control variables
! -------------------------------

  LOGICAL                     ::  &
    lsoilinit_dfi = .FALSE.         ! initialize soil after dfi forward launching

!==============================================================================

END MODULE data_soil     
