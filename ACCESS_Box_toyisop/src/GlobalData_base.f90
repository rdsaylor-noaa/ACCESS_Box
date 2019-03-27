!======================================================================================================================!
!                                                                                                                      !
!     Program:      ACCESS                                                                                             !
!                   Atmospheric Chemistry and Canopy Exchange Simulation System                                        !
!                   Box model version                                                                                  !
!                                                                                                                      !
!     Version:      3.1.0                                                                                              !
!                                                                                                                      !
!======================================================================================================================!
!                                                                                                                      !
!     Last Update:  Mar 2019                                                                                           !
!                                                                                                                      !
!     Contact:      Rick D. Saylor, PhD                                                                                !
!                   Physical Scientist                                                                                 !
!                   U. S. Department of Commerce                                                                       !
!                   National Oceanic and Atmospheric Administration                                                    !
!                   Air Resources Laboratory                                                                           !
!                   Atmospheric Turbulence and Diffusion Division                                                      !
!                   456 S. Illinois Ave                                                                                !
!                   Oak Ridge, TN 37830                                                                                !
!                   email: Rick.Saylor@noaa.gov                                                                        !
!                                                                                                                      !
!**********************************************************************************************************************!
!                   NOTE: See Legal Notice in Main.f90                                                                 !
!**********************************************************************************************************************!
!                                                                                                                      !
!     Module:       GlobalData                                                                                         !
!                                                                                                                      !
!     Description:  Contains all global variables used by the model                                                    !
!                                                                                                                      !
!======================================================================================================================!
module GlobalData
  implicit none
  public
  save

  ! double precision kind parameter
  integer, parameter :: dp=kind(0.d0)
  ! 4 bytes integer
  integer, parameter :: i4=kind(1)

  ! tstart - time of simulation start
  real(kind=dp) :: tstart
  ! tend = time of simulation end
  real(kind=dp) :: tend
  ! dt - integration time step (seconds) 
  real(kind=dp) :: dt
  ! dthalf - half of integration time step (seconds)
  real(kind=dp) :: dthalf
  ! dtenv - time step for input environmental data (seconds)
  real(kind=dp) :: dtenv
  ! ntenv - total number of environmental time steps for simulation
  integer(kind=i4) :: ntenv
  ! dtout - time step for output (seconds)
  !       MUST be even multiple of dt
  real(kind=dp) :: dtout
  ! ntout - total number of output time steps
  integer(kind=i4) :: ntout
  ! zadeg - zenith angle at current meteorological data time step (degrees)
  real(kind=dp) :: zadeg
  ! zarad - zenith angle at current meteorological data time step (radians)
  real(kind=dp) :: zarad

  ! t - current simulation time
  real(kind=dp) :: t
  ! nt - current output time step number
  integer(kind=i4) :: nt
  ! nte - current met data time step number
  integer(kind=i4) :: nte
  ! zcurrent - current z value (only valid within chemistry integration!) (m)
  real(kind=dp) :: zcurrent
  ! zaltm - z-altitude for box model simulation (m)
  real(kind=dp) :: zaltm
  ! zpt - current vertical grid point (only valid within chemistry integration!)
  integer(kind=i4) :: zpt

  ! number of grid nodes
  integer(kind=i4), parameter :: npts = 1
  ! number of grid nodes in canopy
  integer(kind=i4)            :: ncnpy

  ! vertical grid definition (cm)
  real(kind=dp), dimension(npts) :: z

!$INSERT MECHANISM-SPECIFIC SPECIES PARAMETERS

  ! Gas-phase photolytic rate constant data
  ! nphoto - number of J values calculated in CalcJPhoto
  integer(kind=4), parameter :: nphoto = 51
  ! kphoto - gas-phase photolytic rate constants
  real(kind=dp), dimension(nphoto) :: kphoto

  ! Photolysis reaction indices
  !
  ! ozone
  !   O3 + hv -> O(1D) + O2
  integer(kind=i4), parameter :: jo3a      = 1
  !   O3 + hv -> O(3P) + O2
  integer(kind=i4), parameter :: jo3b      = 2
  ! hydrogen peroxide
  !   H2O2 + hv -> 2OH
  integer(kind=i4), parameter :: jh2o2     = 3
  ! nitrogen dioxide
  !   NO2 + hv -> NO + O(3P)
  integer(kind=i4), parameter :: jno2      = 4
  ! nitrate radical
  !   NO3 + hv -> NO + O2
  integer(kind=i4), parameter :: jno3a     = 5
  !   NO3 + hv -> NO2 + O(3P)
  integer(kind=i4), parameter :: jno3b     = 6
  ! nitrous acid
  !   HONO + hv -> NO + OH
  integer(kind=i4), parameter :: jhono     = 7
  ! nitric acid
  !   HNO3 + hv -> NO2 + OH
  integer(kind=i4), parameter :: jhno3     = 8
  ! peroxynitric acid
  !   HO2NO2 + hv -> NO2 + HO2
  integer(kind=i4), parameter :: jho2no2   = 9
  ! formaldehyde
  !   HCHO + hv -> HCO + H
  integer(kind=i4), parameter :: jhchoa    = 10
  !   HCHO + hv -> CO + H2
  integer(kind=i4), parameter :: jhchob    = 11
  ! acetaldehyde
  !   CH3CHO + hv -> CH3 + HCO
  integer(kind=i4), parameter :: jch3cho   = 12
  ! propionaldehyde
  !   C2H5CHO + hv -> C2H5 + HCO
  integer(kind=i4), parameter :: jc2h5cho  = 13
  ! n-butyraldehyde
  !   n-C3H7CHO + hv -> n-C3H7 + HCO
  integer(kind=i4), parameter :: jc3h7choa = 14
  !   n-C3H7CHO + hv -> CH3CHO + CH2=CH2
  integer(kind=i4), parameter :: jc3h7chob = 15
  ! isobutyraldehyde
  !   i-C3H7CHO + hv -> i-C3H7 + HCO
  integer(kind=i4), parameter :: jic3h7cho = 16
  ! methacrolein
  !   MACR + hv -> CH2=C(CH3) + HCO
  integer(kind=i4), parameter :: jmacra    = 17
  !   MACR + hv -> CH2+C(CH3)CO + H
  integer(kind=i4), parameter :: jmacrb    = 18
  ! acetone
  !   CH3COCH3 + hv -> CH3CO + CH3
  integer(kind=i4), parameter :: jacet     = 19
  ! methyl ethyl ketone
  !   MEK + hv -> CH3CO + CH2CH3
  integer(kind=i4), parameter :: jmek      = 20
  ! methyl vinyl ketone
  !   MVK + hv -> CH3CH=CH2 + CO
  integer(kind=i4), parameter :: jmvka     = 21
  !   MVK + hv -> CH3CO + CH2=CH2
  integer(kind=i4), parameter :: jmvkb     = 22
  ! glyoxal
  !   CHOCHO + hv -> 2CO + H2
  integer(kind=i4), parameter :: jglyoxa   = 23
  !   CHOCHO + hv -> HCHO + CO
  integer(kind=i4), parameter :: jglyoxb   = 24
  !   CHOCHO + hv -> 2HCO
  integer(kind=i4), parameter :: jglyoxc   = 25
  ! methyl gloxal
  !   MGLYOX + hv -> CH3CO + HCO
  integer(kind=i4), parameter :: jmglyox   = 26
  ! biacetyl
  !   BIACET + hv -> 2CH3CO
  integer(kind=i4), parameter :: jbiacet   = 27
  ! methyl hydroperoxide
  !   CH3OOH + hv -> CH3O + OH
  integer(kind=i4), parameter :: jch3ooh   = 28
  ! methyl nitrate
  !   CH3NO3 + hv -> CH3O + NO2
  integer(kind=i4), parameter :: jch3no3   = 29
  ! ethyl nitrate
  !   C2H5NO3 + hv -> C2H5O + NO2
  integer(kind=i4), parameter :: jc2h5no3  = 30
  ! n-propyl nitrate
  !   n-C3H7NO3 + hv -> n-C3H7O + NO2
  integer(kind=i4), parameter :: jnc3h7no3 = 31
  ! isopropyl nitrate
  !   i-C3H7NO3 + hv -> i-C3H7O + NO2
  integer(kind=i4), parameter :: jic3h7no3 = 32
  ! t-butyl nitrate
  !   t-C4H9NO3 + hv -> t-C4H9O + NO2
  integer(kind=i4), parameter :: jtc4h9no3 = 33
  ! 1-hydroxy-2-propanone nitrate
  !   CH3C(O)CH2NO3 + hv -> CH3C(O)CH2O + NO2
  integer(kind=i4), parameter :: jnoaa     = 34
  !   CH3C(O)CH2NO3 + hv -> CH3CO + HCHO + NO2
  integer(kind=i4), parameter :: jnoab     = 35
  ! glycolaldehyde
  !   HOCH2CHO + hv -> HCO + HCHO + HO2
  integer(kind=i4), parameter :: jhoch2cho = 36
  ! hydroxyacetone (acetol)
  !   ACETOL + hv -> CH3CO + HCHO + HO2
  integer(kind=i4), parameter :: jacetol   = 37
  ! peroxyacetic acid
  !   CH3C(O)OOH + hv -> CH3O2 + OH
  integer(kind=i4), parameter :: jpaa      = 38
  ! pyruvic acid
  !   CH3COCOOH + hv -> CH3CO3 + HO2
  integer(kind=i4), parameter :: jpyra     = 39
  ! hydroperoxyaldehyde
  !   HPALD + hv -> products
  integer(kind=i4), parameter :: jhpald    = 40
  ! acrolein
  !   CH2=CHCHO + hv -> products
  integer(kind=i4), parameter :: jacro     = 41
  ! hydroxy methyl hydroperoxide
  !   HOCH2OOH + hv -> HOCH2O + OH
  integer(kind=i4), parameter :: jhmhp     = 42
  ! nitroxy ethanol
  !   CH2(OH)CH2(ONO2) + hv -> CH2(OH)CH2(O) + NO2
  integer(kind=i4), parameter :: jneth     = 43
  ! nitroxy acetone
  !   CH3COCH2(ONO2) + hv -> CH3COCH2(O) + NO2
  integer(kind=i4), parameter :: jnacet    = 44
  ! peroxy acetyl nitrate
  !   PAN + hv -> CH3CO(OO) + NO2
  integer(kind=i4), parameter :: jpana     = 45
  ! peroxy acetyl nitrate
  !   PAN + hv -> CH3CO(O) + NO3
  integer(kind=i4), parameter :: jpanb     = 46
  ! peroxy propionyl nitrate
  !   PPN + hv -> CH3CH2CO(OO) + NO2
  integer(kind=i4), parameter :: jppna     = 47
  ! peroxy propionyl nitrate
  !   PPN + hv -> CH3CH2CO(O) + NO3
  integer(kind=i4), parameter :: jppnb     = 48
  ! dinitrogen pentoxide
  !   N2O5 + hv -> NO3 + NO + O(3P)
  integer(kind=i4), parameter :: jn2o5a    = 49
  ! dinitrogen pentoxide
  !   N2O5 + hv -> NO3 + NO2
  integer(kind=i4), parameter :: jn2o5b    = 50
  ! benzaldehyde
  !   BenzALD + hv -> HCO + HO2 + CO
  integer(kind=i4), parameter :: jbald     = 51

  ! q - emissions rate of all species at the current time (molecules/cm3-s)
  real(kind=dp), dimension(npts, ninteg) :: q
  ! vd - dry deposition exchange coefficient (cm/s)
  real(kind=dp), dimension(npts, ninteg) :: vd
  ! rb - leaf boundary layer resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rb
  ! rm - leaf mesophyll resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rm
  ! rc - leaf cuticular resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rc
  ! rs - leaf stomatal resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rs

  ! gp - dry deposition compensation point :: (molecules/cm3)
  real(kind=dp), dimension(npts, ninteg) :: gp

  ! ctot - total solution array for current integration timestep (molecules/cm3)
  real(kind=dp), dimension(npts, ntotal) :: ctot
  ! ctotp - total solution array for previous integration timestep (molecules/cm3)
  real(kind=dp), dimension(npts, ntotal) :: ctotp

  ! cint - integrated species solution array for current integration timestep (??)
  real(kind=dp), dimension(npts, ninteg) :: cint
  ! cintp - integrated species solution array for previous integration timestep (??)
  real(kind=dp), dimension(npts, ninteg) :: cintp

  ! cfxed - fixed species array (molecules/cm3)
  real(kind=dp), dimension(npts, nfixed) :: cfxed

  ! css - steady state species solution array for current integration time step (molecules/cm3)
  real(kind=dp), dimension(npts, nsstate) :: css
  ! cssp - steady state species solution array for previous integration time step (molecules/cm3)
  real(kind=dp), dimension(npts, nsstate) :: cssp

  ! cout - saved total solution array over entire simulation (molecules/cm3)
  real(kind=dp), allocatable :: cout(:,:,:)
  ! vdout - saved deposition velocities array over entire simulation (cm/s)
  real(kind=dp), allocatable :: vdout(:,:,:)

  ! rbout - saved leaf boundary layer resistance (s/cm)
  real(kind=dp), allocatable :: rbout(:,:,:)
  ! rcout - saved leaf cuticular resistance (s/cm)
  real(kind=dp), allocatable :: rcout(:,:,:)
  ! rmout - saved leaf mesophyll resistance (s/cm)
  real(kind=dp), allocatable :: rmout(:,:,:)
  ! rsout - saved leaf stomatal resistance (s/cm)
  real(kind=dp), allocatable :: rsout(:,:,:)
  ! rsoillout - saved soil resistance of chemical species (s/cm)
  real(kind=dp), allocatable :: rsoillout(:,:)

  ! qout - saved emissions array over entire simulation (molecules/cm3-s)
  real(kind=dp), allocatable :: qout(:,:,:)
  ! emtout - saved budget emissions rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: emtout(:,:,:)
  ! cntout - saved constrained species rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: cntout(:,:,:)
  ! dptout - saved budget deposition rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: dptout(:,:,:)
  ! chtout - saved budget chemical production/loss rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: chtout(:,:,:)
  ! rxnout - saved reaction rates for entire simulation (molecules/cm3-s)
  real(kind=dp), allocatable :: rxnout(:,:,:)
  ! ksout - saved reaction rate coefficients for entire simulation (molecules-cm-s units)
  real(kind=dp), allocatable :: ksout(:,:,:)
  ! timeout - saved t values corresponding to output saved in cout
  real(kind=dp), allocatable :: timeout(:)
  ! sdtout - saved sdatetime values corresponding to output saved in cout
  character(len=19), allocatable :: sdtout(:)

  ! Calculated met data saved for output
  ! tkout - saved temperature profiles (K)
  real(kind=dp), allocatable :: tkout(:,:)
  ! pmbout - saved pressure profiles (mb)
  real(kind=dp), allocatable :: pmbout(:,:)
  ! qhout - saved specific humidity profiles (g/kg)
  real(kind=dp), allocatable :: qhout(:,:)
  ! rhout - saved relative humidity profiles (%)
  real(kind=dp), allocatable :: rhout(:,:)
  ! fjout - saved attentuated radiative flux profiles ()
  real(kind=dp), allocatable :: fjout(:,:)
  ! cairout - saved cair profiles (molecules/cm3)
  real(kind=dp), allocatable :: cairout(:,:)
  ! h2oout - saved h2o profiles (molecules/cm3)
  real(kind=dp), allocatable :: h2oout(:,:)

  ! Meteorological data

  ! tk - vertical profile of temperature at current simulation time (K)
  real(kind=dp), dimension(npts) :: tk
  ! tkn - vertical profile of temperature at next dtenv time step (K)
  real(kind=dp), dimension(npts) :: tkn
  ! tkp - vertical profile of temperature at previous dtenv time step (K)
  real(kind=dp), dimension(npts) :: tkp
  ! qh - vertical profile of specific humidity at current simulation time (g/kg)
  real(kind=dp), dimension(npts) :: qh 
  ! qhn - vertical profile of specific humidity at next dtenv time step (g/kg)
  real(kind=dp), dimension(npts) :: qhn
  ! qhp - vertical profile of specific humidity at previous dtenv time step (g/kg)
  real(kind=dp), dimension(npts) :: qhp
  ! pmb - vertical profile of pressure at current simulation time  (mb)
  real(kind=dp), dimension(npts) :: pmb
  ! pmbn - vertical profile of pressure at next dtenv time step (mb)
  real(kind=dp), dimension(npts) :: pmbn
  ! pmbp - vertical profile of pressure at previous dtenv time step (mb)
  real(kind=dp), dimension(npts) :: pmbp
  ! fj - vertical profile of actinic flux at current simulation time (??)
  ! TODO: initially this will only be the actinic flux normalized to the
  ! canopy top (i.e., dimensionless factor applied to Jphoto rates)
  real(kind=dp), dimension(npts) :: fj
  ! fjn - vertical profile of actinic flux at next dtenv time step (??)
  real(kind=dp), dimension(npts) :: fjn
  ! fjp - vertical profile of actinic flux at previous dtenv time step (??)
  real(kind=dp), dimension(npts) :: fjp

  ! Calculated meteorological data
  ! cair - vertical air density at current simulation time (molec/cm3)
  real(kind=dp), dimension(npts) :: cair
  ! h2o - vertical water vapor concentration at current simulation time (molec/cm3)
  real(kind=dp), dimension(npts) :: h2o

  ! Constrained species data
  ! kcns - mixing rate for constrained species (1/s)
  real(kind=dp)    :: kcns
  ! ncns - number of constrained species
  integer(kind=i4) :: ncns
  ! cnssp - indices of constrained species
  integer(kind=i4), allocatable :: cnssp(:)
  ! cnsu - flag for ppbv to molec/cm3 conversion
  integer(kind=i4), allocatable :: cnsu(:)
  ! ccns - constrained species concentrations (molecules/cm3)
  real(kind=dp), allocatable    :: ccns(:,:)

  ! Emissions fluxes data
  ! tlk - vertical profile of leaf temperature at current simulation time (K)
  real(kind=dp), dimension(npts) :: tlk

  ! TODO:  Figure out how to make these dynamically set at runtime
  ! fevergreen - fraction of forest which is composed of evergreen trees
  real(kind=dp), parameter    :: fevergreen=0.1
  ! scale_bvoc - scaling factor for bvoc emissions
  real(kind=dp), parameter    :: scale_bvoc=0.2

  ! Physical-Chemical data
  ! mdiffstp - molecular diffusivities of species in air at 0degC and 1 atm (cm2/s)
  real(kind=dp), dimension(ntotal) :: mdiffstp
  ! mdiffstp_default - default value of mdiffstp (cm2/s) for species with no reliable data
  real(kind=dp), parameter         :: mdiffstp_default=0.100
  ! hstar - effective Henry's Law coefficients (M/atm)
  real(kind=dp), dimension(ntotal) :: hstar
  ! hstar_default - default value of hstar (M/atm)
  real(kind=dp), parameter         :: hstar_default=1.000
  ! f0 - reactivity parameter used in resistance calculation (dimensionless)
  real(kind=dp), dimension(ntotal) :: f0
  ! f0_default - default value of f0
  real(kind=dp), parameter         :: f0_default=0.0
  ! molecmass - molecular mass of simulated species
  real(kind=dp), dimension(ntotal) :: molecmass

  ! simulation control file unit number
  integer(kind=i4), parameter :: UCTRL=11
  ! simulation control file name
  character(len=14), parameter :: filectrl='accessCTRL.dat'
  ! initial conditions file unit number
  integer(kind=i4), parameter :: UICS=12
  ! initial conditions file name
  character(len=16) :: icfile
  ! met data file unit number
  integer(kind=i4), parameter :: UENV=13
  ! met data file name
  character(len=16) :: envfile

  ! simdescr - descriptive simulation identification
  character(len=100) :: simdescr
  ! simname - unique simulation name
  character(len=16) :: simname

  ! slat - latitude (deg)
  real(kind=dp) :: slat
  ! slon - longitude (deg)
  real(kind=dp) :: slon

  ! year - year of simulation
  integer(kind=i4) :: year
  ! month - month of simulation (at start)
  integer(kind=i4) :: month
  ! daymonth - day of the month (at start)
  integer(kind=i4) :: daymonth
  ! doy - day of year at beginning of simulation
  integer(kind=i4) :: doy
  ! hz, mz, sz - local time at beginning of simulation
  integer(kind=i4) :: hz, mz, sz
  ! tzone - time zone differential
  integer(kind=i4) :: tzone

  ! ODESOLVER selection
  integer(kind=i4)            :: ODESOLVER
  integer(kind=i4), parameter :: ODEDVODE   =1     ! DVODE - Default solver
  character(len=8), dimension(1) :: sode
  data sode /'DVODE  '/

  ! unit number for output CTRL simulation summary
  integer(kind=i4), parameter :: USUMM=20
  ! name of CTRL simulation summary file
  character(len=32) :: simsummary
  ! unit number for runtime output file
  integer(kind=i4), parameter :: URUN=21
  ! name of runtime output file
  character(len=32) :: simrunfile
  ! begining output file unit number for species files
  integer(kind=i4), parameter :: UOUT=30

  ! stdsp - indices of species printed to STDOUT
  integer(kind=i4) :: nstdsp
  integer(kind=i4), allocatable :: stdsp(:)

  ! options
  ! OPT_SIMTYPE - dynamic or steady-state simulation
  integer(kind=i4)            :: OPT_SIMTYPE
  integer(kind=i4), parameter :: DYNAMIC=1      ! solar time changes with t
  integer(kind=i4), parameter :: SSTATE=2       ! solar time fixed at tstart
  ! AQPHASE - option for a aqueous chemistry
  logical :: AQPHASE
  ! CHEMISTRY - option for chemistry
  logical :: CHEMISTRY
  ! CONSTRAIN - option for constrained species
  logical :: CONSTRAIN 
  ! EMISSION - option for emissions
  logical :: EMISSION
  ! DRYDEPOS - option for dry deposition
  logical :: DRYDEPOS

  ! rtol, atol - ode relative and absolute tolerances
  real(kind=dp), parameter :: rtol=1.0D-02
  real(kind=dp), parameter :: atol=1.0D+01

  ! Avogadro's number
  real(kind=dp), parameter :: navo=6.022D+23

  ! pi
  real(kind=dp), parameter :: pi=4.0*atan(1.0)

end module GlobalData
!======================================================================================================================!
