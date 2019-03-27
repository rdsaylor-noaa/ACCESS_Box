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
!     Module:       Initialize                                                                                         !
!                                                                                                                      !
!     Description:  Contains routines related to model initialization                                                  !
!                                                                                                                      !
!======================================================================================================================!
module Initialize
  use GlobalData
  use EnvironData
  use PhysChemData
  use DryDep
  use Utils
  use Output
  implicit none

  private SetInitialConditions, SetSimulationData
  public InitializeModel

contains

!**********************************************************************************************************************!
! InitializeModel - performs all the model initialization steps 
!                   Called from Main.f90
!**********************************************************************************************************************!
subroutine InitializeModel()
  real(kind=dp)    :: delta, leap, mjd, hrlocal, hrutc, time21
  integer(kind=i4) :: m
  integer(kind=i4), dimension(12) :: months

  data months /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

  ! read CTRL file
  call SetSimulationData()

  ! calculate doy (day-of-year)
  doy = 0.0
  do m=1,month-1
    doy = doy + dble(months(m))
  end do
  doy = doy + dble(daymonth)

  ! calculate decimal UTC hour
  hrlocal = dble(hz) + (dble(mz)/60.) + (dble(sz)/3600.)
  hrutc = hrlocal - tzone

  ! years from 1949
  delta = dble(year)-1949.

  ! leap days over this period (note correction below)
  leap = aint(delta/4.)

  ! modified Julian Day
  mjd=32916.5+delta*365.+leap+dble(doy)+hrutc/24.

  ! the last year of century is not a leap year unless divisible by 400
  if (dmod(dble(year),100.0D+00) .eq. 0.0 .and. dmod(dble(year),400.0D+00) .ne. 0.0) then
    mjd=mjd-1.
  end if

  ! time in days since Jan 1, 2000 Greenwich Noon
  time21 = mjd-51545.0 
 
  ! time in seconds since Jan 1, 2000 Greenwich Noon
  ! simulation time is kept in seconds
  tstart=time21*24.*3600.

  zadeg = SolarZenithAngle(time21)
  zarad = zadeg*pi/180.

  ! ending time
  tend=tstart+dtout*dble(ntout)
! tend=tstart+dtenv*dble(ntenv)
! ntout = int((tend-tstart)/dtout)
! print *, 'ntout=',ntout

  dthalf = 0.5_dp*dt                   ! half time step
  t=tstart                             ! simulation time start
  nt=0                                 ! output number start
  nte=0                                ! met data number start

  ! allocate storage for final output arrays
  allocate(cout(npts,ntotal,0:ntout))
  allocate(vdout(npts,ninteg,0:ntout))
  allocate(qout(npts,ninteg,0:ntout))
  allocate(emtout(npts,ninteg,0:ntout))
  allocate(cntout(npts,ninteg,0:ntout))
  allocate(dptout(npts,ninteg,0:ntout))
  allocate(chtout(npts,ninteg,0:ntout))
  allocate(ksout(npts,nrxn,0:ntout))
  allocate(rxnout(npts,nrxn,0:ntout))
  allocate(timeout(0:ntout))

  allocate(rbout(npts,ninteg,0:ntout))
  allocate(rmout(npts,ninteg,0:ntout))
  allocate(rcout(npts,ninteg,0:ntout))
  allocate(rsout(npts,ninteg,0:ntout))
  allocate(rsoillout(ninteg,0:ntout))

  allocate(sdtout(0:ntout))

  allocate(tkout(npts,0:ntout))
  allocate(pmbout(npts,0:ntout))
  allocate(qhout(npts,0:ntout))
  allocate(fjout(npts,0:ntout))
  allocate(cairout(npts,0:ntout))
  allocate(h2oout(npts,0:ntout))
  allocate(rhout(npts,0:ntout))

  ! make sure tendency arrays are zeroed initially
  emtout=0.0_dp
  cntout=0.0_dp
  dptout=0.0_dp
  chtout=0.0_dp
  ksout=0.0_dp
  rxnout=0.0_dp
  tkout=0.0_dp
  pmbout=0.0_dp
  qhout=0.0_dp
  fjout=0.0_dp
  cairout=0.0_dp
  h2oout=0.0_dp
  rhout=0.0_dp

  rbout=0.0_dp
  rmout=0.0_dp
  rcout=0.0_dp
  rsout=0.0_dp

  ! set all physical-chemical data
  call SetPhysChemData()

  ! read first time slice of environmental data file
  call ReadEnvironData()

  ! read initial conditions file
  call SetInitialConditions()

  return
end subroutine InitializeModel

!**********************************************************************************************************************!
! SetSimulationData - reads input data from CTRL file for simulation
!                     and writes to STDOUT and the simulation summary file
!**********************************************************************************************************************!
subroutine SetSimulationData()
  integer(kind=i4)  :: j, isp
  character(len=35) :: outdir
  character(len=35) :: strmkdir
  logical           :: eoutdir
  integer(kind=i4), parameter   :: nhlines=7

  ! read CTRL file
  open(unit=UCTRL,file=('./ctrl/' // filectrl))
  ! skip header lines
  do j=1,nhlines
    read(UCTRL,*)
  end do
  read(UCTRL,101) simdescr
  read(UCTRL,*) 
  read(UCTRL,102) simname 
  read(UCTRL,*) OPT_SIMTYPE
  read(UCTRL,*) slat
  read(UCTRL,*) slon
  read(UCTRL,*) year
  read(UCTRL,*) month
  read(UCTRL,*) daymonth
  read(UCTRL,*) hz, mz, sz
  read(UCTRL,*) tzone
  read(UCTRL,*) dt
  read(UCTRL,*) dtout
  read(UCTRL,*) ntout
  read(UCTRL,*) nstdsp
  allocate(stdsp(nstdsp))
  read(UCTRL,*) (stdsp(isp),isp=1,nstdsp)
  read(UCTRL,*) CHEMISTRY
  read(UCTRL,*) AQPHASE
  read(UCTRL,*) CONSTRAIN
  read(UCTRL,*) EMISSION
  read(UCTRL,*) DRYDEPOS
  read(UCTRL,*) ODESOLVER
  read(UCTRL,*)
  read(UCTRL,*) icfile
  read(UCTRL,*)
  read(UCTRL,*) envfile
  close(UCTRL)

  ! block unimplemented options
  if (AQPHASE) then
    write(6,1010)
    stop
  end if

  if (EMISSION) then
    write(6,1020)
    stop
  end if

  if (DRYDEPOS) then
    write(6,1030)
    stop
  end if

  if (ODESOLVER .ne. 1) then
    write(6,1040)
    stop
  end if

1010 format(/'AQPHASE option is not implemented yet!' / 'Stopping ...')
1020 format(/'EMISSION option is not implemented yet!' / 'Stopping ...')
1030 format(/'DRYDEPOS option is not implemented yet!' / 'Stopping ...')
1040 format(/'ODESOLVER option other than 1 (DVODE) is not available yet!' / 'Stopping ...')

100 format(a6)
101 format(a100)
102 format(a16)

  ! write CTRL data to STDOUT
  write(6,300)
  write(6,*) simdescr
  write(6,2010) simname 
  write(6,2011) gasmech
  if (OPT_SIMTYPE .eq. DYNAMIC) then
    write(6,2015)
  else
    write(6,2016)
  end if
  write(6,2020) slat
  write(6,2030) slon
  write(6,2031) year
  write(6,2040) month
  write(6,2041) daymonth
  write(6,2050) hz, mz, sz
  write(6,2051) tzone
  write(6,2060) dt
  write(6,2072) dtout
  write(6,2080) nstdsp
  write(6,2090) (trim(sspc(stdsp(isp))),isp=1,nstdsp)
  write(6,2200) 
  write(6,2110) CHEMISTRY
  write(6,2120) AQPHASE
  write(6,2140) CONSTRAIN
  write(6,2160) EMISSION
  write(6,2170) DRYDEPOS
  write(6,2179) sode(ODESOLVER)
  write(6,2210) 
  write(6,2180) icfile
  write(6,2190) envfile
  write(6,300)

  ! if output directory in 'out' does not exist, create it
  ! and all necessary subdirectories
  outdir='./out/' // trim(simname) // '/.'
  inquire(file=outdir, exist=eoutdir)
  if(eoutdir .eqv. .FALSE.) then
    strmkdir = 'mkdir ./out/' // trim(simname)
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/sp'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/emis'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/ks'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/rates'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/met'
    call system(strmkdir)
  end if

  ! open simulation runtime file to capture STDOUT
  ! runtime file is normally closed by Utils:CleanUp
  ! at the end of the simulation
  simrunfile='./out/' // trim(simname) // '/ACCESS_STD.out'
  simrunfile=trim(simrunfile)
  open(unit=URUN,file=simrunfile)

  ! write header for runtime file
  write(URUN,300)
  write(URUN,*) simdescr
  write(URUN,2010) simname 
  write(URUN,300)

  ! write CTRL data to USUMM 
  simsummary='./out/' // trim(simname) // '/ACCESS_SUMM.out'
  simsummary=trim(simsummary)
  open(unit=USUMM,file=simsummary)
  write(USUMM,300)
  write(USUMM,*) simdescr
  write(USUMM,2010) simname 
  write(USUMM,2011) gasmech
  if (OPT_SIMTYPE .eq. DYNAMIC) then
    write(USUMM,2015)
  else
    write(USUMM,2016)
  end if
  write(USUMM,2020) slat
  write(USUMM,2030) slon
  write(USUMM,2031) year
  write(USUMM,2040) month
  write(USUMM,2041) daymonth
  write(USUMM,2050) hz, mz, sz
  write(USUMM,2051) tzone
  write(USUMM,2060) dt
  write(USUMM,2072) dtout
  write(USUMM,2080) nstdsp
  write(USUMM,2090) (trim(sspc(stdsp(isp))),isp=1,nstdsp)
  write(USUMM,2200) 
  write(USUMM,2110) CHEMISTRY
  write(USUMM,2120) AQPHASE
  write(USUMM,2140) CONSTRAIN
  write(USUMM,2160) EMISSION
  write(USUMM,2170) DRYDEPOS
  write(USUMM,2179) sode(ODESOLVER)
  write(USUMM,2210) 
  write(USUMM,2180) icfile
  write(USUMM,2190) envfile
  write(USUMM,300)
  close(USUMM)

2000 format(/'Starting ACCESS v', a, ' simulation...'/)
2010 format(' Short sim name = ', a)
2011 format(' Chemistry mech = ', a)
2015 format(' OPT_SIMTYPE    = DYNAMIC')
2016 format(' OPT_SIMTYPE    = SSTATE')
2020 format(' Latitude       = ', f7.2, ' deg')
2030 format(' Longitude      = ', f7.2, ' deg')
2029 format(' z altitude     = ', f10.1, ' m')
2031 format(' Year           = ', i4)
2040 format(' Month          = ', i4)
2041 format(' Day            = ', i4)
2050 format(' Start time     = ', i2.2, ':', i2.2, ':', i2.2, ' LT')
2051 format(' Time zone diff = ', i3)
2060 format(' Integration dt = ', f7.1, ' s')
2070 format(' Met data dt    = ', f7.1, ' s')
2071 format(' # metdata stps = ', i4)
2072 format(' Output dt      = ', f7.1, ' s')
2080 format(/' Species to STDOUT = ', i4)
2090 format(100(1x,a))
2110 format(' CHEMISTRY      = ', l2)
2120 format(' AQPHASE        = ', l2)
2140 format(' CONSTRAIN      = ', l2)
2160 format(' EMISSION       = ', l2)
2170 format(' DRYDEPOS       = ', l2)
2179 format(' ODESOLVER      = ', a)
2180 format(' IC file name   = ', a)
2190 format(' ENVDATA file name  = ', a)
2200 format(/' Model Options:')
2210 format(/' Input Files:')
300 format(80('='))
301 format(' ACCESS v',a)
8000 format('  !!!!ERROR!!!!!!'/'  ACCESS version ',6a /'  CTRL version   ',6a &
           /'  Stopping ...')  

  return
end subroutine SetSimulationData

!**********************************************************************************************************************!
! SetInitialConditions - reads initial conditions file and sets ICs
!**********************************************************************************************************************!
subroutine SetInitialConditions()
  real(kind=dp)     :: ic_z0, ic_zi
  integer(kind=i4)  :: i, l, j, m
  integer(kind=i4)  :: ic_npts
  integer(kind=i4)  :: nics
  character(len=16) :: ic_gasmech
  integer(kind=i4), parameter   :: niclines=9
  integer(kind=i4), allocatable :: icsp(:), icu(:), icns(:)
  real(kind=dp), allocatable :: ival(:)

  ! initially set all concentrations to zero
  ! only those specified in IC file are non-zero
  do l=1,ninteg
    cint(npts,l) = 0.0D+00
  end do

  ! read IC file
  open(unit=UICS,file=('./data/' // icfile))
  ! skip header lines
  do j=1,niclines
    read(UICS,*)
  end do
  ! read consistency data and check
  read(UICS,101) ic_gasmech
  ic_gasmech=trim(ic_gasmech)
  if (ic_gasmech .ne. gasmech) then
    write(*,201) icfile
    write(*,202) ic_gasmech, gasmech
    close(UICS)
    stop
  end if

  ! made it this far, so read the data
  read(UICS,*)
  read(UICS,*)
  read(UICS,*)
  read(UICS,*) nics
  read(UICS,*)
  read(UICS,*)
  allocate(icsp(nics))
  allocate(icu(nics))
  allocate(ival(nics))
  allocate(icns(nics))
  do l=1,nics
    read(UICS,*) icsp(l), ival(l), icu(l), icns(l)
  end do

  close(UICS)

  ncns = 0
  do l=1,nics
    if(icns(l) .eq. 1) ncns=ncns+1
  end do
  allocate(ccns(npts,ncns))
  allocate(cnssp(ncns))
  allocate(cnsu(ncns))

  ! convert the ppbv species (where icu(l)=1) to molec/cm3
  do l=1,nics
    cint(npts,icsp(l)) = ival(l)
    if (icu(l) .eq. 1) then
      cint(npts,icsp(l)) = ConvertPPBToMolecCC(cint(npts,icsp(l)),npts)
    end if
  end do
  cintp = cint

  ! constrained species
  m = 0
  do l=1,nics
    if (icns(l) .eq. 1) then
      m=m+1
      cnssp(m) = icsp(l)
      ccns(npts,m) = cint(npts,icsp(l))
    end if
  end do

  deallocate(icsp)
  deallocate(icu)
  deallocate(ival)
  deallocate(icns)

101 format(a16)
201 format('***Initial conditions file ', a, ' is inconsistent with current ACCESS configuration!')
202 format('***IC gasmech = ', a /'***ACCESS gasmech = ', a)

  return
end subroutine SetInitialConditions

end module Initialize
!======================================================================================================================!
