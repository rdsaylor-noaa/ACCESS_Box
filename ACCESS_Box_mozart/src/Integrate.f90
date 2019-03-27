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
!     Module:       Integrate                                                                                          !
!                                                                                                                      !
!     Description:  Contains all top-level routines to integrate the model                                             !
!                                                                                                                      !
!======================================================================================================================!
module Integrate
  use GlobalData
  use EnvironData
  use Emissions
  use DryDep
  use Chem
  use GasChem
  use Utils
  implicit none

  private IntegrateEmissions, &
          IntegrateConstrained, &
          IntegrateChemistry, &
          IntegrateDryDeposition, &
          UpdateSolutionArrays
  public IntegrateModel

contains

!**********************************************************************************************************************!
! subroutine IntegrateModel - handles overall model integration for one output
!                             time step
!**********************************************************************************************************************!
subroutine IntegrateModel(tin, tout)
  real(kind=dp), intent(in) :: tin, tout
  real(kind=dp) :: ts             ! current integration time [tin, tout]

  ts = tin

  zadeg = SolarZenithAngle(ts/(24.*3600.))
  zarad = zadeg*pi/180.

  ! read the next time slice of met data
  ! if interpolation is not enabled, met data at tin is used over the
  ! entire time step
  call ReadEnvironData()

  ! loop through integration sequence until tout
  ! dt is the integration time step
  ! Strang operator splitting is used
  do

    if (EMISSION  .eqv. .true.) then
      call IntegrateEmissions(ts, dt)                  ! emissions 
    end if

    if (CONSTRAIN .eqv. .true.) then
      call IntegrateConstrained(ts, dthalf)            ! constrained species
    end if

    if (CHEMISTRY .eqv. .true.) then
      call IntegrateChemistry(ts, dt)                  ! chemistry of all phases
    end if

    if (CONSTRAIN .eqv. .true.) then
      call IntegrateConstrained(ts+dthalf, dthalf)     ! constrained species
    end if

    if (DRYDEPOS  .eqv. .true.) then
      call IntegrateDryDeposition(ts, dt)              ! dry deposition
    end if

    ts = ts+dt

    call UpdateSolutionArrays()

    if (ts >= tout) exit
 
  end do

  return
end subroutine IntegrateModel

!**********************************************************************************************************************!
! subroutine IntegrateEmissions - adds emissions of all emitted species over tdstep
!                                 emissions are allowed at all vertical grid points
!**********************************************************************************************************************!
subroutine IntegrateEmissions(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  real(kind=dp)    :: dzi, dz1, dzn
  integer(kind=i4) :: i, l
   
  ! obtains emissions, q (molecules/cm3-s)
  call GetEmissions()

  ! adds emissions for all integrated species
  do l=1,ninteg
    do i=1,npts
      ! integration over dtstep
      cint(i,l) = cint(i,l) + q(i,l)*dtstep

      ! accumulating budget emissions (ppbv)
      emtout(i,l,nt)=emtout(i,l,nt) + q(i,l)*dtstep*1.0D+9/cair(i)
    end do
  end do

  return
end subroutine IntegrateEmissions

!**********************************************************************************************************************!
! subroutine IntegrateConstrained - integrates constrained species term
!**********************************************************************************************************************!
subroutine IntegrateConstrained(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  integer(kind=i4) :: i, l

  kcns=0.1

  ! integration for dtstep
  do l=1,ncns
    do i=1,npts
      ! integration over dtstep
      cint(i,cnssp(l)) = cint(i,cnssp(l)) - kcns*dtstep*(cint(i,cnssp(l))-ccns(i,l))
      cint(i,cnssp(l)) = max(cint(i,cnssp(l)), 0.0)
      
      ! accumulating budget constrained species process (ppbv)
      cntout(i,cnssp(l),nt)=cntout(i,cnssp(l),nt) &
                           - kcns*dtstep*(cint(i,cnssp(l))-ccns(i,l))*1.0D+09/cair(i)
    end do
  end do

  return
end subroutine IntegrateConstrained

!**********************************************************************************************************************!
! subroutine IntegrateChemistry - loop over each vertical grid point and integrate
!                                  chemistry in each layer from ts to ts+dtstep
!**********************************************************************************************************************!
subroutine IntegrateChemistry(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  integer(kind=i4) :: ierror
  integer(kind=i4) :: i, l, j
  real(kind=dp), parameter         :: rmin = 0.0D0
! real(kind=dp), parameter         :: rmin = 1.0D-20
  real(kind=dp), dimension(ninteg) :: y
  real(kind=dp), dimension(nrxn)   :: k
  real(kind=dp), dimension(nrxn)   :: r
   
  ! loop over each vertical grid point
  do i=1,npts

    ! set current z value (in meters) and current vertical grid point
    zcurrent = zaltm
    zpt = i
   
    ! at grid point "i", map species concentrations to array "y"
    ! Units:  molecules/cm3
    do l=1,ninteg
      y(l)=cint(i,l)

      ! accumulating budget chemistry process (ppbv)
      chtout(i,l,nt)=chtout(i,l,nt) - y(l)*1.0D+9/cair(i)
    end do

    ! integrate chemistry at grid point "i" from "ts" to "ts+dtstep"
    call IntegChemOneLevel(ninteg, y, ts, ts+dtstep, ierror)

    ! get reaction rates for output later
    call GasRateCoeffs(ts, k)
    call RxnRates(ninteg, k, y, r)
    do j=1,nrxn
      ksout(i,j,nt)=k(j)                ! molecules-cm-s units
      rxnout(i,j,nt)=max(r(j), rmin)    ! molecules/cm3-s
    end do

    ! map "y" back to "cint" 
    do l=1,ninteg
      cint(i,l) = y(l)

      ! accumulating budget chemistry process (ppbv)
      ! equivalent to adding increment (yafter-ybefore)
      chtout(i,l,nt)=chtout(i,l,nt) + y(l)*1.0D+9/cair(i)
    end do

  end do   ! end of vertical grid loop

  return
end subroutine IntegrateChemistry

!**********************************************************************************************************************!
! subroutine IntegrateDryDeposition - for each species removes mass via dry deposition
!**********************************************************************************************************************!
subroutine IntegrateDryDeposition(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  integer(kind=i4) :: i, l
  real(kind=dp) :: cnew
  real(kind=dp) :: ladbox

  ! get deposition velocities
  call GetDryDepExCoeffs()

  ladbox=0.0
  ! loop over all species and levels, removing mass via dry deposition
  ! do not allow concentrations to go < 0 (can only remove what's there to start)
  do l=1,ninteg
    do i=2,npts-1
      cnew = cint(i,l) - (cint(i,l)-gp(i,l))*vd(i,l)*ladbox*dtstep
      cnew = max(cnew, 0.0)

      ! accumulating budget dry deposition process (ppbv)
      dptout(i,l,nt)=dptout(i,l,nt) + (cnew - cint(i,l))*1.0D+09/cair(i)
      cint(i,l) = cnew
    end do
  end do

  return
end subroutine IntegrateDryDeposition

!**********************************************************************************************************************!
! subroutine UpdateSolutionArrays - updates all storage arrays after successful t-step
!**********************************************************************************************************************!
subroutine UpdateSolutionArrays()
  integer(kind=i4) :: i, l

  do l=1,ninteg
    do i=1,npts
      ctot(i,l) = cint(i,l)
    end do
  end do

  return
end subroutine UpdateSolutionArrays

end module Integrate
!======================================================================================================================!
