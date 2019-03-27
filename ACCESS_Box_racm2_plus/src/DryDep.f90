!======================================================================================================================!
!                                                                                                                      !
!     Program:      ACCESS                                                                                             !
!                   Atmospheric Chemistry and Canopy Exchange Simulation System                                        !
!                   Full BVOC chemistry version                                                                        !
!                                                                                                                      !
!     Version:      3.1.0                                                                                              !
!                                                                                                                      !
!======================================================================================================================!
!                                                                                                                      !
!     Last Update:  Dec 2017                                                                                           !
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
!     Module:       DryDep                                                                                             !
!                                                                                                                      !
!     Description:  contains dry deposition velocities algorithms                                                      !
!                                                                                                                      !
!======================================================================================================================!
module DryDep
  use GlobalData
  use PhysChemData
  use Utils
  implicit none

  private  vdl
  public   GetDryDepExCoeffs

contains

!**********************************************************************************************************************!
! subroutine GetDryDepExCoeffs - calculate leaf-scale dry deposition velocities
!                                and compensation points
!**********************************************************************************************************************!
subroutine GetDryDepExCoeffs()
  integer(kind=i4)  :: i, l

  do i=1,npts
  do l=1,ninteg
    rb(i,l)=0.0_dp
    rc(i,l)=0.0_dp
    rm(i,l)=0.0_dp
    rs(i,l)=0.0_dp
    vd(i,l)=0.0_dp
    gp(i,l)=0.0_dp
  end do
  end do

  return
end subroutine GetDryDepExCoeffs

!**********************************************************************************************************************!
! function vdl - calculate dry deposition velocity
!
! Uses formulation suggested by Wolfe (2012) personal communication
!
!**********************************************************************************************************************!
function vdl(rb, rm, rc, rs)
  real(kind=dp), intent(in) :: rb    ! leaf boundary layer resistance (s/cm)
  real(kind=dp), intent(in) :: rm    ! mesophyll resistance (s/cm)
  real(kind=dp), intent(in) :: rc    ! cuticle resistance (s/cm)
  real(kind=dp), intent(in) :: rs    ! stomatal resistance (s/cm)
  real(kind=dp)             :: vdl   ! deposition velocity (cm/s)
  real(kind=dp)             :: rnum, rden, rl
  rnum = rc*(rs+rm)
  rden = rc+2.0*(rs+rm)
  rl = rb + (rnum/rden)
  vdl = 1.0/rl
  return 
end function vdl

end module DryDep
!======================================================================================================================!
