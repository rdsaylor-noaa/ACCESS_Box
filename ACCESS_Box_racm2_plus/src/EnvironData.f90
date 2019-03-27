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
!     Last Update:  July 2018                                                                                          !
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
!     Module:       EnvironData                                                                                        !
!                                                                                                                      !
!     Description:  contains algorithms related to reading and analyzing environmental data                            !
!                                                                                                                      !
!======================================================================================================================!
module EnvironData
  use GlobalData
  use PhysChemData
  use Utils
  implicit none

  public ReadEnvironData

  private

contains

!**********************************************************************************************************************!
! ReadEnvironData ... read one time step of environmental data
!                     NOTE: ReadEnvironData must be adapted to each new environmental data set, depending
!                           on availability of usable data at each simulated site
!
!**********************************************************************************************************************!
  subroutine ReadEnvironData()
    integer(kind=i4)          :: j, n
    integer(kind=i4), parameter :: nmhlines=15
    character(len=19)         :: sdt      ! time slice datetime as a string
    character(len=10)         :: sdate    ! time slice date as a string
    character(len=8)          :: stime    ! time slice time as a string
    real(kind=dp)             :: tac      ! air temperature at zref (C)
    real(kind=dp)             :: tak      ! air temperature at zref (K)
    real(kind=dp)             :: pmbm     ! air pressure at measurement height (mb)
    real(kind=dp)             :: rhm      ! relative humidity at measurement height (%)
    real(kind=dp)             :: fjm      ! actinic flux attenuation due to clouds ()
    real(kind=dp)             :: time21   ! simulation time (seconds since midnight Jan 1, 2000)
    logical                   :: initcall
    data initcall /.TRUE./
    save initcall

    ! initial call
    if (initcall) then
      open(unit=UENV, file=('./data/' // envfile))
      do j=1,nmhlines
        read(UENV,*)
      end do
      initcall = .FALSE.
    end if

    ! read environmental data for the next time slice
    read(UENV,*) sdate, stime, tac, pmbm, rhm, fjm

    sdt = sdate // ' ' // stime
    time21 = DateTime2SimTime(sdate, stime)
    zadeg = SolarZenithAngle(time21)
    zarad = zadeg*pi/180.

    sdtout(nt) = sdt

    ! unit conversions
    tak    = tac + 273.15

    ! temperature (K)
    tk(npts) = tak

    ! pressure (mb)
    pmb(npts) = pmbm

    ! estimate altitude above sea level (m) with the
    ! hypsometric equation
    zaltm = 153.846*tak*(((1013.25/pmbm)**(0.1902))-1.0)

    ! actinic flux attenuation 
    fj(npts) = fjm

    ! qh based on rhm and tk(npts)
    qh(npts) = SpecificHumidity(rhm, tk(npts), pmb(npts))

    ! Calculate air density (molecules/cm3) based on pmb and tk
    cair(npts) = pmb(npts)*7.2428D+18/tk(npts)

    ! Convert specific humidity to molecules/cm3
    h2o(npts) = Convert_qh_to_h2o(qh(npts), cair(npts))

    return
  end subroutine ReadEnvironData

end module EnvironData
