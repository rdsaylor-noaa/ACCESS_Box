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
!     Module:       Output                                                                                             !
!                                                                                                                      !
!     Description:  contains routines related to output of results                                                     !
!                                                                                                                      !
!======================================================================================================================!
module Output
  use GlobalData
  use Utils
  implicit none

  private PrinttoSTDOUT, SaveResults
  public OutputResult, PrintFinaltoFile, PrintCPUtime
contains

!**********************************************************************************************************************!
! OutputResult - prints results to STDOUT and stores for later
!**********************************************************************************************************************!
subroutine OutputResult()

  ! print defined results to STDOUT and simulation runtime file
  call PrinttoSTDOUT()

  ! save results for printing at end of simulation
  call SaveResults()

  nt = nt + 1

  return
end subroutine OutputResult

!**********************************************************************************************************************!
! PrinttoSTDOUT - print selected results to STDOUT
!                 and simultaneously to simulation runtime file
!**********************************************************************************************************************!
subroutine PrinttoSTDOUT()
  integer(kind=i4)  :: i
  integer(kind=i4)  :: l
  character(len=2)  :: ncols
  character(len=35) :: f0
  character(len=30) :: f1
  character(len=30) :: f2
  character(len=4)  :: col1
  logical           :: initcall
  data initcall /.TRUE./
  save initcall

  ! dynamically create format string
  write(ncols,'(i2)') nstdsp
  f0='(1x, a, 2x, a, e12.3, 2x, a, e12.3)'
  f1='(4x,a4,1x,' // ncols // '(3x,a8,3x))'
  f2='(f8.1,' // ncols // 'e14.6)'
  col1='t-t0'

  ! array stdsp maps species indices for STDOUT output
  if (initcall) then
    write(*,fmt=f1) col1, (sspc(stdsp(l)),l=1,nstdsp)
    write(URUN,fmt=f1) col1, (sspc(stdsp(l)),l=1,nstdsp)
    initcall = .FALSE.
  end if

  ! -ConvertOneToOutputFormat determines the appropriate output units
  !  from data in CTRL file
  do i=npts,1,-1
    write(*,f2) (t-tstart), (ConvertOneToOutputFormat(cint(i,stdsp(l)),stdsp(l),i),l=1,nstdsp) 

    write(URUN,f2) (t-tstart), (ConvertOneToOutputFormat(cint(i,stdsp(l)),stdsp(l),i),l=1,nstdsp) 
  end do

  return
end subroutine PrinttoSTDOUT

!**********************************************************************************************************************!
! SaveResults - save all results to storage array
!**********************************************************************************************************************!
subroutine SaveResults()
  integer(kind=i4) :: l, i

  do l=1,ninteg
 
    ! concentrations, emissions and deposition
    do i=1,npts
      cout(i,l,nt) = cint(i,l)
      rsout(i,l,nt) = rs(i,l)
      rbout(i,l,nt) = rb(i,l)
      rmout(i,l,nt) = rm(i,l)
      rcout(i,l,nt) = rc(i,l)
      vdout(i,l,nt) = vd(i,l)
      qout(i,l,nt)  = q(i,l)
    end do
  end do

  do i=1,npts
    ! met data
    tkout(i,nt)    = tk(i)
    pmbout(i,nt)   = pmb(i)
    qhout(i,nt)    = qh(i)
    cairout(i,nt)  = cair(i)
    h2oout(i,nt)   = h2o(i) 
    rhout(i,nt)    = RelativeHumidity(tk(i), pmb(i), qh(i))
  end do

  timeout(nt)=(t-tstart)     ! leave as seconds since tstart

  return
end subroutine SaveResults

!**********************************************************************************************************************!
! PrintFinaltoFile - prints selected results to individual species files
!                    output species are defined in CTRL file
!**********************************************************************************************************************!
subroutine PrintFinaltoFile()
  integer(kind=i4) :: i, l, m, nstp, me, j
  character(len=3)  :: ncols
  character(len=30) :: f1, f2, f3, f4, f5, f6
  character(len=8)  :: srxn
  character(len=45) :: ofname
  character(len=40) :: hdr

  ! dynamically create format strings for output
  if( (ntout+1) < 100) then
    write(ncols,'(i2)') ntout+1
  else
    write(ncols,'(i3)') ntout+1
  end if
  f1='(6x,' // trim(ncols) // '(f8.1,5x))'
  f2='(' // trim(ncols) // '(e13.5))'
  f3='(a)'
  f4='(a,i5.5)'
  f5='(6a10)'
  f6='(2i10,4f10.4)'

  ! number of met data time steps per output time step
  nstp=1

  ! output all species profiles
  ofname='./out/' // trim(simname) // '/sp/sp.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    do i=1,npts
      write(UOUT,fmt=f2) (ConvertOneToOutputFormat(cout(i,l,m),l,i),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output emissions
  ofname='./out/' // trim(simname) // '/emis/emis.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    do i=1,npts
      write(UOUT,fmt=f2) (qout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output instantaneous reactions rate coefficients over simulation
  ofname='./out/' // trim(simname) // '/ks/ks.out'
  open(UOUT,file=ofname)
  do j=1,nrxn
    write(srxn,fmt=f4) 'rxn', j
    write(UOUT,fmt=f3) srxn
    do i=1,npts
      write(UOUT,fmt=f2) (ksout(i,j,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output instantaneous reactions rates over simulation
  ofname='./out/' // trim(simname) // '/rates/rates.out'
  open(UOUT,file=ofname)
  do j=1,nrxn
    write(srxn,fmt=f4) 'rxn', j
    write(UOUT,fmt=f3) srxn
    do i=1,npts
      write(UOUT,fmt=f2) (rxnout(i,j,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output temperature profiles over simulation
  ofname='./out/' // trim(simname) // '/met/tk.dat'
  open(UOUT,file=ofname)
  do i=1,npts
    write(UOUT,fmt=f2) (tkout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output pressure profiles over simulation
  ofname='./out/' // trim(simname) // '/met/pmb.dat'
  open(UOUT,file=ofname)
  do i=1,npts
    write(UOUT,fmt=f2) (pmbout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output specific humidity profiles over simulation
  ofname='./out/' // trim(simname) // '/met/qh.dat'
  open(UOUT,file=ofname)
  do i=1,npts
    write(UOUT,fmt=f2) (qhout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output radiation attenuation profiles over simulation
  ofname='./out/' // trim(simname) // '/met/fj.dat'
  open(UOUT,file=ofname)
  do i=1,npts
    write(UOUT,fmt=f2) (fjout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output air density profiles over simulation
  ofname='./out/' // trim(simname) // '/met/cair.dat'
  open(UOUT,file=ofname)
  do i=1,npts
    write(UOUT,fmt=f2) (cairout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output h2o profiles over simulation
  ofname='./out/' // trim(simname) // '/met/h2o.dat'
  open(UOUT,file=ofname)
  do i=1,npts
    write(UOUT,fmt=f2) (h2oout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output rh profiles over simulation
  ofname='./out/' // trim(simname) // '/met/rh.dat'
  open(UOUT,file=ofname)
  do i=1,npts
    write(UOUT,fmt=f2) (rhout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output elapsed hour/datetime key file
  ofname='./out/' // trim(simname) // '/ACCESS_timekey.dat'
  open(UOUT,file=ofname)
  hdr = 't-t0         datetime'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    me = m*nstp
    write(UOUT,1002) timeout(m), sdtout(me)
  end do
  close(UOUT)

  ! output species units key file
  ofname='./out/' // trim(simname) // '/ACCESS_ppbv.dat'
  open(UOUT,file=ofname)
  hdr = 'ispec        sspc'           
  write(UOUT,1000) trim(hdr)
  do i=1,noutppb
    write(UOUT,1003) outppb(i), sspc(outppb(i))
  end do
  close(UOUT)

1000 format(8x, a)
1001 format(f8.1, 5x, e13.5)
1002 format(f12.1, 5x, a)
1003 format(3x, i4, 10x, a)

  return
end subroutine PrintFinaltoFile

!**********************************************************************************************************************!
! subroutine PrintCPUtime - prints simulation CPU time to STDOUT and summary
!**********************************************************************************************************************!
subroutine PrintCPUtime(tcpu)
   real(kind=4) :: tcpu

   ! to STDOUT
   write(*,1001) tcpu/3600.,tcpu

   ! to USUMM for archive
   open(USUMM,file=simsummary,status='old',action='write',position='append')
   write(USUMM,1001) tcpu/3600.,tcpu
   close(USUMM)

1001 format('Total CPU time: ', f12.3,' hrs'/'                ',f12.3,' sec')

end subroutine PrintCPUtime

end module Output
!======================================================================================================================!
