!########################################################
!     Program  : HMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
!     AUTHORS  : Adriano Amaricci
!########################################################
MODULE COMMON
  USE DMFT_IPT
  USE OPTIMIZE
  implicit none
  real(8)                             :: xmu0,n,n0,z,docc
  complex(8),dimension(:),allocatable :: sigma,fg,fg0,gamma,sold

contains

  function funcv(x)
    implicit none
    real(8),dimension(:),intent(in)  ::  x
    real(8),dimension(size(x))       ::  funcv
    real(8)                          ::  zn0
    xmu0=x(1)
    !Hartree corrected WF is: \tilde{\calG0}^-1 = \calG0^-1 +xmu -xmu0 -Uloc*n
    fg0 = one/(one/gamma +xmu-xmu0-uloc*(n-0.5d0))
    n0  = ipt_measure_dens_matsubara(fg0)
    funcv(1)=n-n0
    write(*,"(3(f13.9))")n,n0,xmu0
  END FUNCTION FUNCV

END MODULE COMMON



program hmmpt_matsubara
  USE DMFT_IPT
  USE DMFT_TOOLS
  USE SCIFOR
  USE COMMON
  implicit none
  real(8)                          :: x(1),D,wmix,energy(3)
  logical                          :: check
  complex(8)                       :: zeta
  logical                          :: converged
  integer                          :: i,iloop
  real(8),dimension(:),allocatable :: wm

  call parse_input_variable(D,"wband","inputIPT.in",default=1d0,comment="half-bandwidth, energy unit")
  call parse_input_variable(wmix,"wmix","inputIPT.in",default=0.75d0,comment="mixing parameter")
  call read_input("inputIPT.in")

  allocate(fg(L))
  allocate(fg0(L))
  allocate(sigma(L))
  allocate(gamma(L))
  allocate(sold(L))
  allocate(wm(L))

  wm = pi/beta*(2*arange(1,L)-1)

  n=0.5d0 ; sigma= Uloc*(n-0.5d0) ; xmu0=xmu 
  iloop=0 ; converged=.false.;D=1.d0
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)")"DMFT-loop",iloop
     do i=1,L
        zeta    = xi*wm(i) + xmu - sigma(i)
        fg(i) = gfbethe(wm(i),zeta,D)
     enddo

     !Get the regular Weiss Field: \Gamma(w+xmu)
     gamma = one/(one/fg + sigma)
     n     = ipt_measure_dens_matsubara(fg)

     !Fix the xmu0 w/ condition n0=n
     x(1)=xmu0
     call broydn(funcv,x,check)
     xmu0=x(1)
     sigma = mpt_solve_matsubara(fg0,n,n0,xmu0)
     if(iloop>1)sigma = wmix*sigma + (1.d0-wmix)*sold
     sold = sigma
     z    = ipt_measure_zeta_matsubara(sigma,fg0)
     docc = ipt_measure_docc_matsubara(sigma,fg0)
     n    = ipt_measure_dens_matsubara(fg)
     call splot("obserbables_all.ipt",n,docc,z,append=.true.)
     converged=check_convergence(sigma,dmft_error,Nsuccess,Nloop)
  enddo
  call splot("G_iw.ipt",wm,fg)
  call splot("G0_iw.ipt",wm,fg0)
  call splot("Sigma_iw.ipt",wm,sigma)
  call splot("obserbables_last.ipt",n,docc,z,energy(1),energy(2),energy(3))

end program hmmpt_matsubara






