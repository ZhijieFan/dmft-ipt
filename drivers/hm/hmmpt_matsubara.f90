!########################################################
!     Program  : HMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
!     AUTHORS  : Adriano Amaricci
!########################################################
MODULE COMMON
  USE GREENFUNX
  USE BROYDEN
  implicit none
  real(8)                :: xmu0,n,n0
  type(matsubara_gf)     :: sigma,fg,fg0,gamma
  complex(8),allocatable :: sold(:)
END MODULE COMMON

FUNCTION FUNCV(x)
  USE COMMON
  USE DMFT_IPT
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  real(8)                          :: zn0
  xmu0=x(1)
  !Hartree corrected WF is: \tilde{\calG0}^-1 = \calG0^-1 +xmu -xmu0 -U*n
  fg0%iw = one/(one/gamma%iw +xmu-xmu0-u*(n-0.5d0))
  call fftgf_iw2tau(fg0%iw,fg0%tau,beta)
  n0=-real(fg0%tau(L))
  funcv(1)=n-n0
  write(*,"(3(f13.9))")n,n0,xmu0
END FUNCTION FUNCV

program hmmpt_matsubara
  USE DMFT_IPT
  USE COMMON
  USE IOTOOLS
  implicit none
  real(8)    :: x(1),z
  logical    :: check
  complex(8) :: zeta
  logical :: converged
  integer :: i
  real(8),dimension(:),allocatable :: wm

  call read_input("inputIPT.in")
  call allocate_gf(fg,L)
  call allocate_gf(fg0,L)
  call allocate_gf(sigma,L)
  call allocate_gf(gamma,L)
  allocate(sold(size(sigma%iw)))
  allocate(wm(L))

  wm = pi/beta*real(2*arange(1,L)-1,8)

  n=0.5d0 ; sigma%iw= U*(n-0.5d0) ; xmu0=xmu 
  iloop=0 ; converged=.false.;D=1.d0
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)")"DMFT-loop",iloop
     do i=1,L
        zeta    = xi*wm(i) + xmu - sigma%iw(i)
        fg%iw(i) = gfbethe(wm(i),zeta,D)
     enddo
     call fftgf_iw2tau(fg%iw,fg%tau,beta)
     n=-real(fg%tau(L))

     !Get the regular Weiss Field: \Gamma(w+xmu)
     gamma%iw= one/(one/fg%iw + sigma%iw)

     !Fix the xmu0 w/ condition n0=n
     x(1)=xmu
     call broydn(x,check)
     xmu0=x(1)
     sigma%iw = solve_mpt_matsubara(fg0%iw,n,n0,xmu0)
     sigma%iw = weigth*sigma%iw + (1.d0-weigth)*sold ; sold=sigma%iw
     converged=check_convergence(sigma%iw,eps_error,Nsuccess,Nloop)
     z=1.d0 - dimag(sigma%iw(1))/wm(1);z=1.d0/z
     call splot("nVSiloop.ipt",iloop,n,append=TT)
     call splot("zetaVSiloop.ipt",iloop,z,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  call close_file("zetaVSiloop.ipt")
  call splot("G_iw.ipt",wm,fg%iw,append=printf)
  call splot("G0_iw.ipt",wm,fg0%iw,append=printf)
  call splot("Sigma_iw.ipt",wm,sigma%iw,append=printf)
end program hmmpt_matsubara






