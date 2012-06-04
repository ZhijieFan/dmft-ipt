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
  real(8)                :: xmu0(2),n(2),n0(2)
  type(matsubara_gf)     :: fg(2),fg0(2)
  complex(8),allocatable :: sold(:,:),sigma(:,:),gamma(:,:)
END MODULE COMMON

FUNCTION FUNCV(x)
  USE COMMON
  USE DMFT_IPT
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0(1)=x(1)
  xmu0(2)=x(2)
  !Hartree corrected WF is: \tilde{\calG0}^-1 = \calG0^-1 +xmu -xmu0 -U*n
  !####CHECK THIS#####
  fg0(1)%iw = one/(one/gamma(1,:) +xmu-xmu0(1)-u*(n(2)-0.5d0))
  fg0(2)%iw = one/(one/gamma(2,:) +xmu-xmu0(2)-u*(n(1)-0.5d0))

  call fftgf_iw2tau(fg0(1)%iw,fg0(1)%tau,beta)
  call fftgf_iw2tau(fg0(2)%iw,fg0(2)%tau,beta)

  n0(1)=-real(fg0(1)%tau(L))
  n0(2)=-real(fg0(2)%tau(L))
  funcv(1)=n(1)-n0(1)
  funcv(2)=n(2)-n0(2)
  write(*,"(6(f13.9))")n(1),n(2),n0(1),n0(2),xmu0(1),xmu0(2)
END FUNCTION FUNCV

program hmmpt_matsubara
  USE DMFT_IPT
  USE COMMON
  USE IOTOOLS
  implicit none
  real(8)    :: x(2),z(2),extH
  logical    :: check
  complex(8) :: zeta(2),zita
  logical :: converged
  integer :: i
  real(8),dimension(:),allocatable :: wm

  call read_input("inputIPT.in")
  call allocate_gf(fg,L)
  call allocate_gf(fg0,L)

  allocate(sigma(2,L))
  allocate(gamma(2,L))
  allocate(sold(2,L))
  allocate(wm(L))

  wm = pi/beta*real(2*arange(1,L)-1,8)

  n=0.5d0 ; xmu0=xmu 
  sigma(1,:)= U*(n(1)-0.5d0)
  sigma(2,:)= U*(n(2)-0.5d0)
  iloop=0 ; converged=.false.;D=1.d0
  extH=0.1d0
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)")"DMFT-loop",iloop
     do i=1,L
        zeta(1)    = xi*wm(i) + xmu +extH - sigma(1,i)
        zeta(2)    = xi*wm(i) + xmu -extH - sigma(2,i)
        zita=sqrt(zeta(1)*zeta(2))
        fg(1)%iw(i) = gfbethe(wm(i),zita,D)
        fg(2)%iw(i) = gfbethe(wm(i),zita,D)
     enddo
     fg(1)%iw(i) =  fg(1)%iw(i)*sqrt(zeta(2)/zeta(1))
     fg(2)%iw(i) =  fg(2)%iw(i)*sqrt(zeta(1)/zeta(2))
     call fftgf_iw2tau(fg(1)%iw,fg(1)%tau,beta)
     call fftgf_iw2tau(fg(2)%iw,fg(2)%tau,beta)
     n(1)=-real(fg(1)%tau(L))
     n(2)=-real(fg(2)%tau(L))
     extH=0.d0

     !Get the regular Weiss Field: \Gamma(w+xmu)
     gamma(1,:)= one/(one/fg(1)%iw + sigma(1,:))
     gamma(2,:)= one/(one/fg(2)%iw + sigma(2,:))

     !Fix the xmu0 w/ condition n0=n
     x(:)=xmu
     call broydn(x,check)
     xmu0=x

     gamma(1,:)=fg0(1)%iw ; gamma(2,:)=fg0(2)%iw
     sigma = solve_mpt_af_matsubara(gamma,n,n0,xmu0)
     sigma = weight*sigma + (1.d0-weight)*sold
     sold=sigma
     converged=check_convergence(sigma,eps_error,Nsuccess,Nloop)
     ! z=1.d0 - dimag(sigma%iw(1))/wm(1);z=1.d0/z
     call splot("nVSiloop.ipt",iloop,n(1),n(2),append=TT)
     ! call splot("zetaVSiloop.ipt",iloop,z,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  ! call close_file("zetaVSiloop.ipt")
  call splot("G_iw.ipt",wm,fg(1)%iw,append=printf)
  call splot("G_iw.ipt",wm,fg(2)%iw,append=TT)
  call splot("G0_iw.ipt",wm,gamma(1,:),append=printf)
  call splot("G0_iw.ipt",wm,gamma(2,:),append=TT)
  call splot("Sigma_iw.ipt",wm,sigma(1,:),append=printf)
  call splot("Sigma_iw.ipt",wm,sigma(2,:),append=TT)
  stop
end program hmmpt_matsubara






