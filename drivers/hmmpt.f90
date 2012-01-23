!########################################################
!     Program  : HMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
!     AUTHORS  : Adriano Amaricci
!########################################################
module COMMON
  USE BROYDEN
  implicit none
  real(8)                :: xmu0,n,n0
  real(8),allocatable    :: wr(:)
  complex(8),allocatable :: sigma(:),fg(:),fg0(:),gamma(:)
end module COMMON

function funcv(x)
  USE COMMON
  USE TOOLS
  USE COMMON_VARS
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0=x(1)
  !Hartree corrected WF is: 
  !\tilde{\calG0}^-1 = \calG0^-1 +xmu -xmu0 -U*n
  fg0 = one/(one/gamma +xmu-xmu0-U*(n-0.5d0))
  n0=sum(fermi(wr,beta)*aimag(fg0))/sum(aimag(fg0))
  funcv(1)=n-n0
  write(*,"(3(f13.9))")n,n0,xmu0
end function funcv

program hmmpt
  USE DMFT_IPT
  USE COMMON
  implicit none
  integer    :: i
  real(8)    :: x(1)
  logical    :: check,converged
  complex(8) :: zeta

  complex(8),allocatable :: sold(:)

  call read_input("inputIPT.in")
  allocate(fg(L),sigma(L),fg0(L),gamma(L))
  allocate(sold(L))
  allocate(wr(L))

  wr=linspace(-wmax,wmax,L,mesh=fmesh)

  n=0.5d0 ; sigma=u*(n-0.5d0) ; xmu0=xmu 
  iloop=0 ; converged=.false. ; sold=sigma
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5,1x)")"DMFT-loop",iloop
     do i=1,L
        zeta  = cmplx(wr(i),eps) + xmu - sigma(i)
        fg(i) = gfbether(wr(i),zeta,2.d0*ts)
     enddo
     n=sum(fermi(wr,beta)*aimag(fg))/sum(aimag(fg))

     !Get the regular Weiss Field: \Gamma(w+xmu)
     gamma= one/(one/fg + sigma) !gamma= cmplx(wr,eps)+xmu-sigma-one/fg

     !Fix the xmu0 w/ condition n0=n
     x(1)=xmu
     call broydn(x,check)
     xmu0=x(1)
     sigma= solve_mpt_sopt(fg0,wr,n,n0,xmu0)
     sigma=weigth*sigma + (1.d0-weigth)*sold
     sold=sigma
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
     call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
     call splot("G_realw.ipt",wr,fg,append=printf)
     call splot("G0_realw.ipt",wr,fg0,append=printf)
     call splot("Sigma_realw.ipt",wr,sigma,append=printf)
  enddo

end program hmmpt







