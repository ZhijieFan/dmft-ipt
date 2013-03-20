!########################################################
!     Program  : HMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
!     AUTHORS  : Adriano Amaricci
!########################################################
module COMMON
  USE DMFT_IPT
  USE ZEROS
  implicit none
  real(8)                :: xmu0,n,n0
  real(8),allocatable    :: wr(:)
  complex(8),allocatable :: sigma(:),fg(:),fg0(:),gamma(:)
end module COMMON

function funcv(x)
  USE COMMON
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


subroutine gamma_func(in,x,fvec,iflag)
  USE COMMON
  implicit none
  integer               :: in ,iflag
  real(8),dimension(in) ::  x
  real(8),dimension(in) ::  fvec
  xmu0=x(1)
  !Hartree corrected WF is: 
  !\tilde{\calG0}^-1 = \calG0^-1 +xmu -xmu0 -U*n
  fg0 = one/(one/gamma +xmu-xmu0-U*(n-0.5d0))
  n0=sum(fermi(wr,beta)*aimag(fg0))/sum(aimag(fg0))
  fvec(1)=n-n0
  write(*,"(3(f13.9))")n,n0,xmu0
end subroutine gamma_func


program hmmpt
  USE COMMON
  USE IOTOOLS
  implicit none
  integer    :: i,iloop
  real(8)    :: x(1),z
  logical    :: check,converged
  complex(8) :: zeta
  complex(8),allocatable :: sold(:)
  real(8),allocatable    :: wm(:)
  type(matsubara_gf)     :: fgm
  external gamma_func
  include "revision.inc"
  call version(revision)

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
     ! x(1)=xmu
     ! call broydn(x,check)
     ! xmu0=x(1)
     call fsolve(gamma_func,x,tol=1.d-15)
     xmu0=x(1)

     sigma= solve_mpt_sopt(fg0,wr,n,n0,xmu0)
     sigma=weight*sigma + (1.d0-weight)*sold
     sold=sigma
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
     call splot("nVSiloop.ipt",iloop,n,append=.true.)
  enddo
  call close_file("nVSiloop.ipt")
  call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
  call splot("G_realw.ipt",wr,fg,append=printf)
  call splot("G0_realw.ipt",wr,fg0,append=printf)
  call splot("Sigma_realw.ipt",wr,sigma,append=printf)


  L=max(8192,8*L)
  call allocate_gf(fgm,L)
  allocate(wm(L))
  wm = pi/beta*real(2*arange(1,L)-1,8)

  call get_matsubara_gf_from_dos(wr,fg,fgm%iw,beta)
  call splot("G_iw.ipt",wm,fgm%iw,append=printf)
  fgm=zero
  call get_matsubara_gf_from_dos(wr,sigma,fgm%iw,beta)
  call splot("Sigma_iw.ipt",wm,fgm%iw,append=printf)
  z=1.d0 - dimag(fgm%iw(1))/wm(1);z=1.d0/z
  call splot("n.z.u.xmu.ipt",n,z,u,xmu,append=printf)

end program hmmpt







