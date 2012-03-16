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

program hmmpt_2dsquare
  USE DMFT_IPT
  USE COMMON
  USE SQUARE_LATTICE
  USE IOTOOLS
  implicit none
  integer    :: i,Lk
  real(8)    :: x(1),z
  logical    :: check,converged
  complex(8) :: zeta

  real(8) :: nread,nerror,ndelta
  integer :: nindex=0

  real(8),allocatable    :: wt(:),epsik(:),wm(:)
  complex(8),allocatable :: sold(:)
  type(matsubara_gf)     :: fgm
  include "revision.inc"
  call version(revision)

  call read_input("inputIPT.in")
  call parse_cmd_variable(nread,"NREAD",default=0.d0)
  call parse_cmd_variable(nerror,"NERROR",default=1.d-4)
  call parse_cmd_variable(ndelta,"NDELTA",default=0.1d0)

  allocate(fg(L),sigma(L),fg0(L),gamma(L))
  allocate(sold(L))
  allocate(wr(L))

  !Build frequency array
  wr=linspace(-wmax,wmax,L,mesh=fmesh)

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)


  n=0.5d0 ; sigma=u*(n-0.5d0) ; xmu0=xmu 
  iloop=0 ; converged=.false. ; sold=zero
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5,1x)")"DMFT-loop",iloop
     do i=1,L
        zeta  = cmplx(wr(i),eps) + xmu - sigma(i)
        fg(i) = sum_overk_zeta(zeta,epsik,wt)
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
     call splot("nVSiloop.ipt",iloop,n,append=TT)
     if(nread/=0.d0)call search_mu(converged)
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

contains

  subroutine search_mu(convergence)
    real(8)               :: naverage
    logical,intent(inout) :: convergence
    real(8)               :: ndelta1
    integer               :: nindex1    
    naverage=2.d0*n
    nindex1=nindex
    ndelta1=ndelta
    if((naverage >= nread+nerror))then
       nindex=-1
    elseif(naverage <= nread-nerror)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0)then !avoid loop forth and back
       ndelta=ndelta1/2.d0 !decreasing the step
    else
       ndelta=ndelta1
    endif
    xmu=xmu+real(nindex,8)*ndelta
    write(*,"(A,f15.12,A,f15.12,A,f20.17,A,f15.12)")" n=",naverage," /",nread,&
         "| shift=",nindex*ndelta,"| xmu=",xmu
    write(*,"(A,f15.12)")"dn=",abs(naverage-nread)
    print*,""
    if(abs(naverage-nread)>nerror)convergence=.false.
    call splot("muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
  end subroutine search_mu

end program hmmpt_2dsquare







