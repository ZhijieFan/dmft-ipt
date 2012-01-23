!########################################################
!     Program  : HMIPT
!     PURPOSE  : Solve the Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt_zerot
  USE IPT_VARS_GLOBAL
  implicit none

  integer                :: i,Lk,ikm
  logical                :: converged
  complex(8)             :: zeta,sqroot
  real(8)                :: sq,sig,ex
  type(real_gf)          :: fg0,sigma
  complex(8),allocatable :: fg(:)
  real(8),allocatable    :: wr(:)


  call read_input("inputIPT.in")
  allocate(fg(2*L))
  call allocate_gf(fg0,L)
  call allocate_gf(sigma,L)

  !grids:
  allocate(wr(2*L))
  fmesh=0.0005d0;wmax=fmesh*L
  forall(i=1:2*L,i<=L)wr(i)=real(i,8)*fmesh
  forall(i=1:2*L,i>L)wr(i)=real(i-2*L-1)*fmesh
  dt = pi/wmax

  D=1.d0 ; sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5,1x)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta  = wr(i) - sigma%w(i)
        sq=real(zeta,8) ; sig=1.d0 ; if(sq<0.d0)sig=-sig
        sqroot   = cdsqrt(zeta**2-one*D**2)
        fg(i)    = 2.d0/(zeta+sig*sqroot)
        fg0%w(i) = one/(one/fg(i) + sigma%w(i))
     enddo
     forall(i=0:L-1)fg0%w(2*L-i)=-fg0%w(i+1)

     call fftgf_rw2rt(fg0%w,fg0%t,L) ; fg0%t=fmesh/pi2*fg0%t
     forall(i=-L:L)sigma%t(i)=(U**2)*(fg0%t(i)**2)*fg0%t(-i)
     call fftgf_rt2rw(sigma%t,sigma%w,L) ; sigma%w= dt*sigma%w

     fg = one/(one/fg0%w - sigma%w)
     converged= check_convergence(sigma%w,eps_error,nsuccess,nloop)
  enddo
  call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
  call splot("Sigma_realw.ipt",wr,sigma%w,append=printf)
  call splot("G0_realw.ipt",wr,fg0%w,append=printf)


end program hmipt_zerot

