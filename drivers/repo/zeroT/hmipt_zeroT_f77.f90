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
  complex(8),allocatable :: fg(:),fg0(:),sigma(:),gc(:)
  complex(8),allocatable :: fg0t(:),sigmat(:),gct(:)
  real(8),allocatable    :: wr(:)


  call read_input("inputIPT.in")
  allocate(fg(2*L),wr(2*L),fg0(2*L),sigma(2*L),gc(2*L))
  allocate(fg0t(-L:L),gct(-L:L),sigmat(-L:L))

  !grids:
  fmesh=0.0005d0;wmax=fmesh*L
  forall(i=1:2*L,i<=L)wr(i)=real(i,8)*fmesh
  forall(i=1:2*L,i>L)wr(i)=real(i-2*L-1)*fmesh
  dt = pi/wmax


  D=1.d0
  sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5,1x)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta  = wr(i) - sigma(i)
        sig=1.d0
        sq=real(zeta,8)
        if(sq<0.d0)sig=-1.d0
        sqroot= cdsqrt(zeta**2-one*D**2)
        fg(i) = 2.d0/(zeta+sig*sqroot)
        fg0(i)  = one/(one/fg(i) + sigma(i))
     enddo
     !
     forall(i=0:L-1)fg0(2*L-i)=-fg0(i+1)
     fg0(2*L)=fg0(2*L-1)

     gc=fg0
     call four1(gc,2*L,1)
     !     starts manipulations of the arrays:
     !     1) [1,2*L=n]---> [-L,L]
     do i=1,2*L
        gct(i-L-1)=fmesh/2.d0/pi*gc(i)
     enddo
     do i=-L,-1
        fg0t(i+L)=gct(i)
     enddo
     do i=0,L-1
        fg0t(i-L)=gct(i)   
     enddo

     call splot("G0_t.ipt",(/(i,i=-L,L )/),fg0t)
     do i=-L+1,L-1
        sigmat(i)=-(U**2)*(fg0t(i)**2)*fg0t(-i)
     enddo
     sigmat(-L)=-(U**2)*(fg0t(L)**2)*fg0t(L-1)


     do i=1,2*L
        gc(i)=sigmat(i-L-1)
     enddo
     CALL four1(gc,2*L,1)
     ex=-1.d0
     do i=1,2*L
        ex=-ex
        sigma(i)=ex*dt*gc(i)
     enddo

     fg = one/(one/fg0 - sigma)
     fg(1)=-fg(2*L)

     converged= check_convergence(sigma,eps_error,nsuccess,nloop)
  enddo

  call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
  call splot("Sigma_realw.ipt",wr,sigma,append=printf)
  call splot("G0_realw.ipt",wr,fg0,append=printf)

end program hmipt_zerot

