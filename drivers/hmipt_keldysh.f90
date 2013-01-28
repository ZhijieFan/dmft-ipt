program hmipt
  !########################################################
  !     Program  : HMIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model using DMFT-IPT
  !     AUTHORS  : Adriano Amaricci
  !########################################################
  !LOCAL:
  USE DMFT_IPT
  USE IOTOOLS
  implicit none

  logical                :: converged
  real(8)                :: n,dw!,A
  integer                :: i,Lk,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)
  real(8),allocatable    :: wr(:)
  !real(8),allcoatable    :: t(:)
  !type(keldysh_equilibrium_gf) :: fg0k,sk

  call read_input("inputIPT.in")
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(wr(L))
  ! allocate(t(-L/2:L/2))
  ! call allocate_gf(fg0k,L)
  ! call allocate_gf(sk,L)

  !build freq. array
  wr  = linspace(-wmax,wmax,L,mesh=dw)
  ! dt  = pi/wmax
  ! t   = linspace(-dt*L/2,dt*L/2,L+1,mesh=dt)

  sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     !GET GLOC:
     fg=zero
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i)
        fg(i) = gfbether(wr(i),zeta,2.d0*ts)
     enddo
     fg0     = one/(one/fg + sigma)

     !Update Keldysh components of the Weiss Field:
     ! do i=1,L
     !    A = -dimag(fg0(i))/pi
     !    fg0k%less%w(i) = pi2*xi*fermi(wr(i),beta)*A
     !    fg0k%gtr%w(i)  = pi2*xi*(fermi(wr(i),beta)-1.d0)*A
     ! enddo
     ! call fftgf_rw2rt(fg0k%less%w,fg0k%less%t,L) ; fg0k%less%t=dw/pi2*fg0k%less%t
     ! call fftgf_rw2rt(fg0k%gtr%w,fg0k%gtr%t,L)   ; fg0k%gtr%t =dw/pi2*fg0k%gtr%t
     ! do i=-L/2,L/2
     !    sk%less%t(i)=(U**2)*(fg0k%less%t(i)**2)*fg0k%gtr%t(-i) 
     !    sk%gtr%t(i) =(U**2)*(fg0k%gtr%t(i)**2)*fg0k%less%t(-i)
     !    sk%ret%t(i) =heaviside(t(i))*(sk%gtr%t(i)-sk%less%t(i))
     ! enddo
     ! if(heaviside(0.d0)==1.d0)sk%ret%t(0)=sk%ret%t(0)/2.d0 
     ! call fftgf_rt2rw(sk%ret%t,sk%ret%w,L) ; sigma=dt*sk%ret%w

     sigma   = solve_ipt_keldysh(fg0,wmax)
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
     call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
     call splot("G0_realw.ipt",wr,fg0,append=printf)
     call splot("Sigma_realw.ipt",wr,sigma,append=printf)
  enddo

end program hmipt
