program hmipt
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  logical                :: converged
  real(8)                :: n,dw,Wband,A,dt
  integer                :: i,L,Lk,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)
  real(8),allocatable    :: wr(:),t(:)
  type(keldysh_gf)       :: kfg

  call parse_input_variable(L,"L",'inputIPT.in',default=10000)
  call parse_input_variable(Wband,'wband','inputIPT.in',default=1d0)
  call read_input("inputIPT.in")

  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(wr(L))

  !build freq. array
  wr  = linspace(-wmax,wmax,L,mesh=dw)
  sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !GET GLOC:
     fg=zero
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i)
        fg(i) = gfbether(wr(i),zeta,Wband)
     enddo
     fg0     = one/(one/fg + sigma)
     sigma   = ipt_solve_keldysh(fg0)
     converged=check_convergence(sigma,dmft_error,nsuccess,nloop)
  enddo
  call splot("G0_realw.ipt",wr,fg0)
  call splot("G_realw.ipt",wr,fg)
  call splot("Sigma_realw.ipt",wr,sigma)


  call allocate_gf(kfg,L)
  kfg%ret%w = fg
  !
  allocate(t(0:L))
  t   = linspace(-pi/wmax*L/2,pi/wmax*L/2,L+1,mesh=dt)
  do i=1,L
     A = -dimag(kfg%ret%w(i))/pi
     kfg%less%w(i) = pi2*xi*fermi(wr(i),beta)*A
     kfg%gtr%w(i)  = pi2*xi*(fermi(wr(i),beta)-1.d0)*A
  enddo
  call splot("Gret_w.ipt",wr,kfg%ret%w)
  call splot("Gless_w.ipt",wr,kfg%less%w)
  call splot("Ggtr_w.ipt",wr,kfg%gtr%w)
  !
  kfg%ret%t  =  fft_rw2rt(kfg%ret%w)*dw/pi2
  kfg%less%t =  fft_rw2rt(kfg%less%w)*dw/pi2
  kfg%gtr%t  =  fft_rw2rt(kfg%gtr%w)*dw/pi2
  call splot("Gret_t.ipt",t,kfg%ret%t)
  call splot("Gless_t.ipt",t,kfg%less%t)
  call splot("Ggtr_t.ipt",t,kfg%gtr%t)


end program hmipt
