program hmipt
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  logical                          :: converged
  real(8)                          :: n,dw,A,dt,ts,depsi
  integer                          :: i,L,Lk,iloop,Ne
  complex(8)                       :: zeta
  complex(8),allocatable           :: Sigma(:),Gfunc(:),Weiss(:)
  real(8),allocatable              :: wr(:),t(:)
  real(8),dimension(:),allocatable :: epsi,Dos2d
  !
  character(len=24)                :: finput  
  type(keldysh_gf)                 :: kfg


  call parse_cmd_variable(finput,"finput",default="inputIPT.conf")
  call parse_input_variable(ts,"ts",finput,default=0.5d0)
  call parse_input_variable(L,"L",finput,default=10000)
  call parse_input_variable(Ne,"Ne",finput,default=1000)
  call read_input(finput)

  allocate(Gfunc(L))
  allocate(Sigma(L))
  allocate(Weiss(L))
  allocate(wr(L))


  allocate(epsi(Ne))
  allocate(dos2d(Ne))
  epsi = linspace(-wmax,wmax,Ne,mesh=depsi)
  do i=1,Ne
     Dos2d(i)=dens_2dsquare(epsi(i),ts)
  enddo
  call splot("Dos2d.ipt",epsi,Dos2d)
  Dos2d=Dos2d*depsi

  !build freq. array
  wr  = linspace(-wmax,wmax,L,mesh=dw)
  Sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !GET GLOC:
     Gfunc=zero
     do i=1,L
        zeta     = cmplx(wr(i),eps) - Sigma(i)
        Gfunc(i) = sum_overk_zeta(zeta,epsi,dos2d)
     enddo
     Weiss     = one/(one/Gfunc + Sigma)
     Sigma   = ipt_solve_keldysh(Weiss)
     converged=check_convergence(Sigma,dmft_error,nsuccess,nloop)
  enddo
  call splot("Weiss_R_realw.ipt",wr,Weiss)
  call splot("Gfunc_R_realw.ipt",wr,Gfunc)
  call splot("Sigma_R_realw.ipt",wr,Sigma)


  !Some manipulation of the Equilibrium Keldysh Green function.
  call allocate_gf(kfg,L)
  kfg%ret%w = Gfunc
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
