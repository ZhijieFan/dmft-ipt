program hmipt
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                :: i,iloop,L
  logical                :: converged
  complex(8)             :: zeta
  real(8)                :: n,wband,wmix,dt,tmax,dw
  real(8),allocatable    :: wr(:),time(:)
  complex(8),allocatable :: sigma(:),fg(:),fg0(:),sold(:),sigt(:),sigw(:)

  call parse_input_variable(L,"L",'inputIPT.in',default=10000)
  call parse_input_variable(Wband,'wband','inputIPT.in',default=1d0)
  call parse_input_variable(wmix,'wmix','inputIPT.in',default=1d0)
  call read_input("inputIPT.in")
  allocate(fg(L),sigma(L),fg0(L),wr(L),sold(L))

  wr=linspace(-wmax,wmax,L,mesh=dw)
  dt= pi/wmax
  tmax=dt*L/2
  allocate(time(L),sigt(L))
  time = linspace(-tmax,tmax-dt,L)

  sigma=zero ; iloop=0 ; converged=.false.       
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta  = cmplx(wr(i),eps) - sigma(i)
        fg(i) = gfbether(wr(i),zeta,wband)
     enddo
     fg0 = one/(one/fg + sigma)
     sold = sigma
     sigma= ipt_solve_real(fg0,wr)
     sigt = f_fft_rw2rt(sigma)*dw/pi2
     sigw = f_fft_rt2rw(sigt)*dt
     call splot("Sigma_t.ipt",time,sigt)
     call splot("Sigma_w.ipt",wr,sigw)
     sigma = wmix*sigma + (1.d0-wmix)*sold
     converged=check_convergence(sigma,dmft_error,nsuccess,nloop)
  enddo
  call splot("G_realw.ipt",wr,-dimag(fg)/pi,dreal(fg))
  call splot("G0_realw.ipt",wr,fg0)
  call splot("Sigma_realw.ipt",wr,sigma)

end program hmipt

