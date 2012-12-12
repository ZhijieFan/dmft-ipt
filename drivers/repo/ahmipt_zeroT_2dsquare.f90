!########################################################
!     Program  : HMIPT
!     PURPOSE  : Solve the Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt_zerot_2dsquare
  USE IPT_VARS_GLOBAL
  USE SQUARE_LATTICE
  implicit none

  integer,parameter      :: M=2048
  integer                :: i,Lk,ik
  logical                :: converged
  complex(8)             :: zeta1,zeta2,det
  real(8)                :: sq,sig,n
  complex(8),allocatable :: fg(:,:),fg0(:,:),w0(:,:)
  type(real_gf)          :: calG11,calG22,calF11,calF22,sigma(2)
  real(8),allocatable    :: wr(:),wt(:),epsik(:),t(:)

  call read_input("inputIPT.in")
  allocate(fg(2,2*L),fg0(2,2*L),w0(2,2*L))
  call allocate_gf(calG11,L)
  call allocate_gf(calG22,L)
  call allocate_gf(calF11,L)
  call allocate_gf(calF22,L)
  call allocate_gf(sigma,L)


  !grids:
  allocate(wr(2*L),t(-L:L))
  wr = linspace(-wmax,wmax,2*L,mesh=fmesh)
  dt = pi/wmax
  t  = linspace(-dt*dble(L),dt*dble(L),2*L+1)

  !
  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)
  call get_free_dos(epsik,wt,wmin=-wmax,wmax=wmax,eps=0.005d0)


  sigma=zero ; sigma(2)%w=-deltasc ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     !call start_loop(iloop,nloop,"DMFT-loop")
     write(*,"(A,i5,1x)",advance="no")"DMFT-loop",iloop

     fg=zero
     do i=1,2*L
        zeta1 =   wr(i) + xmu - sigma(1)%w(i)
        zeta2 =   wr(i) - xmu + (sigma(1)%w(2*L+1-i))
        !sq=real(zeta1,8) 
        sig=1.d0 ; if(wr(i)<0.d0)sig=-sig
        zeta1=zeta1 + sig*xi*eps
        zeta2=zeta2 + sig*xi*eps
        do ik=1,Lk
           det = (zeta1-epsik(ik))*(zeta2+epsik(ik)) - (sigma(2)%w(i))**2
           fg(1,i)=fg(1,i) +  wt(ik)*(zeta2+epsik(ik))/det
           fg(2,i)=fg(2,i) -  wt(ik)*(sigma(2)%w(i))/det
        enddo
     enddo
     n       =   sum(aimag(fg(1,:))*fermi(wr,beta=1000.d0))*fmesh/pi
     deltasc =-u*sum(aimag(fg(2,:))*fermi(wr,beta=1000.d0))*fmesh/pi
     call splot("G_realw.ipt",wr,fg(1,:),append=printf)
     call splot("DOS.ipt",wr,abs(aimag(fg(1,:)))/pi,append=printf)
     call splot("F_realw.ipt",wr,fg(2,:),append=printf)


     do i=1,2*L
        det     = fg(1,i)*(fg(1,2*L+1-i)) + (fg(2,i))**2
        fg0(1,i)= (fg(1,2*L+1-i))/det + sigma(1)%w(i)
        fg0(2,i)= (fg(2,i))/det + sigma(2)%w(i) - deltasc
     end do
     call splot("g0_realw.ipt",wr,fg0(1,:),append=printf)
     call splot("f0_realw.ipt",wr,fg0(2,:),append=printf)

     do i=1,2*L
        det    =  (fg0(1,i)*(fg0(1,2*L+1-i))) - abs(fg0(2,i))**2
        w0(1,i)=  (fg0(1,2*L+1-i))/det
        w0(2,i)=  -conjg(fg0(2,i))/det
     enddo
     call splot("w011_realw.ipt",wr,w0(1,:),append=printf)
     call splot("w022_realw.ipt",wr,w0(2,:),append=printf)

     forall(i=1:2*L)
        calG11%w(i) = w0(1,i)
        calG22%w(i) = -w0(1,2*L+1-i)
        calF11%w(i) = w0(2,i)
        calF22%w(i) = (w0(2,i))
     end forall
     call splot("calG11_realw.ipt",wr,calG11%w,append=printf)
     call splot("calG22_realw.ipt",wr,calG22%w,append=printf)
     call splot("calF11_realw.ipt",wr,calF11%w,append=printf)
     call splot("calF22_realw.ipt",wr,calF22%w,append=printf)


     call fftgf_rw2rt(calG11%w,calG11%t,L) ; calG11%t=fmesh/pi2*calG11%t
     call fftgf_rw2rt(calG22%w,calG22%t,L) ; calG22%t=fmesh/pi2*calG22%t
     call fftgf_rw2rt(calF11%w,calF11%t,L) ; calF11%t=fmesh/pi2*calF11%t
     call fftgf_rw2rt(calF22%w,calF22%t,L) ; calF22%t=fmesh/pi2*calF22%t
     call splot("calG11_t.ipt",t,calG11%t,append=printf)
     call splot("calG22_t.ipt",t,calG22%t,append=printf)
     call splot("calF11_t.ipt",t,calF11%t,append=printf)
     call splot("calF22_t.ipt",t,calF22%t,append=printf)


     forall(i=-L:L)sigma(1)%t(i)= U**2*(calG11%t(i)*calG22%t(i) - calF11%t(i)*(calF22%t(i)))*calG22%t(-i)
     forall(i=-L:L)sigma(2)%t(i)= U**2*(calF11%t(i)*(calF22%t(i)) - calG11%t(i)*calG22%t(i))*calF11%t(-i)
     call fftgf_rt2rw(sigma(1)%t,sigma(1)%w,L) ; sigma(1)%w= dt*sigma(1)%w
     call fftgf_rt2rw(sigma(2)%t,sigma(2)%w,L) ; sigma(2)%w= dt*sigma(2)%w

     write(*,"(2(f16.12),1x)",advance="no")n,deltasc
     sigma(2)%w = sigma(2)%w - deltasc
     converged= check_convergence(sigma(1)%w+sigma(2)%w,eps_error,nsuccess,nloop)


     call splot("Sigma_realw.ipt",wr,sigma(1)%w,append=printf)
     call splot("Self_realw.ipt",wr,sigma(2)%w,append=printf)
     call splot("G0_realw.ipt",wr,fg0(1,:),append=printf)
     call splot("F0_realw.ipt",wr,fg0(2,:),append=printf)
     !call end_loop()
  enddo

end program hmipt_zerot_2dsquare

