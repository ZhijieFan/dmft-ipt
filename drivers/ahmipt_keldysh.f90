!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt
  USE DMFT_IPT 
  USE IOTOOLS
  implicit none
  integer    :: i,ik,Lk,iloop,Lm
  logical    :: converged
  complex(8) :: det,zeta1,zeta2
  real(8)    :: delta,n,A,B,nf,nfL,nfR
  !
  complex(8),allocatable          :: sigma(:,:),fg(:,:),wf0(:,:),calG(:,:),sold(:,:)
  complex(8),allocatable          :: zeta(:)
  type(matsubara_gf),dimension(2) :: gf,sf
  !
  real(8),allocatable             :: wt(:),epsik(:),wr(:),t(:)
  !
  type(keldysh_equilibrium_gf)    :: fg0k(2),sk(2),calG11,calG22,calF12,calF21
  !
  real(8) :: vbias

  include "revision.inc"
  call version(revision)
  call read_input("inputIPT.in")
  call parse_cmd_variable(vbias,'VBIAS',default=0.d0)

  allocate(fg(2,L))
  allocate(sigma(2,L))
  allocate(wf0(2,L))
  allocate(calG(2,L))
  allocate(zeta(L))
  allocate(sold(2,L))
  Lm=L/2
  call allocate_gf(fg0k(1),Lm);  call allocate_gf(fg0k(2),Lm)
  call allocate_gf(sk(1),Lm)  ;  call allocate_gf(sk(2),Lm)
  call allocate_gf(calG11,Lm)
  call allocate_gf(calG22,Lm)
  call allocate_gf(calF12,Lm)
  call allocate_gf(calF21,Lm)

  allocate(wr(L),t(-L/2:L/2))
  wr = linspace(-wmax,wmax,L,mesh=fmesh)
  dt  = pi/wmax
  t   = linspace(-dt*L/2,dt*L/2,L+1,mesh=dt)

  Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,2.d0*ts)

  call get_initial_sigma

  iloop=0    ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     !GET GLOC:
     fg=zero
     zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
     do i=1,L
        zeta1 = zeta(i)
        zeta2 = conjg(zeta(L+1-i))
        do ik=1,Lk
           det = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + conjg(sigma(2,L+1-i))*sigma(2,i)
           fg(1,i)=fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det
           fg(2,i)=fg(2,i) - wt(ik)*conjg(sigma(2,L+1-i))/det
        enddo
     enddo
     delta=-u*sum(dimag(fg(2,:))*fermi(wr,beta))*fmesh/pi
     n    =-sum(dimag(fg(1,:))*fermi(wr,beta))*fmesh/pi

     !GET THE WEISS FIELD \calG0^-1(w)
     ! wf(1)=calG0^-1 = G*(-w)/(G(w)G*(-w) + F(w)F(-w)) + Sigma(w)
     ! wf(2)=calF0^-1 =-F(-w) /(G(w)G*(-w) + F(w)F(-w)) + S(w)
     do i=1,L
        det     = fg(1,i)*conjg(fg(1,L+1-i)) + conjg(fg(2,L+1-i))*fg(2,i)
        wf0(1,i)= conjg(fg(1,L+1-i))/det  + sigma(1,i)   +  u*(n-0.5d0)
        wf0(2,i)= fg(2,i)/det       + conjg(sigma(2,L+1-i))   + delta
     end do
     do i=1,L
        det      =  wf0(1,i)*conjg(wf0(1,L+1-i)) + conjg(wf0(2,L+1-i))*wf0(2,i)
        calG(1,i)=  conjg(wf0(1,L+1-i))/det
        calG(2,i)=  conjg(wf0(2,L+1-i))/det
     end do


     fg0k(1)=zero ; fg0k(2)=zero
     do i=1,L
        A = -dimag(calG(1,i))/pi
        B = -dimag(calG(2,i))/pi
        nfL= fermi(wr(i)-vbias/2.d0,beta)
        nfR= fermi(wr(i)+vbias/2.d0,beta)
        nf= (nfL+nfR)/2.d0
        fg0k(1)%less%w(i) = pi2*xi*nf*A
        fg0k(1)%gtr%w(i)  = pi2*xi*(nf-1.d0)*A
        fg0k(2)%less%w(i) = pi2*xi*nf*B!*delta/abs(wr(i))
        fg0k(2)%gtr%w(i)  = pi2*xi*(nf-1.d0)*B!*delta/abs(wr(i))
     enddo

     calG11%less%w = fg0k(1)%less%w    
     calG11%gtr%w = fg0k(1)%gtr%w    

     forall(i=1:L)calG22%less%w(i) = -(fg0k(1)%gtr%w(L+1-i))
     forall(i=1:L)calG22%gtr%w(i) = -(fg0k(1)%less%w(L+1-i))

     calF12%less%w   = -conjg(fg0k(2)%less%w)
     calF12%gtr%w   =  fg0k(2)%gtr%w

     calF21%less%w   = fg0k(2)%less%w
     calF21%gtr%w   = -conjg(fg0k(2)%gtr%w)


     call fftgf_rw2rt(calG11%less%w,calG11%less%t,Lm) ; calG11%less%t=fmesh/pi2*calG11%less%t
     call fftgf_rw2rt(calG11%gtr%w,calG11%gtr%t,Lm)   ; calG11%gtr%t =fmesh/pi2*calG11%gtr%t
     call fftgf_rw2rt(calG22%less%w,calG22%less%t,Lm) ; calG22%less%t=fmesh/pi2*calG22%less%t
     call fftgf_rw2rt(calG22%gtr%w,calG22%gtr%t,Lm)   ; calG22%gtr%t =fmesh/pi2*calG22%gtr%t
     call fftgf_rw2rt(calF12%less%w,calF12%less%t,Lm) ; calF12%less%t=fmesh/pi2*calF12%less%t
     call fftgf_rw2rt(calF12%gtr%w,calF12%gtr%t,Lm)   ; calF12%gtr%t =fmesh/pi2*calF12%gtr%t
     call fftgf_rw2rt(calF21%less%w,calF21%less%t,Lm) ; calF21%less%t=fmesh/pi2*calF21%less%t
     call fftgf_rw2rt(calF21%gtr%w,calF21%gtr%t,Lm)   ; calF21%gtr%t =fmesh/pi2*calF21%gtr%t

     do i=-Lm,Lm 
        sk(1)%less%t(i) = U**2*(calG11%less%t(i)*calG22%less%t(i) - calF12%less%t(i)*calF21%less%t(i))*calG22%gtr%t(-i)
        sk(1)%gtr%t(i) =  U**2*(calG11%gtr%t(i)*calG22%gtr%t(i) - calF12%gtr%t(i)*calF21%gtr%t(i))*calG22%less%t(-i)
        !
        sk(2)%less%t(i) = U**2*(calF12%less%t(i)*calF21%less%t(i) - calG11%less%t(i)*calG22%less%t(i))*calF12%gtr%t(-i)
        sk(2)%gtr%t(i) =  U**2*(calF12%gtr%t(i)*calF21%gtr%t(i)  - calG11%gtr%t(i)*calG22%gtr%t(i))*calF12%less%t(-i)
        !
        sk(1)%ret%t(i) =heaviside(t(i))*(sk(1)%gtr%t(i)-sk(1)%less%t(i))
        sk(2)%ret%t(i) =heaviside(t(i))*(sk(2)%gtr%t(i)-sk(2)%less%t(i))
     enddo
     if(heaviside(0.d0)==1.d0)sk(1)%ret%t(0)=sk(1)%ret%t(0)/2.d0 
     if(heaviside(0.d0)==1.d0)sk(2)%ret%t(0)=sk(2)%ret%t(0)/2.d0
     call fftgf_rt2rw(sk(1)%ret%t,sk(1)%ret%w,Lm) ;      sk(1)%ret%w=dt*sk(1)%ret%w
     call fftgf_rt2rw(sk(2)%ret%t,sk(2)%ret%w,Lm) ;      sk(2)%ret%w=dt*sk(2)%ret%w


     sigma(1,:) = sk(1)%ret%w
     sigma(2,:) = -delta + sk(2)%ret%w

     write(*,"(2f14.9)",advance="no")2.d0*n,delta
     !sold=sigma
     ! sigma =  solve_ipt_sc_sopt(calG(1:2,:),wr,delta,L)
     ! sigma = weight*sigma + (1.d0-weight)*sold
     !converged = check_convergence(sigma(1,:)+sigma(2,:),eps=eps_error,N1=Nsuccess,N2=Nloop)

     ! call splot("DOS.ipt",wr,-dimag(fg(1,:))/pi,append=.true.)
     ! call splot("G_realw.ipt",wr,fg(1,:),append=.true.)
     ! call splot("F_realw.ipt",wr,fg(2,:),append=.true.)
     ! call splot("Sigma_realw.ipt",wr,sigma(1,:),append=.true.)
     ! call splot("Self_realw.ipt",wr,sigma(2,:),append=.true.)
     ! call splot("calG0_realw.ipt",wr,calG(1,:),append=.true.)
     ! call splot("calF0_realw.ipt",wr,calG(2,:),append=.true.)
     call splot("observables.ipt",vbias,delta,xmu,u,n,beta,dble(iloop),append=.true.)
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=eps_error,N1=Nsuccess,N2=Nloop)
  enddo

  call splot("observables.last",vbias,delta,xmu,u,n,beta,dble(iloop),append=printf)
  call splot("DOS.last",wr,-dimag(fg(1,:))/pi,append=printf)
  call splot("G_realw.last",wr,fg(1,:),append=printf)
  call splot("F_realw.last",wr,fg(2,:),append=printf)
  call splot("Sigma_realw.last",wr,sigma(1,:),append=printf)
  call splot("Self_realw.last",wr,sigma(2,:),append=printf)
  call splot("calG0_realw.last",wr,calG(1,:),append=printf)
  call splot("calF0_realw.last",wr,calG(2,:),append=printf)

contains 

  subroutine get_initial_sigma()
    logical                :: check1,check2,check
    inquire(file="Sigma_realw.restart",exist=check1)
    inquire(file="Self_realw.restart",exist=check2)
    check=check1*check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_realw.restart",wr,sigma(1,:))
       call sread("Self_realw.restart",wr,sigma(2,:))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 ;  delta=deltasc 
       sigma(2,:)=-delta ; sigma(1,:)=zero ; sold=sigma
    endif
  end subroutine get_initial_sigma

end program hmipt
