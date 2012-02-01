!########################################################
!     Program  : AHMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the attractive Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program ahmmpt
  USE DMFT_IPT
  USE IOTOOLS
  implicit none
  integer                :: i,ik,Lk
  logical                :: converged
  complex(8)             :: zeta1,zeta2,det
  real(8)                :: n,delta,n0,delta0
  !
  complex(8),allocatable :: sigma(:,:),fg(:,:)
  complex(8),allocatable :: wf0(:,:),calG(:,:)
  complex(8),allocatable :: sold(:,:),zeta(:)
  !
  real(8),allocatable    :: wt(:),epsik(:),wr(:)

  call read_input("inputIPT.in")
  allocate(fg(2,L),sigma(2,L))
  allocate(wf0(2,L),calG(2,L))
  allocate(sold(2,L),zeta(L))
  !
  allocate(wr(L))
  wr = linspace(-wmax,wmax,L,mesh=fmesh)
  !
  D=2.d0*ts
  Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D_=D,eps_=1.d-5)


  n=0.5d0 ;  delta=deltasc 
  sigma(2,:)=-delta ; sigma(1,:)=zero ;sold=sigma
  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     fg=zero
     zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
     do i=1,L
        zeta1 = zeta(i)
        zeta2 = conjg(zeta(L-i+1))
        do ik=1,Lk
           det     = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + sigma(2,i)*sigma(2,L-i+1)
           fg(1,i) =fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det
           fg(2,i) =fg(2,i) - wt(ik)*sigma(2,L-i+1)/det
        enddo
     enddo
     n    = -sum(aimag(fg(1,:))*fermi(wr,beta))*fmesh/pi
     delta= -(u*sum(aimag(fg(2,:))*fermi(wr,beta))*fmesh/pi)

     !Hartree corrected WF is: xmu=xmu0
     !\tilde{\calG0} = [G^-1 + Sigma - \Sigma_HFB]^-1
     do i=1,L
        det     = fg(1,i)*conjg(fg(1,L-i+1)) - fg(2,i)*fg(2,L-i+1)
        wf0(1,i)= conjg(fg(1,L-i+1))/det  + sigma(1,i)   - u*(n-0.5d0)
        wf0(2,i)= fg(2,L-i+1)/det         + sigma(2,i)   + delta
     end do
     do i=1,L
        det      =  wf0(1,i)*conjg(wf0(1,L-i+1)) - wf0(2,i)*wf0(2,L-i+1)
        calG(1,i)=  conjg(wf0(1,L-i+1))/det
        calG(2,i)=  wf0(2,L-i+1)/det
     end do
     n0    =  -sum(aimag(calG(1,:))*fermi(wr,beta))*fmesh/pi
     delta0=  -(u*sum(aimag(calG(2,:))*fermi(wr,beta))*fmesh/pi)
     write(*,"(4(f16.12))",advance="no"),n,n0,delta,delta0

     sigma =  solve_mpt_sc_sopt(calG,wr,n,n0,delta,delta0)
     sigma=weigth*sigma + (1.d0-weigth)*sold ; sold=sigma
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps_error,Nsuccess,Nloop)
     call splot("nVSiloop.ipt",iloop,n,append=TT)
     call splot("deltaVSiloop.ipt",iloop,delta,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  call close_file("deltaVSiloop.ipt")
  call splot("DOS.ipt",wr,-aimag(fg(1,:))/pi,append=printf)
  call splot("G_realw.ipt",wr,fg(1,:),append=printf)
  call splot("F_realw.ipt",wr,fg(2,:),append=printf)
  call splot("Sigma_realw.ipt",wr,sigma(1,:),append=printf)
  call splot("Self_realw.ipt",wr,sigma(2,:),append=printf)
  call splot("calG0_realw.ipt",wr,calG(1,:),append=printf)
  call splot("calF0_realw.ipt",wr,calG(2,:),append=printf)
  call splot("n.delta_realw.ipt",n,delta,append=printf)

end program ahmmpt
