!########################################################
!     Program  : AHMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the attractive Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
module COMMON
  USE BROYDEN
  implicit none
  !Put here vars in common with the BROYDN function
  real(8) :: xmu0,n,n0,delta,delta0
  complex(8)             :: zeta1,zeta2,det
  complex(8),allocatable :: fg0(:,:),calG(:,:)
end module COMMON

program ahmmpt
  USE DMFT_IPT
  USE COMMON
  implicit none
  real(8)    :: x(1)
  logical    :: check
  integer                :: ik,Lk
  logical                :: converged
  complex(8),allocatable :: sigma(:,:),fg(:,:),sold(:,:)
  real(8),allocatable    :: wt(:),epsik(:)

  call read_input("inputIPT.in")
  allocate(fg(2,-L:L),sigma(2,-L:L),fg0(2,-L:L),calG(2,-L:L))
  allocate(sold(2,-L:L))

  Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(D,wt,epsik,Lk)


  n=0.5d0 ;  delta=deltasc ;xmu0=xmu
  sigma=zero ; sigma(2,:)=-delta 
  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="yes")"DMFT-loop",iloop

     fg=zero
     do i=-L,L
        zeta1 = cmplx(wr(i),eps,8) + xmu - sigma(1,i)
        zeta2 = cmplx(wr(i),eps,8) - xmu + sigma(1,-i)
        do ik=1,Lk
           det     = (zeta1-epsik(ik))*(zeta2+epsik(ik)) - sigma(2,i)**2
           fg(1,i) =fg(1,i) + wt(ik)*(zeta2+epsik(ik))/det
           fg(2,i) =fg(2,i) - wt(ik)*sigma(2,i)/det
        enddo
     enddo
     where(aimag(fg(1,:))>0.d0)fg(1,:)=real(fg(1,:),8)
     n    = -sum(aimag(fg(1,:))*fermi(wr,beta))*fmesh/pi
     delta=  u*sum(aimag(fg(2,:))*fermi(wr,beta))*fmesh/pi
     call splot("gloc.ipt",wr,fg(1,:))
     call splot("floc.ipt",wr,fg(2,:))

     !Get the Weiss Field
     !G0^-1 = Gloc^-1 + Sigma
     do i=-L,L
        det     = fg(1,i)*fg(1,-i) - fg(2,i)**2
        fg0(1,i)= fg(1,-i)/det + sigma(1,i)! + u*(n-0.5d0)
        fg0(2,i)=-fg(2,i)/det  + sigma(2,i)! - delta
     end do
     call splot("g0.ipt",wr,fg0(1,:))
     call splot("f0.ipt",wr,fg0(2,:))

     x(1)=xmu
     call broydn(x,check)
     xmu0=x(1)

     !write(*,"(4(f16.12))",advance="no"),n,n0,delta,delta0
     sigma =  solve_mpt_sc_sopt(calG,fg,n,n0,delta,delta0)
     where(aimag(sigma(1,:))>0.d0)sigma(1,:)=real(sigma(1,:),8)
     sigma=weigth*sigma + (1.d0-weigth)*sold ; sold=sigma
     converged = check_convergence(sigma(1,:)+sigma(2,:))
  enddo

end program ahmmpt

function funcv(x)
  USE COMMON
  USE DMFT_IPT
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0=x(1)

  !Hartree corrected WF is: xmu=xmu0
  !\tilde{\calG0}^-1 = G0^-1 +xmu -xmu0 - \Sigma_HFB
  do i=-L,L
     det     =   fg0(1,i)*fg0(1,-i) - fg0(2,i)**2
     calG(1,i)=  fg0(1,-i)/det + xmu-xmu0 + U*(n-0.5d0)
     calG(2,i)= -fg0(2,i)/det  - delta
  end do
  !where(aimag(calG(1,:))>0.d0)calG(1,:)=real(calG(1,:),8)
  n0    =  -sum(aimag(calG(1,:))*fermi(wr,beta))*fmesh/pi
  delta0= u*sum(aimag(calG(2,:))*fermi(wr,beta))*fmesh/pi
  call splot("calG.ipt",wr,calG(1,:))
  call splot("calF.ipt",wr,calG(2,:))

  funcv(1)=n-n0
  write(*,"(3(f13.9))")n,n0,xmu0
end function funcv
