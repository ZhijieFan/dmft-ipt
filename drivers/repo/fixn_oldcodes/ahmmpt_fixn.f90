!########################################################
!     Program  : AHMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the attractive Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
module COMMON
  USE DMFT_IPT
  USE BROYDEN
  USE LATTICE
  implicit none
  complex(8),allocatable :: sigma(:,:),fg(:,:),fg0(:,:),calG(:,:),sold(:,:)
  real(8)                :: n,delta,n0,delta0,ntarget
end module COMMON



program ahmmpt_fixn
  USE COMMON
  implicit none
  real(8)    :: x(1)
  logical    :: check
  logical                :: converged
  integer                :: ik
  complex(8)             :: zeta1,zeta2,det
  call read_input("inputIPT.in")
  allocate(fg(2,-L:L),sigma(2,-L:L),fg0(2,-L:L),calG(2,-L:L))
  allocate(sold(2,-L:L))

  call  build_BetheLattice(Nx,D,Nk=Lk)
  allocate(epsik(Lk)) ; call get_epsik_bethe(epsik,D)

  ntarget=xmu ; xmu=0.d0
  n=0.5d0 ;  delta=deltasc
  sigma=zero ; sigma(2,:)=-delta ; sold=sigma

  iloop=0 ; converged=.false. ; sold=sigma
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="yes")"DMFT-loop",iloop
     x(1)=xmu
     call broydn(x,check)
     do i=-L,L
        det      = fg(1,i)*fg(1,-i) - fg(2,i)*(fg(2,i))
        calG(1,i)= fg(1,-i)/det + sigma(1,i) - u*(n-0.5d0)
        calG(2,i)= fg(2,i)/det  + sigma(2,i) - delta
     end do
     do i=-L,L
        det     =  calG(1,i)*calG(1,-i) - calG(2,i)*calG(2,i)
        fg0(1,i)=  calG(1,-i)/det
        fg0(2,i)=  calG(2,i)/det
     end do
     n0    =  sum(aimag(fg0(1,:))*fermi(wr,beta))/sum(dimag(fg0(1,:))) !*fmesh/pi
     delta0= -u*sum(aimag(fg0(2,:))*fermi(wr,beta))*fmesh/pi
     delta = -u*sum(aimag(fg(2,:))*fermi(wr,beta))*fmesh/pi
     sigma =  solve_mpt_sc_sopt(fg0,fg,n,n0,delta,delta0)
     sigma = weigth*sigma + (1.d0-weigth)*sold ; sold=sigma
     converged = check_convergence(sigma(1,:)+sigma(2,:))
  end do
end program ahmmpt_fixn



function funcv(x)
  USE COMMON
  implicit none
  real(8),intent(in) :: x(:)
  real(8)            :: funcv(size(x))
  integer                :: ik
  complex(8)             :: zeta1,zeta2,det
  xmu=x(1)
  fg=zero
  do i=-L,L
     zeta1 = cmplx(wr(i),eps,8) + xmu - sigma(1,i)
     zeta2 = cmplx(wr(i),eps,8) - xmu + sigma(1,-i)
     do ik=1,Lk
        det     = (zeta1-epsik(ik))*(zeta2+epsik(ik)) - sigma(2,i)*sigma(2,i)
        fg(1,i) =fg(1,i) + wt(ik)*(zeta2+epsik(ik))/det
        fg(2,i) =fg(2,i) + wt(ik)*sigma(2,i)/det
     enddo
  enddo
  n    =  sum(aimag(fg(1,:))*fermi(wr,beta))/sum(dimag(fg(1,:))) !*fmesh/pi
  if(xmu>3.d0 .OR. xmu<-3.d0)xmu=0.d0
  funcv(1)=2.d0*n-ntarget
  write(*,"(3(f13.9))")2.d0*n,ntarget,xmu
end function funcv
