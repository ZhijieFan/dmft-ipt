!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt_matsubara
  USE DMFT_IPT
  USE IOTOOLS
  implicit none

  integer                :: i,ik,Lk,iloop
  logical                :: converged
  complex(8)             :: zeta,zeta1,zeta2,cdet,x1,x2,zsqrt
  real(8)                :: n,delta
  !
  complex(8),allocatable :: fg(:,:),fg0(:,:),sigma(:,:),calG(:,:),det(:)
  real(8),allocatable    :: fgt(:,:)
  !
  real(8),allocatable    :: wt(:),epsik(:),wm(:)

  call read_input("inputIPT.in")

  allocate(wm(L))
  wm  = pi/beta*real(2*arange(1,L)-1,8)

  allocate(fg(2,L),fgt(2,0:L))
  allocate(fg0(2,L),sigma(2,L),calG(2,L),det(L))

  D=2.d0*ts
  Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D)

  call get_initial_sigma

  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     fg=zero
     do i=1,L
        zeta =  xi*wm(i) + xmu - sigma(1,i)
        fg(:,i)=zero
        do ik=1,Lk
           cdet = abs(zeta-epsik(ik))**2 + (sigma(2,i))**2
           fg(1,i)=fg(1,i) + wt(ik)*(conjg(zeta)-epsik(ik))/cdet
           fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/cdet
        enddo
     enddo
     call fftgf_iw2tau(fg(1,:),fgt(1,0:L),beta)
     call fftgf_iw2tau(fg(2,:),fgt(2,0:L),beta,notail=.true.)
     n=-fgt(1,L) ; delta= -u*fgt(2,L)


     det       = abs(fg(1,:))**2 + (fg(2,:))**2
     fg0(1,:) =  conjg(fg(1,:))/det + sigma(1,:)
     fg0(2,:) =  fg(2,:)/det        + sigma(2,:) + delta

     det       =  abs(fg0(1,:))**2 + (fg0(2,:))**2
     calG(1,:) =  conjg(fg0(1,:))/det
     calG(2,:) =  fg0(2,:)/det


     write(*,"(4(f16.12))",advance="no")delta
     sigma =  solve_ipt_sc_matsubara(calG,delta)
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=eps_error,N1=Nsuccess,N2=Nloop)
  enddo

  call splot("Sigma_iw.last",wm,sigma(1,:),append=printf)
  call splot("Self_iw.last",wm,sigma(2,:),append=printf)
  call splot("G_iw.last",wm,fg(1,:),append=printf)
  call splot("F_iw.last",wm,fg(2,:),append=printf)
  call splot("calG_iw.last",wm,calG(1,:),append=printf)
  call splot("calF_iw.last",wm,calG(2,:),append=printf)
  call splot("observables.last",u,beta,n,delta,append=printf)

  call get_sc_internal_energy

contains

  include "internal_energy_ahm_matsubara.f90"

  subroutine get_initial_sigma()
    logical :: check1,check2,check
    inquire(file="Sigma_iw.restart",exist=check1)
    inquire(file="Self_iw.restart",exist=check2)
    check=check1.AND.check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_iw.restart",wm,sigma(1,:))
       call sread("Self_iw.restart",wm,sigma(2,:))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 ; delta=deltasc
       sigma(2,:)=-delta ; sigma(1,:)=zero
    endif
  end subroutine get_initial_sigma

end program hmipt_matsubara
