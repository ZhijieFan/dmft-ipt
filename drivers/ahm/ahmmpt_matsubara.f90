!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmmpt_matsubara
  USE DMFT_IPT
  USE IOTOOLS
  implicit none
  integer                         :: i,ik,esp,Lk
  logical                         :: converged,check1,check2,check
  complex(8)                      :: zeta,cdet
  real(8)                         :: n,delta,n0,delta0
  !
  complex(8),allocatable          :: fg(:,:),fg0(:,:),sigma(:,:),calG(:,:)
  real(8),allocatable             :: fgt(:,:),fg0t(:,:)
  complex(8),allocatable          :: det(:),sold(:,:),sconvergence(:)
  !
  real(8),allocatable             :: wt(:),epsik(:),wm(:),tau(:)


  include "revision.inc"
  call version(revision)
  call read_input("inputIPT.in")

  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)
  write(*,"(A,I9,A)")"Using ",L," frequencies"
  !
  allocate(fg(2,L),fgt(2,0:L))
  allocate(fg0(2,L),fg0t(2,0:L))
  allocate(sigma(2,L),calG(2,L))
  allocate(det(L),Sold(2,L))
  allocate(sconvergence(2*L))


  D=2.d0*ts; Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D)

  call get_initial_sigma


  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     fg=zero
     do i=1,L
        zeta =  xi*wm(i) + xmu - sigma(1,i)
        do ik=1,Lk
           cdet = abs(zeta-epsik(ik))**2 + (sigma(2,i))**2 
           fg(1,i)=fg(1,i) + wt(ik)*(conjg(zeta)-epsik(ik))/cdet

           fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/cdet
        enddo
     enddo
     call fftgf_iw2tau(fg(1,:),fgt(1,0:L),beta)
     call fftgf_iw2tau(fg(2,:),fgt(2,0:L),beta,notail=.true.)
     n=-fgt(1,L) ; delta= -u*fgt(2,L)

     !calcola calG0^-1, calF0^-1 (WFs)
     det      =  abs(fg(1,:))**2 + fg(2,:)**2
     fg0(1,:) =  conjg(fg(1,:))/det + sigma(1,:) + u*(n-0.5d0)
     fg0(2,:) =  fg(2,:)/det        + sigma(2,:) +  delta

     det       =  abs(fg0(1,:))**2 + fg0(2,:)**2
     calG(1,:) =  conjg(fg0(1,:))/det
     calG(2,:) =  fg0(2,:)/det

     call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
     call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.) !; fg0t(2,:)=-fg0t(2,:)
     n0=-fg0t(1,L) ; delta0= -u*fg0t(2,L)
     write(*,"(4(f16.12))",advance="no")n,n0,delta,delta0

     sigma     =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
     sigma     =  weight*sigma + (1.d0-weight)*sold;sold=sigma

     !this is an idea of Massimo, instead of checking error on the sum, do it 
     !as a single array. 
     sconvergence(1:L)=sigma(1,1:L) ; sconvergence(L+1:2*L)=sigma(2,1:L)
     converged = check_convergence(sconvergence(:),eps=eps_error,N1=Nsuccess,N2=nloop)

     if(nread/=0.d0)call search_mu(converged)

     call splot("nVSiloop.ipt",iloop,n,append=TT)
     call splot("deltaVSiloop.ipt",iloop,delta,append=TT)

     !This stays here for the time being, until the code is bomb proof
     call splot("G_iw.ipt",wm,fg(1,:),append=printf)
     call splot("F_iw.ipt",wm,fg(2,:),append=printf)
     call splot("G_tau.ipt",tau,fgt(1,:),append=printf)
     call splot("F_tau.ipt",tau,fgt(2,:),append=printf)
     call splot("calG_iw.ipt",wm,calG(1,:),append=printf)
     call splot("calF_iw.ipt",wm,calG(2,:),append=printf)
     call splot("calG_tau.ipt",tau,fg0t(1,:),append=printf)
     call splot("calF_tau.ipt",tau,fg0t(2,:),append=printf)
     call splot("Sigma_iw.ipt",wm,sigma(1,:),append=printf)
     call splot("Self_iw.ipt",wm,sigma(2,:),append=printf)
     call splot("observables.ipt",xmu,u,n,n0,delta,delta0,beta,dble(iloop),append=printf)
  enddo

  call close_file("nVSiloop.ipt")
  call close_file("deltaVSiloop.ipt")

  call splot("G_iw.last",wm,fg(1,:),append=FF)
  call splot("F_iw.last",wm,fg(2,:),append=FF)
  call splot("G_tau.last",tau,fgt(1,:),append=FF)
  call splot("F_tau.last",tau,fgt(2,:),append=FF)
  call splot("calG_iw.last",wm,calG(1,:),append=FF)
  call splot("calF_iw.last",wm,calG(2,:),append=FF)
  call splot("calG_tau.last",tau,fg0t(1,:),append=FF)
  call splot("calF_tau.last",tau,fg0t(2,:),append=FF)
  call splot("Sigma_iw.last",wm,sigma(1,:),append=FF)
  call splot("Self_iw.last",wm,sigma(2,:),append=FF)
  call splot("observables.last",xmu,u,n,n0,delta,delta0,beta,dble(iloop),append=FF)

  call get_sc_internal_energy

contains

  include "internal_energy_ahm_matsubara.f90"

  subroutine search_mu(convergence)
    integer, save         ::nindex
    integer               ::nindex1
    real(8)               :: naverage,ndelta1
    logical,intent(inout) :: convergence
    naverage=n
    nindex1=nindex
    ndelta1=ndelta
    if((naverage >= nread+nerror))then
       nindex=-1
    elseif(naverage <= nread-nerror)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0)then !avoid loop forth and back
       ndelta=real(ndelta1/2.d0,8) !decreasing the step
       xmu=xmu+real(nindex,8)*ndelta
    else
       ndelta=ndelta1
       xmu=xmu+real(nindex,8)*ndelta
    endif
    write(*,"(A,2f15.12,A,f15.12,A,I3,f15.12)")"mu,n=",xmu,naverage,"/",nread,"| ",nindex,ndelta
    if(abs(naverage-nread)>nerror)convergence=.false.
    call splot("muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
  end subroutine search_mu


  subroutine get_initial_sigma()
    inquire(file="Sigma_iw.last",exist=check1)
    inquire(file="Self_iw.last",exist=check2)
    check=check1.AND.check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_iw.last",wm,sigma(1,:))
       call sread("Self_iw.last",wm,sigma(2,:))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 ; delta=deltasc
       sigma(2,:)=-delta ; sigma(1,:)=zero
       sold=sigma
    endif
  end subroutine get_initial_sigma

end program hmmpt_matsubara
