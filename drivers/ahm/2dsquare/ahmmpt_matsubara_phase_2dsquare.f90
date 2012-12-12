!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program ahmmpt_matsubara_phase_2dsquare
  USE DMFT_IPT
  USE IOTOOLS
  USE RANDOM
  USE CHRONOBAR
  USE SQUARE_LATTICE
  implicit none
  integer                         :: i,ik,esp,Lk,M
  logical                         :: converged,check1,check2,check
  complex(8)                      :: zeta,cdet
  complex(8)                      :: delta,delta0
  real(8)                         :: n,n0,w
  real(8)                         :: sc_phase,sc_mod
  real(8)                         :: sc0_phase,sc0_mod

  !
  complex(8),allocatable          :: fg(:,:),fg0(:,:),sigma(:,:),calG(:,:)
  real(8),allocatable             :: fgt(:)
  complex(8),allocatable          :: fft(:)
  complex(8),allocatable          :: det(:),sold(:,:)
  !
  real(8),allocatable             :: wt(:),epsik(:),wm(:),tau(:)


  include "revision.inc"
  call version(revision)
  call read_input("inputIPT.in")


  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)

  !
  allocate(fg(2,L),fg0(2,L),calG(2,L))
  allocate(sigma(2,L),sold(2,L))
  allocate(det(L))
  allocate(fgt(0:L),fft(0:L))


  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)
  call get_initial_sigma

  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="yes")"DMFT-loop",iloop

     call start_timer
     fg=zero
     do i=1,L
        zeta =  xi*wm(i) + xmu - sigma(1,i)
        do ik=1,Lk
           cdet = abs(zeta-epsik(ik))**2 + abs(sigma(2,i))**2
           fg(1,i)=fg(1,i) + wt(ik)*conjg(zeta-epsik(ik))/cdet
           fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/cdet !inverse is transpose!!!
        enddo
     enddo

     call fftgf_iw2tau(fg(1,:),fgt(0:),beta)
     call fftff_iw2tau(fg(2,:),fft(0:),beta)
     n=-fgt(L) ; delta= -u*fft(L)

     !calcola calG0^-1, calF0^-1 (WFs)
     det      =  abs(fg(1,:))**2     + abs(fg(2,:))**2
     fg0(1,:) =  conjg(fg(1,:))/det  + sigma(1,:) + u*(n-0.5d0)
     fg0(2,:) =  fg(2,:)/det         + sigma(2,:) +  delta

     det       =  abs(fg0(1,:))**2   + abs(fg0(2,:))**2
     calG(1,:) =  conjg(fg0(1,:))/det
     calG(2,:) =  fg0(2,:)/det

     call fftgf_iw2tau(calG(1,:),fgt(0:),beta)
     call fftff_iw2tau(calG(2,:),fft(0:),beta)
     n0=-fgt(L) ; delta0= -u*fft(L)

     sc_mod=abs(delta)
     sc0_mod=abs(delta0)
     sc_phase=atan2(dimag(delta),real(delta,8))
     sc0_phase=atan2(dimag(delta0),real(delta,8))

     write(*,"(6(f16.12))",advance="no")n,n0,sc_mod,sc_phase,sc0_mod,sc0_phase

     sigma =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
     sigma =  weight*sigma + (1.d0-weight)*sold;sold=sigma

     !this is an idea of Massimo, check error as a single array.
     converged = check_convergence_scalar(sc_mod,eps=eps_error,N1=Nsuccess,N2=nloop)

     if(nread/=0.d0)call search_mu(converged)

     call splot("nVSiloop.ipt",iloop,n,n0,append=TT)
     call splot("deltaVSiloop.ipt",iloop,sc_mod,sc0_mod,append=TT)
     call splot("phaseVSiloop.ipt",iloop,sc_phase,sc0_phase,append=TT)

     ! This stays here for the time being, until the code is bomb proof
     call splot("G_iw.ipt",wm,fg(1,:),append=printf)
     call splot("F_iw.ipt",wm,fg(2,:),append=printf)
     call splot("G_tau.ipt",tau,fgt,append=printf)
     call splot("F_tau.ipt",tau,fft,append=printf)
     call splot("calG_iw.ipt",wm,calG(1,:),append=printf)
     call splot("calF_iw.ipt",wm,calG(2,:),append=printf)
     call splot("calG_tau.ipt",tau,fgt,append=printf)
     call splot("calF_tau.ipt",tau,fft,append=printf)
     call splot("Sigma_iw.ipt",wm,sigma(1,:),append=printf)
     call splot("Self_iw.ipt",wm,sigma(2,:),append=printf)
     call splot("observables.ipt",beta,xmu,u,n,n0,sc_mod,sc0_mod,sc_phase,sc0_phase,append=printf)
     call stop_timer
  enddo

  call close_file("nVSiloop.ipt")
  call close_file("deltaVSiloop.ipt")
  call close_file("phaseVSiloop.ipt")

  call splot("G_iw.last",wm,fg(1,:),append=FF)
  call splot("F_iw.last",wm,fg(2,:),append=FF)
  call splot("G_tau.last",tau,fgt,append=FF)
  call splot("F_tau.last",tau,fft,append=FF)
  call splot("calG_iw.last",wm,calG(1,:),append=FF)
  call splot("calF_iw.last",wm,calG(2,:),append=FF)
  call splot("Sigma_iw.last",wm,sigma(1,:),append=FF)
  call splot("Self_iw.last",wm,sigma(2,:),append=FF)
  call splot("observables.last",beta,xmu,u,n,n0,sc_mod,sc0_mod,sc_phase,sc0_phase,append=FF)

  !call get_sc_internal_energy

contains  

  subroutine search_mu(convergence)
    integer, save         :: nindex
    integer               :: nindex1
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
       ndelta=ndelta1/2.d0    !decreasing the step
    else
       ndelta=ndelta1
    endif
    xmu=xmu+real(nindex,8)*ndelta

    write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",nread,"| shift=",nindex*ndelta,"| mu=",xmu

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
       call init_random_number
       call random_number(sc_phase)
       sc_phase=sc_phase*pi2
       n=0.5d0
       delta=abs(deltasc)*exp(xi*sc_phase)
       sigma(2,:)=-delta ; sigma(1,:)=zero
       sold=sigma
    endif
  end subroutine get_initial_sigma

end program ahmmpt_matsubara_phase_2dsquare
