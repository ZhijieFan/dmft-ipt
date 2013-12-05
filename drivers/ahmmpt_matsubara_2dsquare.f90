!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmmpt_2dsquare_matsubara
  USE DMFT_IPT
  USE SCIFOR_VERSION
  !USE SQUARE_LATTICE
  USE FUNCTIONS
  USE INTEGRATE
  USE IOTOOLS
  implicit none
  integer                         :: i,ik,esp,Lk,iloop
  logical                         :: converged,check1,check2,check
  complex(8)                      :: zeta,cdet
  real(8)                         :: n,delta,n0,delta0,dtau,wband,de
  !
  complex(8),allocatable          :: fg(:,:),fg0(:,:),sigma(:,:),calG(:,:)
  real(8),allocatable             :: fgt(:,:),fg0t(:,:)
  complex(8),allocatable          :: det(:),sold(:,:),sconvergence(:)
  !
  real(8),allocatable             :: wt(:),epsik(:),wm(:),tau(:)


  include "revision.inc"
  call version(revision)

  call read_input("inputIPT.in")

  !allocate grids
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)


  !allocate functions
  allocate(fg(2,L),fgt(2,0:L))
  allocate(fg0(2,L),fg0t(2,0:L))
  allocate(sigma(2,L),calG(2,L))
  allocate(det(L),Sold(2,L))
  allocate(sconvergence(2*L))


  !build square lattice structure:
  !Lk   = square_lattice_dimension(Nx)
  !allocate(wt(Lk),epsik(Lk))
  !wt   = square_lattice_structure(Lk,Nx)
  !epsik= square_lattice_dispersion_array(Lk,ts)
  Lk=Nx
  allocate(epsik(Lk),wt(Lk))
  print*,"Using ",Lk," points for the e-integral"
  wband=4.d0*ts+0.5d0
  epsik = linspace(-wband,wband,Lk,mesh=de)
  do i=1,Lk
     wt(i)=dens_2dsquare(epsik(i),ts)
  enddo
  wt=wt/trapz(de,wt)
  call splot("DOS2d.ipt",epsik,wt)
  wt = wt*de

  iloop=0 ; converged=.false.

  call get_initial_sigma

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

     if(iloop>1)calG =  weight*calG + (1.d0-weight)*sold
     sold=calG

     call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
     call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.) !; fg0t(2,:)=-fg0t(2,:)
     n0=-fg0t(1,L) ; delta0= -u*fg0t(2,L)
     write(*,"(4(f16.12))",advance="no")n,n0,delta,delta0

     sigma =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)

     !this is an idea of Massimo, instead of checking error on the sum, use single array. 
     sconvergence(1:L)=fg(1,1:L) ; sconvergence(L+1:2*L)=fg(2,1:L)
     converged = check_convergence(sconvergence(:),eps=eps_error,N1=Nsuccess,N2=nloop)

     if(nread/=0.d0)call search_mu(converged)
     if(iloop==nloop)converged=.true.

     call splot("nVSiloop.ipt",iloop,n,append=.true.)
     call splot("deltaVSiloop.ipt",iloop,delta,append=.true.)

     !This stays here for the time being, until the code is bomb proof
     if(printf)then
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
        call splot("observables.ipt",iloop,beta,xmu,u,n,n0,delta,delta0,append=printf)
     endif
  enddo

  call close_file("nVSiloop.ipt")
  call close_file("deltaVSiloop.ipt")

  call splot("G_iw.last",wm,fg(1,:),append=.false.)
  call splot("F_iw.last",wm,fg(2,:),append=.false.)
  call splot("G_tau.last",tau,fgt(1,:),append=.false.)
  call splot("F_tau.last",tau,fgt(2,:),append=.false.)
  call splot("calG_iw.last",wm,calG(1,:),append=.false.)
  call splot("calF_iw.last",wm,calG(2,:),append=.false.)
  call splot("calG_tau.last",tau,fg0t(1,:),append=.false.)
  call splot("calF_tau.last",tau,fg0t(2,:),append=.false.)
  call splot("Sigma_iw.last",wm,sigma(1,:),append=.false.)
  call splot("Self_iw.last",wm,sigma(2,:),append=.false.)
  call splot("observables.last",iloop,beta,xmu,u,n,n0,delta,delta0,append=.false.)

contains


  subroutine search_mu(convergence)
    integer, save         ::nindex,be_far=0
    integer               ::nindex1
    real(8)               :: naverage,ndelta1
    logical,intent(inout) :: convergence
    naverage=2.d0*n
    nindex1=nindex
    ndelta1=ndelta
    if((naverage >= nread+nerror))then
       nindex=-1
    elseif(naverage <= nread-nerror)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
       ndelta=ndelta1/2.d0
    else
       ndelta=ndelta1
       ! if(abs(nindex+nindex1)==2)then
       !    be_far=be_far+1
       ! else
       !    be_far=0
       ! endif
       ! if(be_far>6)then
       !    ndelta=2.d0*ndelta
       !    be_far=0
       ! endif
    endif
    xmu=xmu+real(nindex,8)*ndelta
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
    endif
  end subroutine get_initial_sigma

end program
