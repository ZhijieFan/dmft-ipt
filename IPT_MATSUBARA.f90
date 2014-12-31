!###############################################################
!     PURPOSE  : SOLVE DMFT-IPT REPULSIVE IN MATSUBARA FREQUENCY
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_MATSUBARA
  USE IPT_VARS_GLOBAL
  USE FFTGF
  implicit none
  private

  interface ipt_solve_matsubara
     module procedure solve_ipt_matsubara_1b!,solve_ipt_matsubara_mb
  end interface ipt_solve_matsubara

  interface mpt_solve_matsubara
     module procedure solve_ipt_matsubara_1b!,solve_ipt_matsubara_mb
  end interface mpt_solve_matsubara

  interface ipt_solve_matsubara_sc
     module procedure solve_ipt_sc_matsubara_r!,solve_ipt_sc_matsubara_c
  end interface ipt_solve_matsubara_sc

  interface mpt_solve_matsubara_sc
     module procedure solve_mpt_sc_matsubara_r!,solve_mpt_sc_matsubara_c
  end interface mpt_solve_matsubara_sc

  !half-filling solver:
  public :: ipt_solve_matsubara
  public :: ipt_solve_matsubara_sc
  !away from h-f solver:
  public :: mpt_solve_matsubara
  public :: mpt_solve_matsubara_sc
  !
  !public :: solve_ipt_matsubara_4th !multi-band !ACTHUNG EXPERIMENTAL!!

  real(8) :: S0,S1,C0,C1,C2,C3

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara normal: 
  ! - half-filling
  !+-------------------------------------------------------------------+
  function solve_ipt_matsubara_1b(fg0_iw) result(sigma_iw)
    complex(8),dimension(:)            :: fg0_iw
    complex(8),dimension(size(fg0_iw)) :: sigma_iw
    real(8),dimension(size(fg0_iw))    :: fg0_tau,sigma_tau
    integer                            :: i,Lf
    real(8)                            :: n
    Lf=size(fg0_iw)
    n = 0.5d0
    S0=Uloc*(n-0.5d0)
    S1=Uloc*Uloc*n*(1d0-n)
    C0=0d0
    C1=1d0
    C2=xmu-S0
    C3=S1+(xmu-S0)*(xmu-S0)
    fg0_tau  = f_fft_gf_iw2tau(fg0_iw,beta,[C0,C1,C2,C3])
    forall(i=1:Lf)sigma_tau(i)=Uloc*Uloc*fg0_tau(i)*fg0_tau(Lf-i+1)*fg0_tau(i)
    sigma_iw = f_fft_sigma_tau2iw(sigma_tau,beta,[S0,S1])
  end function solve_ipt_matsubara_1b



  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara normal: 
  ! - away-from-half-filling
  !+-------------------------------------------------------------------+
  function solve_mpt_matsubara(fg0_iw,n,n0,xmu0) result(sigma_iw)
    complex(8),dimension(:)            :: fg0_iw
    complex(8),dimension(size(fg0_iw)) :: sigma_iw
    real(8),dimension(size(fg0_iw))    :: fg0_tau,sigma_tau
    real(8)                            :: n,n0,xmu0
    real(8)                            :: A,B,A1,A2,B1,B2
    integer                            :: i,Lf
    Lf=size(fg0_iw)
    S0=Uloc*(n-0.5d0)
    S1=Uloc*Uloc*n*(1d0-n)
    C0=0d0
    C1=1d0
    C2=xmu-S0
    C3=S1+(xmu-S0)*(xmu-S0)
    fg0_tau  = f_fft_gf_iw2tau(fg0_iw,beta,[C0,C1,C2,C3])
    forall(i=1:Lf)sigma_tau(i)=Uloc*Uloc*fg0_tau(i)*fg0_tau(Lf-i+1)*fg0_tau(i)
    sigma_iw = f_fft_sigma_tau2iw(sigma_tau,beta,[S0,S1])
    A1= n*(1.d0-n)
    A2= n0*(1.d0-n0)
    A = A1/A2
    B1 = (xmu0-xmu) + Uloc*(1.d0-2.d0*n)
    B2 = n0*(1.d0-n0)*Uloc*Uloc
    B  = B1/B2
    sigma_iw = Uloc*(n-0.5d0) + A*sigma_iw/(1.d0-B*sigma_iw)
  end function solve_mpt_matsubara






  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara attractive: 
  ! - half-filling, REAL order parameter.
  ! - away-from-half-filling, REAL order parameter.
  ! - half-filling, CMPLX order parameter.
  ! - away-from-half-filling, CMPLX order parameter.
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara_r(fg0_iw,delta) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    real(8)                                             :: delta
    complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
    real(8),dimension(:),allocatable                    :: calG11t,calG22t,calFt
    real(8),dimension(:),allocatable                    :: sigmat,selft
    integer                                             :: i,LM
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_ipt_sc_matsubara_r: size(input,1)!= 2"
    allocate(calG11(LM),calG11t(LM))
    allocate(calG22(LM),calG22t(LM))
    allocate(calF(LM),calFt(LM))
    allocate(sigmat(LM),selft(LM))
    !
    !GEt all components of the HFB-corrected Weiss-Fields:
    calG11  =  fg0_iw(1,:)
    calG22  = -conjg(fg0_iw(1,:))
    calF    =  fg0_iw(2,:)
    calG11t = f_fft_gf_iw2tau(calG11,beta,[0d0,1d0,0d0,0d0])
    calG22t = f_fft_gf_iw2tau(calG22,beta,[0d0,1d0,0d0,0d0])
    calFt   = f_fft_gf_iw2tau(calF,beta,[0d0,0d0,0d0,0d0])
    !Get the 2nd-order Sigma:
    forall(i=1:LM)
       sigmat(i)= Uloc*Uloc*(calG11t(i)*calG22t(i) - calFt(i)**2)*calG22t(LM-i+1)
       selft(i) =-Uloc*Uloc*(calFt(i)**2           - calG11t(i)*calG22t(i))*calFt(i)
    end forall
    sigma_iw(1,:) = f_fft_sigma_tau2iw(sigmat,beta,[0d0,0d0])
    sigma_iw(2,:) = f_fft_sigma_tau2iw(selft,beta,[0d0,0d0]) - delta
    !
    open(11,file="Sigma_tau.ipt",access="append")
    open(12,file="Self_tau.ipt",access="append")
    do i=1,LM
       write(11,*)(i-1)*beta/dble(LM-1),sigmat(i)
       write(12,*)(i-1)*beta/dble(LM-1),selft(i)
    enddo
    close(11);close(12)
    !
    deallocate(calG11,calG11t)
    deallocate(calG22,calG22t)
    deallocate(calF,calFt)
    deallocate(sigmat,selft)
  end function solve_ipt_sc_matsubara_r




  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order modified perturbation theory in Matsubara for 
  ! the attractice model away, REAL order parameter case
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_matsubara_r(fg0_iw,n,n0,delta,delta0) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
    real(8),dimension(:),allocatable                    :: calG11t,calG22t,calFt
    real(8),dimension(:),allocatable                    :: sigmat,selft
    integer                                             :: i,LM
    real(8)                                             :: n,n0
    real(8)                                             :: delta,delta0
    real(8)                                             :: A,B
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_mipt_sc_matsubara_r: size(input,1)!= 2"
    allocate(calG11(LM),calG11t(LM))
    allocate(calG22(LM),calG22t(LM))
    allocate(calF(LM),calFt(LM))
    allocate(sigmat(LM),selft(LM))
    !
    !GEt all components of the HFB-corrected Weiss-Fields:
    calG11  =  fg0_iw(1,:)
    calG22  = -conjg(fg0_iw(1,:))
    calF    =  fg0_iw(2,:)
    calG11t = f_fft_gf_iw2tau(calG11,beta,[0d0,1d0,0d0,0d0])
    calG22t = f_fft_gf_iw2tau(calG22,beta,[0d0,1d0,0d0,0d0])
    calFt   = f_fft_gf_iw2tau(calF,beta,[0d0,0d0,0d0,0d0])
    !Get the 2nd-order Sigma:
    forall(i=1:LM)
       sigmat(i)= Uloc*Uloc*(calG11t(i)*calG22t(i) - calFt(i)**2)*calG22t(LM-i+1)
       selft(i)= -Uloc*Uloc*(calFt(i)**2 - calG11t(i)*calG22t(i))*calFt(i)
    end forall
    sigma_iw(1,:) = f_fft_sigma_tau2iw(sigmat,beta,[0d0,0d0])
    sigma_iw(2,:) = f_fft_sigma_tau2iw(selft,beta,[0d0,0d0])
    !
    A=Uloc*Uloc*n*(1.d0-n)-delta**2
    B=Uloc*Uloc*n0*(1.d0-n0)-delta0**2
    sigma_iw(1,:) =-Uloc*(n-0.5d0) + sigma_iw(1,:)*A/B
    sigma_iw(2,:) =-delta       + sigma_iw(2,:)*A/B
    !
    open(11,file="Sigma_tau.ipt",access="append")
    open(12,file="Self_tau.ipt",access="append")
    do i=1,LM
       write(11,*)(i-1)*beta/dble(LM-1),sigmat(i)
       write(12,*)(i-1)*beta/dble(LM-1),selft(i)
    enddo
    close(11);close(12)
    !
    deallocate(calG11,calG11t)
    deallocate(calG22,calG22t)
    deallocate(calF,calFt)
    deallocate(sigmat,selft)
  end function solve_mpt_sc_matsubara_r






  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: Solve 2nd order perturbation theory in Matsubara for 
  ! ! the attractive model at half-filling, REAL order parameter case
  ! !+-------------------------------------------------------------------+
  ! function solve_ipt_sc_matsubara_c(fg0_iw,delta) result(sigma_iw)
  !   complex(8),dimension(:,:)                           :: fg0_iw
  !   complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
  !   complex(8)                                          :: delta
  !   complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
  !   real(8),dimension(:),allocatable                    :: calG11t,calG22t
  !   real(8),dimension(:),allocatable                    :: sigmat
  !   complex(8),dimension(:),allocatable                 :: calFt,selft
  !   integer                                             :: i,LM
  !   LM=size(fg0_iw,2)
  !   if(size(fg0_iw,1)/=2)stop "solve_ipt_sc_matsubara_c: size(input,1)!= 2"
  !   allocate(calG11(LM),calG11t(0:LM))
  !   allocate(calG22(LM),calG22t(0:LM))
  !   allocate(calF(LM),calFt(0:LM))
  !   allocate(sigmat(0:LM),selft(0:LM))
  !   !
  !   !Get the HF-corrected Weiss-Fields:
  !   calG11 =  fg0_iw(1,:)
  !   calG22 = -conjg(fg0_iw(1,:))
  !   calF   =  fg0_iw(2,:)
  !   call fftgf_iw2tau(calG11,calG11t(0:),beta)
  !   call fftgf_iw2tau(calG22,calG22t(0:),beta)
  !   call fftff_iw2tau(calF,calFt(0:),beta)
  !   !Get the 2nd-order Sigma:
  !   forall(i=0:LM)
  !      sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - abs(calFt(i))**2)*calG22t(LM-i)
  !      selft(i) = -U**2*(abs(calFt(i))**2 - calG11t(i)*calG22t(i))*calFt(i) !ACTHUNG HERE: time inversion simmetry
  !   end forall
  !   call fftgf_tau2iw(sigmat(0:),sigma_iw(1,:),beta)
  !   call fftff_tau2iw(selft(0:),sigma_iw(2,:),beta)
  !   !
  !   sigma_iw(2,:)=sigma_iw(2,:) - delta
  !   !
  !   open(11,file="Sigma_tau.ipt",access="append")
  !   open(12,file="Self_tau.ipt",access="append")
  !   do i=0,LM
  !      write(11,*)i*beta/dble(LM),sigmat(i)
  !      write(12,*)i*beta/dble(LM),dreal(selft(i)),dimag(selft(i))
  !   enddo
  !   close(11);close(12)
  !   !
  !   deallocate(calG11,calG11t)
  !   deallocate(calG22,calG22t)
  !   deallocate(calF,calFt)
  !   deallocate(sigmat,selft)
  ! end function solve_ipt_sc_matsubara_c



  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: Solve 2nd order modified perturbation theory in Matsubara for 
  ! ! the attractice model away, REAL order parameter case
  ! !+-------------------------------------------------------------------+
  ! function solve_mpt_sc_matsubara_c(fg0_iw,n,n0,delta,delta0) result(sigma_iw)
  !   complex(8),dimension(:,:)                           :: fg0_iw
  !   complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
  !   complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
  !   real(8),dimension(:),allocatable                    :: calG11t,calG22t
  !   real(8),dimension(:),allocatable                    :: sigmat
  !   complex(8),dimension(:),allocatable                 :: calFt,selft
  !   integer                                             :: i,LM
  !   real(8)                                             :: n,n0
  !   complex(8)                                          :: delta,delta0
  !   real(8)                                             :: A,B
  !   LM=size(fg0_iw,2)
  !   if(size(fg0_iw,1)/=2)stop "solve_ipt_sc_matsubara_c: size(input,1)!= 2"
  !   allocate(calG11(LM),calG11t(0:LM))
  !   allocate(calG22(LM),calG22t(0:LM))
  !   allocate(calF(LM),calFt(0:LM))
  !   allocate(sigmat(0:LM),selft(0:LM))
  !   !
  !   !Get the HF-corrected Weiss-Fields:
  !   calG11 =  fg0_iw(1,:)
  !   calG22 = -conjg(fg0_iw(1,:))
  !   calF   =  fg0_iw(2,:)
  !   call fftgf_iw2tau(calG11,calG11t(0:),beta)
  !   call fftgf_iw2tau(calG22,calG22t(0:),beta)
  !   call fftff_iw2tau(calF,calFt(0:),beta)
  !   !Get the 2nd-order Sigma:
  !   forall(i=0:LM)
  !      sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - abs(calFt(i))**2)*calG22t(LM-i)
  !      selft(i) = -U**2*(abs(calFt(i))**2 - calG11t(i)*calG22t(i))*calFt(i) !ACTHUNG HERE: time inversion simmetry
  !   end forall
  !   call fftgf_tau2iw(sigmat(0:),sigma_iw(1,:),beta)
  !   call fftff_tau2iw(selft(0:),sigma_iw(2,:),beta)
  !   !
  !   !This is not obvious, just a guess now but I need to check!!
  !   A=U**2*n*(1.d0-n)-abs(delta)**2
  !   B=U**2*n0*(1.d0-n0)-abs(delta0)**2
  !   sigma_iw(1,:) =-U*(n-0.5d0) + sigma_iw(1,:)*A/B
  !   sigma_iw(2,:) =-delta       + sigma_iw(2,:)*A/B
  !   !
  !   open(11,file="Sigma_tau.ipt",access="append")
  !   open(12,file="Self_tau.ipt",access="append")
  !   do i=0,LM
  !      write(11,*)i*beta/dble(LM),sigmat(i)
  !      write(12,*)i*beta/dble(LM),dreal(selft(i)),dimag(selft(i))
  !   enddo
  !   close(11);close(12)
  !   !
  !   deallocate(calG11,calG11t)
  !   deallocate(calG22,calG22t)
  !   deallocate(calF,calFt)
  !   deallocate(sigmat,selft)
  ! end function solve_mpt_sc_matsubara_c







  !##################################################################
  !##################################################################
  !##################################################################   
  ! DEBUG MULTI-BAND && 4th ORDER
  !##################################################################
  !##################################################################
  !##################################################################     

  ! function solve_ipt_matsubara_mb(fg0_iw) result(sigma_iw)
  !   complex(8),dimension(:,:)                           :: fg0_iw
  !   complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
  !   real(8),dimension(size(fg0_iw,1),0:size(fg0_iw,2))  :: fg0_tau,sigma_tau
  !   integer                                             :: i,ib,Lf,Nb
  !   Nb=size(fg0_iw,1)
  !   !<DEBUG
  !   if(Nb>2)stop "Nband>2 is not yet implemented... ask developer"
  !   !>DEBUG
  !   Lf=size(fg0_iw,2)
  !   do ib=1,Nb
  !      call fftgf_iw2tau(fg0_iw(ib,:),fg0_tau(ib,0:),beta)
  !   enddo
  !   !Get contribution of the first diagram (1s,1sbar,1sbar)/(2s,2sbar,2sbar)
  !   forall(ib=1:Nb,i=0:Lf)sigma_tau(ib,i)=U*U*fg0_tau(ib,i)*fg0_tau(ib,Lf-i)*fg0_tau(ib,i)
  !   !Get contribution of the second class of 
  !   !diagrams (1s,2s,2s)+(1s,2sbar,2sbar)/(2s,1s,1s)+(2s,1sbar,1sbar)
  !   if(Ust/=0.d0)then
  !      do ib=1,Nb
  !         do i=0,Lf
  !            sigma_tau(ib,i)=sigma_tau(ib,i)+&
  !                 2.d0*Ust*Ust*fg0_tau(ib,i)*fg0_tau(3-ib,Lf-i)*fg0_tau(3-ib,i)
  !         enddo
  !      enddo
  !   endif
  !   do ib=1,Nb
  !      call fftgf_tau2iw(sigma_tau(ib,0:),sigma_iw(ib,:),beta)
  !   enddo
  !   open(100,file="Sigma_tau.ipt")
  !   do i=0,Lf
  !      write(100,"(5F20.12)")i*beta/dble(Lf),(sigma_tau(ib,i),ib=1,Nb)
  !   enddo
  !   close(100)
  ! end function solve_ipt_matsubara_mb


  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: Solve 4th order perturbation theory in Matsubara for 
  ! ! the repulsive model at half-filling
  ! !+-------------------------------------------------------------------+
  ! function solve_ipt_matsubara_4th(fg0_iw) result(sigma_iw)
  !   complex(8),dimension(:)            :: fg0_iw
  !   complex(8),dimension(size(fg0_iw)) :: sigma_iw
  !   real(8),dimension(:),allocatable   :: fg0_tau
  !   real(8),dimension(:),allocatable   :: sigma_tau
  !   real(8),dimension(:),allocatable   :: sigma_tau_a,sigma_tau_b,sigma_tau_c,sigma_tau_d
  !   real(8),dimension(:),allocatable   :: ker,intx,kerx,kery,chi_a,chi_b
  !   integer :: itau,ix,iy
  !   real(8) :: dtau
  !   L=size(fg0_iw)
  !   dtau = beta/dble(L)
  !   !Get G_0(tau)
  !   allocate(fg0_tau(0:L))
  !   allocate(sigma_tau(0:L))
  !   call fftgf_iw2tau(fg0_iw,fg0_tau(0:),beta)
  !   allocate(kerx(0:L),kery(0:L),intx(0:L),ker(0:L))
  !   allocate(sigma_tau_a(0:L),chi_a(0:L))
  !   !get \Sigma^(4a):
  !   do itau=0,L
  !      do ix=0,itau
  !         kerx(ix) =  fg0_tau(itau-ix)**2*fg0_tau(ix)**2
  !      enddo
  !      do ix=itau,L
  !         kerx(ix) = (-fg0_tau(L+itau-ix))**2*fg0_tau(ix)**2
  !      enddo
  !      chi_a(itau) = trapz(dtau,kerx(0:L))
  !   enddo
  !   do itau=0,L
  !      do ix=0,itau
  !         kerx(ix) =  fg0_tau(itau-ix)**2*chi_a(ix)
  !      enddo
  !      do ix=itau,L
  !         kerx(ix) = (-fg0_tau(L+itau-ix))**2*chi_a(ix)
  !      enddo
  !      sigma_tau_a(itau) = U*U*U*U*fg0_tau(itau)*trapz(dtau,kerx(0:L))
  !   enddo
  !   !get \Sigma^(4b):
  !   allocate(sigma_tau_b(0:L),chi_b(0:L))    
  !   do itau=0,L
  !      do ix=0,itau
  !         kerx(ix) =  fg0_tau(itau-ix)**3*fg0_tau(ix)
  !      enddo
  !      do ix=itau,L
  !         kerx(ix) = (-fg0_tau(L+itau-ix))**3*fg0_tau(ix)
  !      enddo
  !      chi_b(itau) = trapz(dtau,kerx(0:L))
  !   enddo
  !   do itau=0,L
  !      do ix=0,itau
  !         kerx(ix) =  fg0_tau(itau-ix)*chi_b(ix)
  !      enddo
  !      do ix=itau,L
  !         kerx(ix) = (-fg0_tau(L+itau-ix))*chi_b(ix)
  !      enddo
  !      sigma_tau_b(itau) = U*U*U*U*fg0_tau(itau)**2*trapz(dtau,kerx(0:L))
  !   enddo
  !   !get \Sigma^(4c):
  !   allocate(sigma_tau_c(0:L))
  !   do itau=0,L
  !      do iy=0,itau
  !         kery(iy) =  fg0_tau(itau-iy)**2*fg0_tau(iy)
  !      enddo
  !      do iy=itau,L
  !         kery(iy) = (-fg0_tau(L+itau-iy))**2*fg0_tau(iy)
  !      enddo
  !      !
  !      do ix=0,itau
  !         kerx(ix) =  fg0_tau(itau-ix)*fg0_tau(ix)**2
  !      enddo
  !      do ix=itau,L
  !         kerx(ix) = (-fg0_tau(L+itau-ix))*fg0_tau(ix)**2
  !      enddo
  !      !
  !      do iy=0,L
  !         do ix=0,iy
  !            ker(ix) =  (-fg0_tau(L+ix-iy))*kerx(ix) !there  was an error itau=>iy
  !         enddo
  !         do ix=iy,L
  !            ker(ix) = fg0_tau(ix-iy)*kerx(ix)
  !         enddo
  !         intx(iy) = trapz(dtau,ker(0:L))
  !      enddo
  !      sigma_tau_c(itau) = U*U*U*U*trapz(dtau,kery(0:L)*intx(0:L))
  !   enddo
  !   !get \Sigma^(4d):
  !   allocate(sigma_tau_d(0:L))
  !   do itau=0,L
  !      do iy=0,itau
  !         kery(iy) =  fg0_tau(itau-iy)*fg0_tau(iy)
  !      enddo
  !      do iy=itau,L
  !         kery(iy) = (-fg0_tau(L+itau-iy))*fg0_tau(iy)
  !      enddo
  !      !
  !      do ix=0,itau
  !         kerx(ix) =  fg0_tau(itau-ix)*fg0_tau(ix)
  !      enddo
  !      do ix=itau,L
  !         kerx(ix) = (-fg0_tau(L+itau-ix))*fg0_tau(ix)
  !      enddo
  !      !
  !      do iy=0,L
  !         do ix=0,iy
  !            ker(ix) =  (-fg0_tau(L+ix-iy))**2*kerx(ix)
  !         enddo
  !         do ix=iy,L
  !            ker(ix) = fg0_tau(ix-iy)**2*kerx(ix)
  !         enddo
  !         intx(iy) = trapz(dtau,ker(0:L))
  !      enddo
  !      sigma_tau_d(itau) = -U*U*U*U*fg0_tau(itau)*trapz(dtau,kery(0:L)*intx(0:L))
  !   enddo
  !   forall(itau=0:L)sigma_tau(itau)=U**2*(fg0_tau(itau))**2*fg0_tau(L-itau)
  !   sigma_tau(0:) = sigma_tau(0:) + 3d0*(sigma_tau_a(0:)+sigma_tau_b(0:)+sigma_tau_c(0:)+sigma_tau_d(0:))
  !   call fftgf_tau2iw(sigma_tau(0:),sigma_iw,beta)
  !   open(100,file="Sigma_tau.ipt")
  !   do itau=0,L
  !      write(100,*)itau*beta/dble(L),sigma_tau(itau)
  !   enddo
  !   close(100)
  ! end function solve_ipt_matsubara_4th






end module IPT_MATSUBARA
