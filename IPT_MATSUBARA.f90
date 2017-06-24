!###############################################################
!     PURPOSE  : SOLVE DMFT-IPT REPULSIVE IN MATSUBARA FREQUENCY
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_MATSUBARA
  USE IPT_VARS_GLOBAL
  USE IPT_GF
  USE SF_LINALG, only: inv
  USE SF_CONSTANTS, only: xi,pi,one
  USE SF_ARRAYS, only:arange
  USE SF_SPECIAL, only: bethe_lattice
  USE SF_IOTOOLS, only: free_unit
  implicit none
  private

  interface ipt_solve_matsubara
     module procedure solve_ipt_matsubara
     module procedure solve_ipt_sc_matsubara
  end interface ipt_solve_matsubara

  interface mpt_solve_matsubara
     module procedure solve_mpt_matsubara
     module procedure solve_mpt_sc_matsubara
  end interface mpt_solve_matsubara

  interface ipt_measure_dens_matsubara
     module procedure ipt_measure_dens_matsubara_SW,ipt_measure_dens_matsubara_G
  end interface ipt_measure_dens_matsubara

  interface ipt_measure_energy_matsubara
     module procedure ipt_measure_energy_matsubara_bethe,ipt_measure_energy_matsubara_hk
  end interface ipt_measure_energy_matsubara


  !half-filling solver:
  public :: ipt_solve_matsubara
  public :: mpt_solve_matsubara

  public :: ipt_measure_potential_energy_matsubara
  public :: ipt_measure_kinetic_energy_matsubara
  public :: ipt_measure_hartree_energy_matsubara
  public :: ipt_measure_dens_matsubara
  public :: ipt_measure_phi_matsubara
  public :: ipt_measure_docc_matsubara
  public :: ipt_measure_zeta_matsubara
  !
  public :: ipt_measure_observables_matsubara
  public :: ipt_measure_energy_matsubara


  ! real(8) :: S0,S1,C0,C1,C2,C3


contains



  !##################################################################
  !##################################################################
  !##################################################################   
  ! SOLVER ROUTINES: NORMAL AND SUPERCONDUCTIVE CASES
  !##################################################################
  !##################################################################
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara normal: 
  ! - half-filling
  ! - away-from-half-filling
  !+-------------------------------------------------------------------+
  function solve_ipt_matsubara(fg0_iw) result(sigma_iw)
    complex(8),dimension(:)            :: fg0_iw
    complex(8),dimension(size(fg0_iw)) :: sigma_iw
    real(8),dimension(0:size(fg0_iw))  :: fg0_tau,sigma_tau
    integer                            :: i,Lf,unit
    real(8)                            :: n
    Lf=size(fg0_iw)
    n = 0.5d0
    call fft_iw2tau(fg0_iw,fg0_tau(0:),beta)
    forall(i=0:Lf)sigma_tau(i)=Uloc(1)*Uloc(1)*fg0_tau(i)*fg0_tau(Lf-i)*fg0_tau(i)
    call fft_tau2iw(sigma_tau(0:),sigma_iw,beta)
    unit=free_unit()
    open(unit,file="Sigma_tau.ipt")
    do i=0,Lf
       write(unit,*)i*beta/Lf,sigma_tau(i)
    enddo
    close(unit)
  end function solve_ipt_matsubara

  function solve_mpt_matsubara(fg0_iw,n,n0,xmu0) result(sigma_iw)
    complex(8),dimension(:)            :: fg0_iw
    complex(8),dimension(size(fg0_iw)) :: sigma_iw
    real(8),dimension(0:size(fg0_iw))  :: fg0_tau,sigma_tau
    real(8)                            :: n,n0,xmu0
    real(8)                            :: A,B,A1,A2,B1,B2
    integer                            :: i,Lf
    Lf=size(fg0_iw)
    call fft_iw2tau(fg0_iw,fg0_tau(0:),beta)
    forall(i=0:Lf)sigma_tau(i)=Uloc(1)*Uloc(1)*fg0_tau(i)*fg0_tau(Lf-i)*fg0_tau(i)
    call fft_tau2iw(sigma_tau(0:),sigma_iw,beta)
    A1= n*(1.d0-n)
    A2= n0*(1.d0-n0)
    A = A1/A2
    B1 = (xmu0-xmu) + Uloc(1)*(1.d0-2.d0*n)
    B2 = n0*(1.d0-n0)*Uloc(1)*Uloc(1)
    B  = B1/B2
    sigma_iw = Uloc(1)*(n-0.5d0) + A*sigma_iw/(1.d0-B*sigma_iw)
    open(11,file="Sigma_tau.ipt")
    do i=0,Lf
       write(11,*)i*beta/Lf,sigma_tau(i)
    enddo
    close(11)
  end function solve_mpt_matsubara







  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara attractive: 
  !   REAL order parameter case
  ! - half-filling, REAL order parameter.
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara(fg0_iw,delta) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    real(8)                                             :: delta
    complex(8),dimension(size(fg0_iw,2))                :: calG11,calG22,calF
    real(8),dimension(0:size(fg0_iw,2))                 :: calG11t,calG22t,calFt
    real(8),dimension(0:size(fg0_iw,2))                 :: sigmat,selft
    integer                                             :: i,LM
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_ipt_sc_matsubara: size(input,1)!= 2"
    calG11  =  fg0_iw(1,:)
    calG22  = -conjg(fg0_iw(1,:))
    calF    =  fg0_iw(2,:)
    call fft_iw2tau(calG11,calG11t(0:),beta)
    call fft_iw2tau(calG22,calG22t(0:),beta)
    call fft_iw2tau(calF,calFt(0:),beta,notail=.true.)
    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)= Uloc(1)*Uloc(1)*(calG11t(i)*calG22t(i) - calFt(i)*calFt(i))*calG22t(LM-i)
       selft(i) =-Uloc(1)*Uloc(1)*(calFt(i)*calFt(i)     - calG11t(i)*calG22t(i))*calFt(i)
    end forall
    call fft_tau2iw(Sigmat(0:),sigma_iw(1,:),beta)
    call fft_tau2iw(Selft(0:),sigma_iw(2,:),beta)
    sigma_iw(2,:)=sigma_iw(2,:)-delta
    !
    open(11,file="Sigma_tau.ipt")
    open(12,file="Self_tau.ipt")
    do i=0,LM
       write(11,*)i*beta/LM,sigmat(i)
       write(12,*)i*beta/LM,selft(i)
    enddo
    close(11);close(12)
  end function solve_ipt_sc_matsubara

  !PURPOSE: Solve 2nd order perturbation theory in Matsubara attractive: 
  !   REAL order parameter case
  ! - away-from-half-filling, REAL order parameter.
  function solve_mpt_sc_matsubara(fg0_iw,n,n0,delta,delta0) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    complex(8),dimension(size(fg0_iw,2))                :: calG11,calG22,calF
    real(8),dimension(0:size(fg0_iw,2))                 :: calG11t,calG22t,calFt
    real(8),dimension(0:size(fg0_iw,2))                 :: sigmat,selft
    integer                                             :: i,LM
    real(8)                                             :: n,n0
    real(8)                                             :: delta,delta0
    real(8)                                             :: A,B
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_mipt_sc_matsubara_r: size(input,1)!= 2"
    !GEt all components of the HFB-corrected Weiss-Fields:
    calG11  =  fg0_iw(1,:)
    calG22  = -conjg(fg0_iw(1,:))
    calF    =  fg0_iw(2,:)
    call fft_iw2tau(calG11,calG11t(0:),beta)
    call fft_iw2tau(calG22,calG22t(0:),beta)
    call fft_iw2tau(calF,calFt(0:),beta,notail=.true.)
    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)= Uloc(1)*Uloc(1)*(calG11t(i)*calG22t(i) - calFt(i)*calFt(i))*calG22t(LM-i)
       selft(i) =-Uloc(1)*Uloc(1)*(calFt(i)*calFt(i)     - calG11t(i)*calG22t(i))*calFt(i)
    end forall
    call fft_tau2iw(Sigmat(0:),sigma_iw(1,:),beta)
    call fft_tau2iw(Selft(0:),sigma_iw(2,:),beta)
    !
    A=Uloc(1)*Uloc(1)*n*(1.d0-n)-delta**2
    B=Uloc(1)*Uloc(1)*n0*(1.d0-n0)-delta0**2
    sigma_iw(1,:) =-Uloc(1)*(n-0.5d0) + sigma_iw(1,:)*A/B
    sigma_iw(2,:) =-delta       + sigma_iw(2,:)*A/B
    !
    open(11,file="Sigma_tau.ipt")
    open(12,file="Self_tau.ipt")
    do i=0,LM
       write(11,*)i*beta/LM,sigmat(i)
       write(12,*)i*beta/LM,selft(i)
    enddo
    close(11);close(12)
    !
  end function solve_mpt_sc_matsubara






  !##################################################################
  !##################################################################
  !##################################################################   
  ! MEASURE OBSERVABLES
  !##################################################################
  !##################################################################
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE: measure observables in Matsubara formalism
  !+-------------------------------------------------------------------+
  function ipt_measure_observables_matsubara(Sigma,Weiss) result(obs)
    complex(8),dimension(:)          :: Sigma
    complex(8),dimension(size(Sigma)):: Weiss
    real(8),dimension(2)             :: obs
    obs(1) = ipt_measure_docc_matsubara(Sigma,Weiss)
    obs(2) = ipt_measure_zeta_matsubara(Sigma,Weiss)
  end function ipt_measure_observables_matsubara

  !PURPOSE: measure density
  function ipt_measure_dens_matsubara_SW(Sigma,Weiss) result(ipt_dens)
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss,Green
    real(8),dimension(0:size(Sigma)) :: Gtau
    real(8)                           :: ipt_dens
    integer :: i
    Green = one/Weiss - Sigma
    Green = one/Green
    call fft_iw2tau(Green,Gtau(0:),beta)
    do i=0,size(Sigma)
       write(450,*)i*beta/size(Sigma),Gtau(i)
    enddo
    ipt_dens = fft_gbeta_minus(Green,beta)
  end function ipt_measure_dens_matsubara_SW
  function ipt_measure_dens_matsubara_G(Green) result(ipt_dens)
    complex(8),dimension(:)           :: Green
    real(8)                           :: ipt_dens
    ipt_dens = fft_gbeta_minus(Green,beta)
  end function ipt_measure_dens_matsubara_G


  function ipt_measure_phi_matsubara(Green) result(ipt_phi)
    complex(8),dimension(:)           :: Green
    real(8)                           :: ipt_phi
    ipt_phi = fft_gbeta_minus(Green,beta,notail=.true.)
  end function ipt_measure_phi_matsubara

  !PURPOSE: measure renormalization constant zeta
  function ipt_measure_zeta_matsubara(Sigma,Weiss) result(ipt_zeta)
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss
    real(8)                           :: ipt_zeta,wm1
    wm1=pi/beta
    ipt_zeta = 1d0 - dimag(Sigma(1))/wm1
    ipt_zeta = 1d0/ipt_zeta
  end function ipt_measure_zeta_matsubara

  !PURPOSE: measure double occupancy
  function ipt_measure_docc_matsubara(Sigma,Weiss) result(ipt_docc)
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss
    real(8)                           :: ipt_docc,epot,ehar
    epot = ipt_measure_potential_energy_matsubara(Sigma,Weiss)
    ehar = ipt_measure_hartree_energy_matsubara(Sigma,Weiss)
    ipt_docc = 0.25d0
    if(uloc(1) > 0d0)ipt_docc = epot/uloc(1) - ehar/uloc(1)
  end function ipt_measure_docc_matsubara

  !PURPOSE: measure all energies for the Bethe lattice
  function ipt_measure_energy_matsubara_bethe(Sigma,Weiss,Lk,D) result(obs)
    real(8),optional                 :: D
    real(8)                          :: D_
    complex(8),dimension(:)          :: Sigma
    complex(8),dimension(size(Sigma)):: Weiss
    real(8),dimension(3)             :: obs
    real(8),dimension(:),allocatable :: Hk,wtk
    integer                          :: Lk,i
    D_=1d0;if(present(D))D_=D
    allocate(Hk(Lk),Wtk(Lk))
    call bethe_lattice(Wtk,Hk,Lk,D_)
    obs(1) = ipt_measure_kinetic_energy_matsubara(Hk,Wtk,Sigma)
    obs(2) = ipt_measure_potential_energy_matsubara(Sigma,Weiss)
    obs(3) = ipt_measure_hartree_energy_matsubara(Sigma,Weiss)
    deallocate(Hk,Wtk)
  end function ipt_measure_energy_matsubara_bethe
  !PURPOSE: measure all energies for a given Hamiltonian Hk
  function ipt_measure_energy_matsubara_hk(Sigma,Weiss,Hk,Wtk) result(obs)
    real(8),dimension(:)              :: Hk
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss
    real(8),dimension(size(Hk))       :: Wtk
    real(8),dimension(3)              :: obs
    obs(1) = ipt_measure_kinetic_energy_matsubara(Hk,Wtk,Sigma)
    obs(2) = ipt_measure_potential_energy_matsubara(Sigma,Weiss)
    obs(3) = ipt_measure_hartree_energy_matsubara(Sigma,Weiss)
  end function ipt_measure_energy_matsubara_hk

  !PURPOSE: measure potential energy
  function ipt_measure_potential_energy_matsubara(Sigma,Weiss) result(ipt_Epot)
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss,Green
    real(8)                           :: ipt_Epot
    Green = one/Weiss - Sigma
    Green = one/Green
    ipt_Epot=sum(dreal(Sigma)*dreal(Green))
    ipt_Epot=ipt_Epot-sum(dimag(Sigma)*dimag(Green))
    ipt_Epot=ipt_Epot/beta*2d0
  end function ipt_measure_potential_energy_matsubara

  !PURPOSE: measure hartree energy term
  function ipt_measure_hartree_energy_matsubara(Sigma,Weiss) result(ipt_Ehartree)
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss,Green
    real(8)                           :: ipt_Ehartree,n
    Green = one/Weiss - Sigma
    Green = one/Green
    n = fft_gbeta_minus(Green,beta)
    ipt_Ehartree = -Uloc(1)*n + Uloc(1)*0.25d0 
  end function ipt_measure_hartree_energy_matsubara

  !PURPOSE: measure kinetic energy
  function ipt_measure_kinetic_energy_matsubara(Hk,Wtk,Sigma) result(ipt_Ekin)
    real(8),dimension(:)                  :: Hk
    complex(8),dimension(:)               :: Sigma
    real(8),dimension(size(Hk))           :: Wtk
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    real(8)                               :: ipt_Ekin
    Sigma_(1,1,:)=Sigma
    Hk_(1,1,:) = Hk
    ipt_Ekin = f_ipt_kinetic_normal(Hk_,Wtk,Sigma_)
  end function ipt_measure_kinetic_energy_matsubara
  function f_ipt_kinetic_normal(Hk,Wtk,Sigma) result(ipt_Ekin)
    integer                                  :: Lk,No,Liw
    integer                                  :: i,ik,iorb
    complex(8),dimension(:,:,:)              :: Hk
    complex(8),dimension(:,:,:)              :: Sigma
    real(8),dimension(:)                     :: Wtk
    !
    real(8),dimension(:,:),allocatable       :: Sigma_HF
    real(8),dimension(:),allocatable         :: wm
    complex(8),dimension(:,:),allocatable    :: Ak,Bk
    complex(8),dimension(:,:),allocatable    :: Ck,Zk
    complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
    real(8)                                  :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                                  :: H0,ipt_Ekin
    !
    No = size(Hk,1)
    Lk = size(Hk,3)
    Liw= size(Sigma,3)
    if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
    if(No/=size(Sigma,1).OR.No/=size(Sigma,2))stop "get_kinetic_energy: size(Sigma,1/2)!=size(Hk,1) [Norb_total]"
    if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
    !
    allocate(wm(Liw))
    allocate(Sigma_HF(No,No))
    allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No))
    !
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    H0=0d0
    Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik) - Sigma(:,:,i)
          select case(No)
          case default
             call inv(Gk)
          case(1)
             Gk = 1d0/Gk
          end select
          Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Ak,Gk - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
       enddo
    enddo
    spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak= Hk(:,:,ik)
       Bk=-Hk(:,:,ik)-Sigma_HF(:,:)
       Ck= matmul(Ak,Bk)
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ipt_Ekin=H0+Tail0+Tail1
    deallocate(wm,Sigma_HF,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
  contains
    function trace_matrix(M,dim) result(tr)
      integer                       :: dim
      complex(8),dimension(dim,dim) :: M
      complex(8) :: tr
      integer                       :: i
      tr=dcmplx(0d0,0d0)
      do i=1,dim
         tr=tr+M(i,i)
      enddo
    end function trace_matrix
  end function f_ipt_kinetic_normal









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







  !##################################################################
  !##################################################################
  !##################################################################   
  ! SUPERCONDUCTING CASE WITH COMPLEX ORDER PARAMETER.
  !##################################################################
  !##################################################################
  !##################################################################     


  ! !+-------------------------------------------------------------------+
  ! !PURPOSE: Solve 2nd order perturbation theory in Matsubara attractive: 
  ! ! - half-filling, CMPLX order parameter.
  ! ! - away-from-half-filling, CMPLX order parameter.
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






end module IPT_MATSUBARA
