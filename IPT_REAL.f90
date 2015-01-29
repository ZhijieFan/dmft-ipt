!###############################################################
! PURPOSE  : Contains common routines for IPT calculations
! using Second Order PErturbation theory in DMFT approximation
! on the real axis:
! IPT_REAL_nornal: normal case
! IPT_REAL_superc: superconducting case
!###############################################################
module IPT_REAL_NORMAL
  USE IPT_VARS_GLOBAL
  USE SF_SPECIAL, only: fermi
  USE SF_INTEGRATE, only: kronig
  implicit none
  private
  real(8),dimension(:),allocatable    :: A0m,A0p,P1,P2
  integer,allocatable,dimension(:,:)  :: iy_m_ix
  integer,save                        :: loop=1
  complex(8),dimension(:),allocatable :: fg0,sigma
  real(8),dimension(:),allocatable    :: wr
  real(8)                             :: n,n0,xmu0,mesh
  integer                             :: MM

  interface ipt_solve_real
     module procedure solve_ipt_sopt
  end interface ipt_solve_real

  interface mpt_solve_real
     module procedure solve_mpt_sopt
  end interface mpt_solve_real

  public :: ipt_solve_real
  public :: mpt_solve_real

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : interface function for the impurity solver (SOPT)
  ! half-filling case
  !+-------------------------------------------------------------------+
  function solve_ipt_sopt(fg0_,wr_) result(sigma_)
    complex(8),dimension(:)          :: fg0_
    complex(8),dimension(size(fg0_)) :: sigma_
    real(8),dimension(size(fg0_))    :: wr_
    MM=size(fg0_)
    if(loop==1)then
       if(.not.allocated(wr))allocate(wr(MM))
       if(.not.allocated(fg0))allocate(fg0(MM))
       if(.not.allocated(sigma))allocate(sigma(MM))
       call get_frequency_index       
    endif
    fg0=fg0_; wr=wr_ ; mesh=abs(wr(2)-wr(1))
    call getAs
    call getPolarization
    call Sopt
    sigma_=sigma
    loop=loop+1
  end function solve_ipt_sopt




  !+-------------------------------------------------------------------+
  !PURPOSE  : interface function for the impurity solver (SOPT)
  ! away from half-filling case
  !+-------------------------------------------------------------------+
  function solve_mpt_sopt(fg0_,wr_,n_,n0_,xmu0_) result(sigma_)
    complex(8),dimension(:)          :: fg0_
    complex(8),dimension(size(fg0_)) :: sigma_
    real(8),dimension(size(fg0_))    :: wr_
    real(8)                          :: A,B,A1,A2,B1,B2,n_,n0_,xmu0_
    MM=size(fg0_)
    if(loop==1) then
       if(.not.allocated(fg0))allocate(fg0(MM))
       if(.not.allocated(sigma))allocate(sigma(MM))
       if(.not.allocated(wr))allocate(wr(MM))
       call get_frequency_index
    endif
    fg0=fg0_; wr=wr_ ; mesh=abs(wr(2)-wr(1))
    n=n_    ; n0=n0_ ; xmu0=xmu0_
    A=0.d0  ; B=0.d0
    call getAs
    call getPolarization
    call Sopt
    A1 = n*(1.d0-n)
    A2 = n0*(1.d0-n0)
    A  = A1/A2
    B1 = (xmu0-xmu) + Uloc(1)*(1.d0-2.d0*n)
    B2 = n0*(1.d0-n0)*Uloc(1)*Uloc(1)
    B  = B1/B2
    sigma_ = Uloc(1)*(n-0.5d0) + A*sigma/(1.d0-B*sigma)
    loop=loop+1
  end function solve_mpt_sopt




  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve 2^nd order perturbation theory
  !+-------------------------------------------------------------------+
  subroutine Sopt
    integer :: ix,iy,iz
    real(8) :: sum1,sum2
    real(8),dimension(MM) :: reS,imS
    do ix=1,MM
       sum1=zero
       sum2=zero
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             sum1=sum1+A0p(MM-iz+1)*P1(iy)*mesh
             sum2=sum2+A0m(MM-iz+1)*P2(iy)*mesh
          end if
       enddo
       imS(ix)=-Uloc(1)*Uloc(1)*(sum1+sum2)*pi
    enddo
    reS = kronig(imS,wr,size(ImS))
    sigma = reS + xi*imS
  end subroutine Sopt

  !PURPOSE  : Get auxiliary array Aplus, Aminus for faster polarization evaluation
  subroutine getAs
    real(8) :: dos(MM)
    dos(:) =-dimag(fg0(:))/pi
    A0p(:) = dos(:)*fermi(wr(:),beta)
    A0m(:) = dos(:)*(1.d0-fermi(wr(:),beta))
  end subroutine getAs

  !PURPOSE  : Get polarization bubbles
  subroutine getPolarization
    integer :: ix,iy,iz    
    P1=zero
    P2=zero
    do ix=1,MM
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             P1(ix)=P1(ix) + A0p(iy)*A0m(iz)*mesh
             P2(ix)=P2(ix) + A0m(iy)*A0p(iz)*mesh
          endif
       enddo
    enddo
  end subroutine getPolarization

  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    if(.not.allocated(iy_m_ix))allocate(iy_m_ix(MM,MM))
    iy_m_ix=0
    do ix=1,MM
       do iy=1,MM
          iz = iy - ix + MM/2 
          if(iz<1 .OR. iz>MM) iz=-1 !out of range-> if(iz>-L)
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0m))allocate(A0m(MM))
    if(.not.allocated(A0p))allocate(A0p(MM))
    if(.not.allocated(P1)) allocate(P1(MM))
    if(.not.allocated(P2)) allocate(P2(MM))
  end subroutine get_frequency_index

END MODULE IPT_REAL_NORMAL





!-----------------------------------------------------------------------
!SUPERCONDUCTING CASE:
!-----------------------------------------------------------------------
MODULE IPT_REAL_SUPERC
  USE IPT_VARS_GLOBAL
  USE SF_SPECIAL, only: fermi
  USE SF_INTEGRATE, only: kronig
  implicit none
  private

  real(8),dimension(:),allocatable      :: A0p11,A0m11,A0p22,A0m22,B0p,B0m,C0p,C0m
  real(8),dimension(:),allocatable      :: P1,P2,Q1,Q2,R1,R2,T1,T2
  real(8),dimension(:),allocatable      :: wr
  integer,allocatable,dimension(:,:)    :: iy_m_ix
  integer,save                          :: loop=1
  complex(8),allocatable,dimension(:,:) :: fg0,sigma
  real(8)                               :: n,n0,delta,delta0,mesh
  integer                               :: Lm

  interface ipt_solve_real_sc
     module procedure solve_ipt_sc_sopt
  end interface ipt_solve_real_sc

  interface mpt_solve_real_sc
     module procedure solve_mpt_sc_sopt
  end interface mpt_solve_real_sc

  public :: ipt_solve_real_sc
  public :: mpt_solve_real_sc  

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_sopt(fg0_,wr_,delta_,L) result(sigma_)
    integer                   :: L
    complex(8),dimension(2,L) :: fg0_
    complex(8),dimension(2,L) :: sigma_
    real(8),dimension(L)      :: wr_
    real(8)                   :: delta_
    Lm=L
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,Lm))
       if(.not.allocated(sigma))allocate(sigma(2,Lm))
       if(.not.allocated(wr))allocate(wr(Lm))
       call get_frequency_index
    endif
    fg0=fg0_ ; delta=delta_ ; wr=wr_ ; mesh=abs(wr(2)-wr(1))
    call getAs
    call getPolarization
    call Sopt
    sigma_(1,:) =          sigma(1,:)
    sigma_(2,:) = -delta + sigma(2,:)
    loop=loop+1
  end function solve_ipt_sc_sopt

  !PURPOSE  : 
  function solve_mpt_sc_sopt(fg0_,wr_,n_,n0_,delta_,delta0_,L) result(sigma_)
    integer                   :: L
    complex(8),dimension(2,L) :: fg0_
    complex(8),dimension(2,L) :: sigma_
    real(8),dimension(L)      :: wr_
    real(8)                   :: n_,n0_,delta_,delta0_
    real(8)                   :: A,B
    Lm=L
    if(loop==1) then
       if(.not.allocated(fg0))allocate(fg0(2,Lm))
       if(.not.allocated(sigma))allocate(sigma(2,Lm))
       if(.not.allocated(wr))allocate(wr(Lm))
       call get_frequency_index  
    endif
    fg0=fg0_  ; wr=wr_ 
    n=n_ ; n0=n0_ ; delta=delta_ ; delta0=delta0_
    mesh=abs(wr(2)-wr(1))
    call getAs
    call getPolarization
    call Sopt
    A = Uloc(1)*Uloc(1)*n*(1.d0-n)-delta**2
    B = Uloc(1)*Uloc(1)*n0*(1.d0-n0)-delta0**2
    sigma_(1,:) = -Uloc(1)*(n-0.5d0)+ sigma(1,:)*A/B
    sigma_(2,:) = -delta         + sigma(2,:)*A/B
    loop=loop+1
  end function solve_mpt_sc_sopt



  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve 2^nd order perturbation theory
  !+-------------------------------------------------------------------+
  subroutine Sopt
    integer :: ix,iy,iz
    real(8) :: sumP1,sumP2
    real(8) :: sumQ1,sumQ2
    real(8) :: sumR1,sumR2
    real(8) :: sumT1,sumT2
    real(8),dimension(Lm) :: reSig,imSig
    real(8),dimension(Lm) :: reS,imS
    do ix=1,Lm
       sumP1=zero
       sumP2=zero
       sumQ1=zero
       sumQ2=zero
       sumR1=zero
       sumR2=zero
       sumT1=zero
       sumT2=zero
       do iy=1,Lm
          iz= iy_m_ix(iy,ix)
          if((iz>=1).and.(iz<=Lm)) then! in questo modo il contributo fuori da dove le polarizzazioni sono definite e' eliminato 
             sumP1=sumP1 + A0m22(Lm+1-iz)*P1(iy)*mesh ! Adriano mi ha spiegato il senso di com'era prima ma  
             sumP2=sumP2 + A0p22(Lm+1-iz)*P2(iy)*mesh ! questo mi sembra MORALMENTE 
             sumQ1=sumQ1 + B0m(Lm+1-iz)*Q1(iy)*mesh   ! piu' chiaro per i posteri per cui... 
             sumQ2=sumQ2 + B0p(Lm+1-iz)*Q2(iy)*mesh
             sumR1=sumR1 + C0m(Lm+1-iz)*R1(iy)*mesh!C0m(iy)*R1(iz)*mesh   !! DOUBLE CHECK
             sumR2=sumR2 + C0p(Lm+1-iz)*R2(iy)*mesh!C0p(iy)*R2(iz)*mesh
             sumT1=sumT1 + A0m22(Lm+1-iz)*T1(iy)*mesh
             sumT2=sumT2 + A0p22(Lm+1-iz)*T2(iy)*mesh
          endif
       enddo
       imSig(ix)= -Uloc(1)*Uloc(1)*(sumP1 + sumP2 - sumQ1 - sumQ2)*pi
       imS(ix)  = -Uloc(1)*Uloc(1)*(sumR1 + sumR2 - sumT1 - sumT2)*pi
    enddo
    reSig = kronig(imSig,wr,size(imSig))
    reS   = kronig(imS,wr,size(imS))
    sigma(1,:) = reSig + xi*imSig
    sigma(2,:) = reS   + xi*imS
  end subroutine Sopt

  !PURPOSE  : Get auxiliary array Aplus, Aminus for faster polarization evaluation
  subroutine getAs
    integer :: i
    real(8) :: w,dos11(Lm),dos22(Lm),dosF1(Lm),dosF2(Lm)
    do i=1,Lm
       dos11(i) = -dimag(fg0(1,i))/pi
       dos22(i) = -dimag(-conjg(fg0(1,Lm+1-i)))/pi
       dosF1(i) = -dimag(fg0(2,i))/pi
       dosF2(i) = -dimag(fg0(2,i))/pi
    enddo
    A0p11(:) = dos11*fermi(wr,beta)
    A0m11(:) = dos11*(1.d0-fermi(wr,beta))
    A0p22(:) = dos22*fermi(wr,beta)
    A0m22(:) = dos22*(1.d0-fermi(wr,beta))
    B0p(:) = dosF1*fermi(wr,beta)
    B0m(:) = dosF1*(1.d0-fermi(wr,beta))
    C0p(:) = dosF2*fermi(wr,beta)
    C0m(:) = dosF2*(1.d0-fermi(wr,beta))
  end subroutine getAs

  !PURPOSE  : Get polarization bubbles
  subroutine getPolarization
    integer :: ix,iy,iz
    P1=zero
    P2=zero
    Q1=zero
    Q2=zero
    R1=zero
    R2=zero
    T1=zero
    T2=zero
    do ix=1,Lm
       do iy=1,Lm
          iz= iy_m_ix(iy,ix)
          if((iz>=1).and.(iz<=Lm)) then
             P2(ix)=P2(ix) + A0p11(iy)*A0m22(iz)*mesh
             P1(ix)=P1(ix) + A0m11(iy)*A0p22(iz)*mesh
             Q2(ix)=Q2(ix) + C0p(iy)*A0m22(iz)*mesh     
             Q1(ix)=Q1(ix) + C0m(iy)*A0p22(iz)*mesh     
             R2(ix)=R2(ix) + B0p(iy)*B0m(iz)*mesh
             R1(ix)=R1(ix) + B0m(iy)*B0p(iz)*mesh
             T2(ix)=T2(ix) + A0p11(iy)*B0m(iz)*mesh
             T1(ix)=T1(ix) + A0m11(iy)*B0p(iz)*mesh
          endif
       enddo
    enddo
  end subroutine getPolarization

  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    allocate(iy_m_ix(Lm,Lm))
    iy_m_ix=0
    do ix=1,Lm
       do iy=1,Lm
          iz = iy - ix + Lm/2
          ! if(iz<-Lm .OR. iz>Lm) iz=-1 !out of range, the same as the old if(iz>-L)
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0p11))allocate(A0p11(Lm))
    if(.not.allocated(A0m11))allocate(A0m11(Lm))
    if(.not.allocated(A0p22))allocate(A0p22(Lm))
    if(.not.allocated(A0m22))allocate(A0m22(Lm))
    if(.not.allocated(B0p))allocate(B0p(Lm))
    if(.not.allocated(B0m))allocate(B0m(Lm))
    if(.not.allocated(C0p))allocate(C0p(Lm))
    if(.not.allocated(C0m))allocate(C0m(Lm))
    if(.not.allocated(P1))allocate(P1(Lm))
    if(.not.allocated(P2))allocate(P2(Lm))
    if(.not.allocated(Q1))allocate(Q1(Lm))
    if(.not.allocated(Q2))allocate(Q2(Lm))
    if(.not.allocated(R1))allocate(R1(Lm))
    if(.not.allocated(R2))allocate(R2(Lm))
    if(.not.allocated(T1))allocate(T1(Lm))
    if(.not.allocated(T2))allocate(T2(Lm))
  end subroutine get_frequency_index

end module IPT_REAL_SUPERC



!-----------------------------------------------------------------------
!FINAL MODULE
!-----------------------------------------------------------------------
MODULE IPT_REAL
  USE IPT_REAL_NORMAL
  USE IPT_REAL_SUPERC
  implicit none
  private
  public :: ipt_solve_real
  public :: ipt_solve_real_sc
  public :: mpt_solve_real
  public :: mpt_solve_real_sc  
END MODULE IPT_REAL
