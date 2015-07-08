module IPT_KELDYSH
  use IPT_GF
  use IPT_VARS_GLOBAL
  use SF_SPECIAL
  use SF_ARRAYS
  USE SF_IOTOOLS

  implicit none
  private

  integer             :: M
  integer,save        :: loop=1 
  type(keldysh_gf)    :: fg0,sigma
  type(keldysh_gf)    :: fg0k(2),sk(2),calG11,calG22,calF12,calF21
  real(8)             :: dt_,dw_
  real(8),allocatable :: t_(:),wr_(:)

  public :: ipt_solve_keldysh
  public :: ipt_solve_keldysh_sc

contains





  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function ipt_solve_keldysh(fg0_) result(sigma_)
    complex(8),dimension(:)          :: fg0_
    complex(8),dimension(size(fg0_)) :: sigma_
    real(8)                          :: tmax_
    integer                          :: i
    real(8)                          :: A
    M=size(fg0_)
    if(loop==1)then
       if(.not.fg0%status)  call allocate_gf(fg0,M)
       if(.not.sigma%status)call allocate_gf(sigma,M)
    endif
    fg0%ret%w = fg0_
    dt_       = pi/wmax
    tmax_     = dt_*M/2
    allocate(wr_(M),t_(0:M))
    wr_  = linspace(-wmax,wmax,M,mesh=dw_)
    t_   = linspace(-tmax_,tmax_,M+1,mesh=dt_)
    !
    !FFT to real time:
    !G0^{<,>}(t) = FT[G0^{<,>}(w)]
    do i=1,M
       A = -dimag(fg0%ret%w(i))/pi
       fg0%less%w(i) = pi2*xi*fermi(wr_(i),beta)*A
       fg0%gtr%w(i)  = pi2*xi*(fermi(wr_(i),beta)-1.d0)*A
    enddo
    ! call splot("G0ret_w.1",wr_,fg0%ret%w)
    ! call splot("G0less_w.1",wr_,fg0%less%w)
    ! call splot("G0gtr_w.1",wr_,fg0%gtr%w)
    ! !
    ! fg0%ret%t  =  f_fft_rw2rt(fg0%ret%w)*dw_/pi2
    fg0%less%t =  f_fft_rw2rt(fg0%less%w)*dw_/pi2
    fg0%gtr%t  =  f_fft_rw2rt(fg0%gtr%w)*dw_/pi2
    ! call splot("G0ret_t.1",t_,fg0%ret%t)
    ! call splot("G0less_t.1",t_,fg0%less%t)
    ! call splot("G0gtr_t.1",t_,fg0%gtr%t)
    ! !
    ! fg0%ret%t(:) =heaviside(t_(:))*(fg0%gtr%t(:)-fg0%less%t(:))
    ! call splot("G0ret_t.2",t_,fg0%ret%t)
    ! !
    do i=0,M
       sigma%less%t(i)=Uloc(1)*Uloc(1)*fg0%less%t(i)*fg0%gtr%t(M-i)*fg0%less%t(i)
       sigma%gtr%t(i) =Uloc(1)*Uloc(1)*fg0%gtr%t(i)*fg0%less%t(M-i)*fg0%gtr%t(i)
       sigma%ret%t(i) =heaviside(t_(i))*(sigma%gtr%t(i)-sigma%less%t(i))
    enddo
    if(heaviside(0.d0)==1.d0)sigma%ret%t(M/2)=sigma%ret%t(M/2)/2.d0
    ! call splot("Sret_t.1",t_,Sigma%ret%t)
    ! call splot("Sless_t.1",t_,Sigma%less%t)
    ! call splot("Sgtr_t.1",t_,Sigma%gtr%t)
    sigma%ret%w = f_fft_rt2rw(sigma%ret%t)*dt_
    ! sigma%less%w = f_fft_rt2rw(sigma%less%t)*dt_
    ! sigma%gtr%w = f_fft_rt2rw(sigma%gtr%t)*dt_
    ! call splot("Sret_w.2",wr_,Sigma%ret%w)
    ! call splot("Sless_w.2",wr_,Sigma%less%w)
    ! call splot("Sgtr_w.2",wr_,Sigma%gtr%w)
    sigma_=sigma%ret%w
    loop=loop+1
    deallocate(wr_,t_)
  end function ipt_solve_keldysh




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function ipt_solve_keldysh_sc(fg0_,delta,vbias) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    real(8),optional                                :: vbias
    real(8)                                         :: delta,tmax_,nf,nfL,nfR
    integer                                         :: i
    real(8)                                         :: A,B
    M=size(fg0_,2)
    if(loop==1)then
       if(.not.fg0k(1)%status)  call allocate_gf(fg0k(1),M)
       if(.not.fg0k(2)%status)  call allocate_gf(fg0k(2),M)
       if(.not.calG11%status)  call allocate_gf(calG11,M)
       if(.not.calG22%status)  call allocate_gf(calG22,M)
       if(.not.calF12%status)  call allocate_gf(calF12,M)
       if(.not.calF21%status)  call allocate_gf(calF21,M)
       if(.not.sk(1)%status)call allocate_gf(sk(1),M)
       if(.not.sk(2)%status)call allocate_gf(sk(2),M)
    endif
    dt_       = pi/wmax
    tmax_     = dt_*M/2
    allocate(wr_(M),t_(0:M))
    wr_  = linspace(-wmax,wmax,M,mesh=dw_)
    t_   = linspace(-tmax_,tmax_,M+1,mesh=dt_)
    !
    fg0k(1)=zero ; fg0k(2)=zero
    do i=1,M
       A = -dimag(fg0_(1,i))/pi
       B = -dimag(fg0_(2,i))/pi
       nf= fermi(wr_(i),beta)
       if(present(vbias))then
          nfL= fermi(wr_(i)-vbias/2,beta)!istep(wr(i)-vbias/2.d0)
          nfR= fermi(wr_(i)+vbias/2,beta)!istep(wr(i)+vbias/2.d0)
          nf= (nfL+nfR)/2
       endif
       fg0k(1)%less%w(i) = pi2*xi*nf*A
       fg0k(1)%gtr%w(i)  = pi2*xi*(nf-1.d0)*A
       fg0k(2)%less%w(i) = pi2*xi*nf*B       
       fg0k(2)%gtr%w(i)  = pi2*xi*(nf-1.d0)*B
    enddo
    !
    calG11%less%w = fg0k(1)%less%w    
    calG11%gtr%w  = fg0k(1)%gtr%w    
    calG22%less%w = -(fg0k(1)%gtr%w(M:1:-1))
    calG22%gtr%w  = -(fg0k(1)%less%w(M:1:-1))
    calF12%less%w = -conjg(fg0k(2)%less%w)
    calF12%gtr%w  =  fg0k(2)%gtr%w
    calF21%less%w = fg0k(2)%less%w
    calF21%gtr%w  = -conjg(fg0k(2)%gtr%w)
    !
    calG11%less%t = f_fft_rw2rt(calG11%less%w)*dw_/pi2
    calG11%gtr%t  = f_fft_rw2rt(calG11%gtr%w)*dw_/pi2
    calG22%less%t = f_fft_rw2rt(calG22%less%w)*dw_/pi2
    calG22%gtr%t  = f_fft_rw2rt(calG22%gtr%w)*dw_/pi2
    calF12%less%t = f_fft_rw2rt(calF12%less%w)*dw_/pi2
    calF12%gtr%t  = f_fft_rw2rt(calF12%gtr%w)*dw_/pi2
    calF21%less%t = f_fft_rw2rt(calF21%less%w)*dw_/pi2
    calF21%gtr%t  = f_fft_rw2rt(calF21%gtr%w)*dw_/pi2

    do i=0,M 
       sk(1)%less%t(i) = Uloc(1)*Uloc(1)*(calG11%less%t(i)*calG22%less%t(i) - calF12%less%t(i)*calF21%less%t(i))*calG22%gtr%t(M-i)
       sk(1)%gtr%t(i)  = Uloc(1)*Uloc(1)*(calG11%gtr%t(i)*calG22%gtr%t(i) - calF12%gtr%t(i)*calF21%gtr%t(i))*calG22%less%t(M-i)
       sk(2)%less%t(i) = Uloc(1)*Uloc(1)*(calF12%less%t(i)*calF21%less%t(i) - calG11%less%t(i)*calG22%less%t(i))*calF12%gtr%t(M-i)
       sk(2)%gtr%t(i)  = Uloc(1)*Uloc(1)*(calF12%gtr%t(i)*calF21%gtr%t(i)  - calG11%gtr%t(i)*calG22%gtr%t(i))*calF12%less%t(M-i)
       sk(1)%ret%t(i)  = heaviside(t_(i))*(sk(1)%gtr%t(i)-sk(1)%less%t(i))
       sk(2)%ret%t(i)  = heaviside(t_(i))*(sk(2)%gtr%t(i)-sk(2)%less%t(i))
    enddo
    if(heaviside(0.d0)==1.d0)sk(1)%ret%t(0)=sk(1)%ret%t(0)/2.d0 
    if(heaviside(0.d0)==1.d0)sk(2)%ret%t(0)=sk(2)%ret%t(0)/2.d0
    !
    sk(1)%ret%w = f_fft_rt2rw(sk(1)%ret%t)*dt_
    sk(2)%ret%w = f_fft_rt2rw(sk(2)%ret%t)*dt_
    !
    sigma_(1,:) = sk(1)%ret%w
    sigma_(2,:) = -delta + sk(2)%ret%w
    !
    deallocate(wr_,t_)
  end function ipt_solve_keldysh_sc



end module IPT_KELDYSH

