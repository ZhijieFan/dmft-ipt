
!########################################################
!     Program  : HMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt_matsuara_2dsquare
  USE DMFT_IPT
  USE SQUARE_LATTICE
  USE IOTOOLS
  USE SPLINE
  implicit none

  logical                :: converged
  real(8)                :: n,z
  integer                :: i,Lk,iloop
  integer :: p,q
  logical :: upmflag
  complex(8)             :: zeta
  type(matsubara_gf)     :: fg,sigma,fg0
  !complex(8),allocatable :: fg0(:)
  real(8),allocatable    :: wm(:),tau(:),wt(:),epsik(:),nk(:),dtau(:)

  call read_input("inputIPT.in")
  call parse_cmd_variable(p,"P",default=10)
  call parse_cmd_variable(q,"Q",default=100)
  call parse_cmd_variable(upmflag,"UPMFLAG",default=.true.)
  L=2*P*Q
  !allocate functions:

  call allocate_gf(fg,L)
  call allocate_gf(sigma,L)
  call allocate_gf(fg0,L)

  !build freq. array
  allocate(wm(L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)


  allocate(tau(0:L),dtau(0:L))
  if(upmflag)then
     tau(0:) = upminterval(0.d0,beta,beta/2.d0,P,Q,type=0)
     call fftgf_iw2tau_upm(wm,(/(one,i=1,L)/),tau(0:),dtau(0:),beta)
  else
     tau(0:) = linspace(0.d0,beta,L+1)
     call fftgf_iw2tau((/(one,i=1,L)/),dtau(0:),beta)
  endif
  call splot("Delta_tau.ipt",tau(0:),dtau(0:))

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  !dmft loop:
  n=0.5d0
  sigma%iw=u*n
  !sigma%iw=u*(n-0.5d0)
  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta = xi*wm(i) + xmu - sigma%iw(i)
        fg%iw(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo

     if(upmflag)then
        call fftgf_iw2tau_upm(wm,fg%iw,tau(0:),fg%tau(0:),beta)
     else
        call fftgf_iw2tau(fg%iw,fg%tau(0:),beta)
     endif
     n   = -fg%tau(L)
     fg0%iw = one/(one/fg%iw + sigma%iw - u*n)
     !fg0%iw = one/(one/fg%iw + sigma%iw - u*(n-0.5d0))
     if(upmflag)then
        call fftgf_iw2tau_upm(wm,fg0%iw,tau(0:),fg0%tau(0:),beta)
     else
        call fftgf_iw2tau(fg0%iw,fg0%tau,beta)
     endif
     forall(i=0:L)sigma%tau(i)=U**2*(fg0%tau(i))**2*fg0%tau(L-i)+u*n*dtau(i)

     if(upmflag)then
        call fftgf_tau2iw_upm(wm,sigma%iw,tau(0:),sigma%tau(0:),beta)
     else
        call fftgf_tau2iw(sigma%tau,sigma%iw,beta)
     endif
     !sigma%iw = sigma%iw  + u*n
     !sigma%iw = sigma%iw + u*(n-0.5d0)

     !! sigma%iw= solve_ipt_matsubara(fg0) + u*n
     !! sigma%iw= solve_ipt_matsubara(fg0) + u*(n-0.5d0)

     converged=check_convergence(sigma%iw,eps_error,nsuccess,nloop)
     z=1.d0 - dimag(sigma%iw(1))/wm(1);z=1.d0/z
     call splot("nVSiloop.ipt",iloop,n,append=TT)
     call splot("zetaVSiloop.ipt",iloop,z,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  call close_file("zetaVSiloop.ipt")
  call splot("G_iw.ipt",wm,fg%iw,append=printf)
  call splot("G0_iw.ipt",wm,fg0%iw,append=printf)
  call splot("Sigma_iw.ipt",wm,sigma%iw,append=printf)
  call splot("G_tau.ipt",tau(0:),fg%tau(0:),append=printf)
  call splot("G0_tau.ipt",tau(0:),fg0%tau(0:),append=printf)
  call splot("Sigma_tau.ipt",tau(0:),sigma%tau(0:),append=printf)
  nk = square_lattice_momentum_distribution(Lk)
  call splot("nkVSepsk.ipt",epsik,nk,append=printf)

contains

  function square_lattice_momentum_distribution(Lk) result(nk)
    integer            :: Lk
    integer            :: ik,i
    real(8)            :: nk(Lk)
    do ik=1,Lk
       fg%iw=one/(xi*wm - epsik(ik) - sigma%iw)
       call fftgf_iw2tau(fg%iw,fg%tau,beta)
       nk(ik)=-fg%tau(L)
    enddo
  end function square_lattice_momentum_distribution


  subroutine fftgf_iw2tau_upm(wm,gw,tm,gt,beta)
    integer                             :: i,j,N,L
    real(8),dimension(:)                :: wm
    complex(8),dimension(size(wm))      :: gw,tmpGw
    real(8),dimension(:)                :: tm
    real(8),dimension(size(tm))         :: gt
    complex(8)                          :: tail,fg
    real(8)                             :: tau,beta,mues,At,foo
    L=size(wm) ; N=size(tm)
    mues =-real(gw(L),8)*wm(L)**2
    tmpGw=(0.d0,0.d0)
    do i=1,L
       tail=-cmplx(mues,wm(i),8)/(mues**2+wm(i)**2)
       fg  = (0.d0,0.d0)
       if(i<=N)fg  = gw(i)-tail
       tmpGw(i)=fg
    enddo
    do i=1,N-1
       tau=tm(i)
       if(mues >= 0.d0)then
          if((mues*beta) > 30.d0)then
             At = -exp(-mues*tau)
          else
             At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       else
          if((mues*beta) < -30.d0)then
             At = -exp(mues*(beta-tau))
          else
             At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       endif
       foo   = sum(cos(wm(:)*tau)*dreal(tmpGw(:))) + sum(sin(wm(:)*tau)*dimag(tmpGw(:)))
       gt(i) = foo*2.d0/beta + At
    enddo
    gt(N)=-(gt(1)+1.d0)
  end subroutine fftgf_iw2tau_upm


  subroutine fftgf_tau2iw_upm(wm,gw,tm,gt,beta)
    integer                             :: i,j,N,L,NN
    real(8),dimension(:)                :: wm
    complex(8),dimension(size(wm))      :: gw
    real(8),dimension(:)                :: tm
    real(8),dimension(size(tm))         :: gt
    real(8),allocatable                 :: tmpGt(:),tmpt(:)
    complex(8)                          :: foo,tail,fg,xi
    real(8)                             :: tau,beta,mues,wn
    L=size(wm) ; N=size(tm)
    NN=min(32*L,2**14)
    xi=cmplx(0.d0,1.d0,8)
    allocate(tmpt(NN),tmpGt(NN))
    tmpt = linspace(0.d0,beta,NN)
    call cubic_spline(gt,tm,tmpGt,tmpt)
    do j=1,L                    !for all Matsubara frequencies:
       wn    = wm(j)!pi/beta*real(2*j-1,8)
       gw(j) = trapz(tmpt(:),exp(xi*wn*tmpt(:))*tmpGt(:))
    enddo
    deallocate(tmpt,tmpGt)
  end subroutine fftgf_tau2iw_upm

end program hmipt_matsuara_2dsquare
