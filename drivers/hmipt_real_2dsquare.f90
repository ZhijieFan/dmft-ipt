program hmipt_2dsquare
  !########################################################
  !     Program  : HMIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model using DMFT-IPT
  !     AUTHORS  : Adriano Amaricci
  !########################################################
  !LOCAL:
  USE DMFT_IPT
  USE VECTORS
  USE SQUARE_LATTICE
  USE IOTOOLS
  USE INTEGRATE
  USE TIMER

  implicit none

  logical                :: converged
  real(8)                :: n,wband,kint,eint,x
  integer                :: i,Lk,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)
  real(8),allocatable    :: wr(:),wt(:),epsik(:),nk(:),dos(:)

  call read_input("inputIPT.in")
  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(wr(L))


  !build freq. array
  wr = linspace(-wmax,wmax,L,mesh=fmesh)

  !build square lattice structure:
  dos=0.d0
  wband=4.d0*ts
  allocate(dos(L))
  do i=1,L
     x=1.d0/wmax*(wr(i)/wband)**2-1.d0
     call comelp(x,kint,eint)
     dos(i)=2.d0/wband/pi**2*kint*heaviside(wband-abs(wr(i)))
  enddo
  dos=dos/sum(dos)/fmesh
  call splot("rho0_realw.ipt",wr,dos)
  dos=dos*fmesh

  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  call get_initial_sigma

  ! call get_optical_conductivity
  ! stop


  !dmft loop:
  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     fg=zero
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i)
        !fg(i) = trapz(fmesh,wt/(zeta-epsik))
        fg(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo
     n   = trapz(fmesh,dimag(fg)*fermi(wr,beta))/trapz(fmesh,dimag(fg))
     fg0 = one/(one/fg + sigma)
     sigma= solve_ipt_sopt(fg0,wr)
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
     call splot("nVSiloop.ipt",iloop,n,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
  call splot("G_realw.ipt",wr,fg,append=printf)
  call splot("G0_realw.ipt",wr,fg0,append=printf)
  call splot("Sigma_realw.ipt",wr,sigma,append=printf)
  ! nk = square_lattice_momentum_distribution(Lk)
  ! call splot("nkVSepsk.ipt",epsik,nk,append=printf)

  call get_optical_conductivity()

contains

  subroutine get_initial_sigma()
    logical :: check
    inquire(file="Sigma_realw.last",exist=check)
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_realw.last",wr,sigma)
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 
       sigma=zero
    endif
  end subroutine get_initial_sigma

  function square_lattice_momentum_distribution(Lk) result(nk)
    integer            :: Lk
    integer,parameter  :: M=1024
    integer            :: ik,i
    type(matsubara_gf) :: gm,sm
    real(8)            :: nk(Lk),wm(M),w
    call allocate_gf(gm,M)
    call allocate_gf(sm,M)
    wm   = pi/beta*real(2*arange(1,M)-1,8)
    call get_matsubara_gf_from_DOS(wr,sigma,sm%iw,beta)
    do ik=1,Lk
       gm%iw=one/(xi*wm - epsik(ik) - sm%iw)
       call fftgf_iw2tau(gm%iw,gm%tau,beta)
       nk(ik)=-gm%tau(M)
    enddo
  end function square_lattice_momentum_distribution


  subroutine get_optical_conductivity()
    integer                               :: ix,iy,iz,ik,i,iw,it,iv
    integer                               :: Nw
    real(8)                               :: Fiw,Fiviw
    type(vect2D)                          :: vk
    real(8),dimension(:),allocatable      :: omega(:)
    real(8),dimension(:),allocatable      :: reoc(:),imoc(:),ock(:)
    complex(8),dimension(:),allocatable   :: oc(:)
    real(8),dimension(:,:),allocatable    :: Ak
    !
    real(8)                               :: t_(0:L/2)
    real(8)                               :: oct(0:L/2)

    Nw=L/2
    allocate(omega(-Nw+1:Nw))
    forall(i=1:L)omega(-Nw+i)=wr(i)

    allocate(Ak(Lk,-Nw+1:Nw))
    forall(ik=1:Lk,i=1:L)Ak(ik,-Nw+i)=-dimag(one/(cmplx(wr(i),eps)-epsik(ik)-sigma(i)))/pi

    allocate(reoc(-Nw+1:Nw),imoc(-Nw+1:Nw),ock(-Nw+1:Nw),oc(-Nw+1:Nw))
    oc   = zero
    reoc = zero
    imoc = zero
    call start_timer
    do ik=1,Lk
       ock = zero
       vk  = square_lattice_velocity(kgrid(ik2ix(ik),ik2iy(ik)))
       do iv=1,Nw
          do iw=-Nw+1,Nw-iv
             Fiw     = fermi(omega(iw),beta)
             Fiviw   = fermi(omega(iv+iw),beta)
             ock(iv) = ock(iv) + ((Fiw-Fiviw)/omega(iv))*Ak(ik,iw)*Ak(ik,iv+iw)*vk%x*vk%y
          enddo                 !end w-loop
       enddo                    !end v-loop
       reoc(:) = reoc(:)+ock(:)*wt(ik)*ts*fmesh
       call eta(ik,Lk)
    enddo
    call stop_timer
    ! forall(i=1:Nw)reoc(-i+1)=-reoc(i)
    ! imoc(1:Nw) = kronig(reoc(1:Nw),omega(1:Nw),Nw)
    ! oc   = cmplx(reoc,imoc,8)
    call splot("OC_realw.ipt",omega,reoc)!,imoc) 

    t_=linspace(0.d0,pi/wmax*Nw,Nw+1)
    oct=zero
    do it=0,Nw
       do iw=1,Nw
          oct(it)=oct(it) + reoc(iw)*exp(xi*omega(iw)*t_(it))*fmesh
       enddo
    enddo
    call splot("OC_t.ipt",t_,oct)

  end subroutine get_optical_conductivity


  subroutine comelp ( hk, ck, ce )
    !*****************************************************************************80
    !
    !! COMELP computes complete elliptic integrals K(k) and E(k).
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
    !    they give permission to incorporate this routine into a user program 
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    07 July 2012
    !
    !  Author:
    !
    !    Shanjie Zhang, Jianming Jin
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) HK, the modulus.  0 <= HK <= 1.
    !
    !    Output, real ( kind = 8 ) CK, CE, the values of K(HK) and E(HK).
    !
    implicit none
    real ( kind = 8 ) ae
    real ( kind = 8 ) ak
    real ( kind = 8 ) be
    real ( kind = 8 ) bk
    real ( kind = 8 ) ce
    real ( kind = 8 ) ck
    real ( kind = 8 ) hk
    real ( kind = 8 ) pk
    pk = 1.0D+00 - hk * hk
    if ( hk == 1.0D+00 ) then
       ck = 1.0D+300
       ce = 1.0D+00
    else
       ak = ((( &
            0.01451196212D+00   * pk &
            + 0.03742563713D+00 ) * pk &
            + 0.03590092383D+00 ) * pk &
            + 0.09666344259D+00 ) * pk &
            + 1.38629436112D+00
       bk = ((( &
            0.00441787012D+00   * pk &
            + 0.03328355346D+00 ) * pk &
            + 0.06880248576D+00 ) * pk &
            + 0.12498593597D+00 ) * pk &
            + 0.5D+00
       ck = ak - bk * log ( pk )
       ae = ((( &
            0.01736506451D+00   * pk &
            + 0.04757383546D+00 ) * pk &
            + 0.0626060122D+00  ) * pk &
            + 0.44325141463D+00 ) * pk &
            + 1.0D+00
       be = ((( &
            0.00526449639D+00   * pk &
            + 0.04069697526D+00 ) * pk &
            + 0.09200180037D+00 ) * pk &
            + 0.2499836831D+00  ) * pk
       ce = ae - be * log ( pk )
    end if
    return
  end subroutine comelp

end program hmipt_2dsquare
