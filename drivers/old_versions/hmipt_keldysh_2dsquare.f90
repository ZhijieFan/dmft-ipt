program hmipt
  !########################################################
  !     Program  : HMIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model using DMFT-IPT
  !     AUTHORS  : Adriano Amaricci
  !########################################################
  !LOCAL:
  USE DMFT_IPT
  USE SQUARE_LATTICE
  USE VECTORS
  USE IOTOOLS
  USE TIMER
  implicit none

  logical                :: converged
  real(8)                :: n,dw
  integer                :: i,Lk,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:),fg0less(:),fg0gtr(:)
  real(8),allocatable    :: wr(:),wt(:),epsik(:),nk(:),t(:)

  call read_input("inputIPT.in")
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0less(L),fg0gtr(L),fg0(L))
  allocate(wr(L))

  !build freq. array
  wr    = linspace(-wmax,wmax,L,mesh=dw)

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)


  sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     !GET GLOC:
     fg=zero
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i)
        fg(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo

     !Update Keldysh components of the Weiss Field:
     fg0     = one/(one/fg + sigma)
     sigma   = solve_ipt_keldysh(fg0,wmax)
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
     call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
     call splot("G0_realw.ipt",wr,fg0,append=printf)
     call splot("Sigma_realw.ipt",wr,sigma,append=printf)
  enddo

  call get_optical_conductivity
  call get_oc2

contains

  subroutine get_optical_conductivity()
    integer                               :: ix,iy,iz,ik,i,j,iw,it,iv
    integer                               :: Nw
    real(8)                               :: Fiw,Fiviw
    type(vect2D)                          :: vk,kt
    real(8),dimension(:),allocatable      :: omega(:)
    real(8),dimension(:),allocatable      :: reoc(:),imoc(:),ock(:)
    complex(8),dimension(:),allocatable   :: oc(:)
    real(8),dimension(:,:),allocatable    :: Ak
    !
    real(8)                               :: t_(0:L/2)
    real(8)                               :: oct(0:L/2)
    real(8)      :: sx,sy,ex,ey
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
       ix  = ik2ix(ik)
       iy  = ik2iy(ik)
       kt = kgrid(ix,iy)
       sx = 2.d0*sin(kt%x)
       sy = 2.d0*sin(kt%y)
       !vk  = square_lattice_velocity(ix,iy)
       do iv=1,Nw
          do iw=-Nw+1,Nw-iv
             Fiw     = fermi(omega(iw),beta)
             Fiviw   = fermi(omega(iv+iw),beta)
             ock(iv) = ock(iv) + ((Fiw-Fiviw)/omega(iv))*Ak(ik,iw)*Ak(ik,iv+iw)*sx*sy!vk%x*vk%y
          enddo                 !end w-loop
       enddo                    !end v-loop
       reoc(:) = reoc(:)+ock(:)*wt(ik)*ts*dw

       call eta(ik,Lk)
    enddo
    call stop_timer
    ! forall(i=1:Nw)reoc(-i+1)=-reoc(i)
    ! imoc(1:Nw) = kronig(reoc(1:Nw),omega(1:Nw),Nw)
    ! oc   = cmplx(reoc,imoc,8)
    call splot("OC_realw.ipt",omega(1:Nw),reoc(1:Nw))!,imoc) 

    t_=linspace(0.d0,pi/wmax*Nw,Nw+1)
    oct=zero
    do it=0,Nw
       do iw=1,Nw
          oct(it)=oct(it) + reoc(iw)*exp(xi*omega(iw)*t_(it))*dw
       enddo
    enddo
    call splot("OC_t.ipt",t_,oct)
  end subroutine get_optical_conductivity





  subroutine get_oc2
    integer                      :: ix,iy,iz,ik,i,j,iw,it,is
    integer                      :: Nw
    real(8)                      :: Fiw,Fiviw,A,dt
    type(vect2D)                 :: kt
    real(8),allocatable          :: omega(:),t(:)
    real(8),allocatable          :: reoc(:),imoc(:),ock(:)
    real(8),allocatable          :: ocij(:,:),chiij(:,:)
    real(8),allocatable          :: oc(:)
    complex(8),allocatable       :: ocw(:)
    type(keldysh_equilibrium_gf) :: Gk
    real(8)                      :: sx,sy,ex,ey

    Nw=L/2

    dt       = pi/wmax
    allocate(t(-Nw:Nw))
    t   = linspace(-dt*real(Nw,8),dt*real(Nw,8),2*Nw+1,mesh=dt)

    allocate(omega(-Nw+1:Nw))
    forall(i=1:L)omega(-Nw+i)=wr(i)

    call allocate_gf(Gk,Nw)

    allocate(ock(-Nw:Nw),oc(-Nw:Nw))
    allocate(ocij(0:Nw,0:Nw),chiij(0:Nw,0:Nw))
    ock=0.d0
    oc =0.d0
    ocij=zero
    chiij=zero
    call start_timer
    do ik=1,Lk
       ix=ik2ix(ik)
       iy=ik2iy(ik)
       kt = kgrid(ix,iy)
       sx = 2.d0*sin(kt%x)
       sy = 2.d0*sin(kt%y)
       ex = 2.d0*cos(kt%x)
       ey = 2.d0*cos(kt%y)
       !Get Gk^x=R,<,>(t)
       Gk%ret%w = one/(cmplx(wr(:),eps)-epsik(ik)-sigma(:))
       do i=1,2*Nw
          A            = -dimag(Gk%ret%w(i))/pi
          Gk%less%w(i) = pi2*xi*fermi(wr(i),beta)*A
          Gk%gtr%w(i)  = pi2*xi*(fermi(wr(i),beta)-1.d0)*A
       enddo
       call fftgf_rw2rt(Gk%less%w,Gk%less%t,Nw) ; Gk%less%t=dw/pi2*Gk%less%t
       call fftgf_rw2rt(Gk%gtr%w,Gk%gtr%t,Nw)   ; Gk%gtr%t =dw/pi2*Gk%gtr%t
       call fftgf_rw2rt(Gk%ret%w,Gk%ret%t,Nw)   ; Gk%ret%t =dw/pi2*Gk%ret%t
       do i=-Nw,Nw
          ock(i)=-2.d0/pi*sx**2*dimag(Gk%ret%t(i)*Gk%less%t(-i))!+2.d0/pi*ex*xi*Gk%less%t(i)
       enddo
       call cfft_1d_ex(one*ock)
       do i=0,Nw
          do j=0,Nw
             chiij(i,j) = chiij(i,j) + wt(ik)*ock(i-j)
          enddo
       enddo
       call eta(ik,Lk)
    enddo
    call stop_timer

    !Get OC from Chi:
    do i=0,Nw
       do j=0,Nw
          do it=j,Nw
             ocij(i,j)=ocij(i,j)-chiij(i,it)*dt
          enddo
       enddo
    enddo

    forall(i=0:Nw,j=0:Nw)oc(i-j)=ocij(i,j)
    call splot("OC2_t.ipt",t(0:),oc(0:))

    allocate(ocw(Nw))
    ocw=0.d0
    do iw=1,Nw
       !do it=0,Nw
       it=Nw
       do is=0,it
          ocw(iw)=ocw(iw) + ocij(it,it-is)*exp(-xi*omega(iw)*t(is))*dt
       enddo
       !enddo
    enddo
    call splot("OC2_realw.ipt",omega(1:Nw),ocw(1:Nw))

  end subroutine get_oc2


end program hmipt
