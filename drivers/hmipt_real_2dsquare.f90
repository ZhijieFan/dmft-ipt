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
  real(8)                :: n
  integer                :: i,Lk,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)
  real(8),allocatable    :: wr(:),wt(:),epsik(:),nk(:)

  call read_input("inputIPT.in")
  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(wr(L))


  ! call get_optical_conductivity
  ! stop

  !build freq. array
  wr = linspace(-wmax,wmax,L,mesh=fmesh)

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  call get_initial_sigma
  !dmft loop:
  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i)
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
  nk = square_lattice_momentum_distribution(Lk)
  call splot("nkVSepsk.ipt",epsik,nk,append=printf)
  call get_optical_conductivity

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
    real(8) :: scond(L/2),ome(L/2),dw
    real(8) :: Ak(Lk,L),FF,AA
    integer :: ix,iy,iz,ik,i,iw
    type(vect2D) :: vk
    integer :: x_p_y(L,L)
    x_p_y=0
    do ix=1,L
       do iy=1,L
          iz = ix+iy!-1!+L/2
          if(iz<1 .OR. iz>L) iz=-1
          x_p_y(ix,iy)=iz
       enddo
    enddo
    forall(ik=1:Lk,i=1:L)Ak(ik,i)=-dimag(one/(cmplx(wr(i),eps)-epsik(ik)-sigma(i)))/pi
    ome=linspace(5.d-2,wmax,L/2)
    call start_timer
    scond=0.d0
    do iw=1,L/2
       do ik=1,Lk
          vk=square_lattice_velocity(kgrid(ik2ix(ik),ik2iy(ik)))
          do ix=1,L-iw+1
             ! iz=x_p_y(iw,ix)
             ! if(iz>0)then
             iz=iw+ix
             AA=Ak(ik,ix)*Ak(ik,iz)
             FF=(fermi(wr(ix),beta)-fermi(wr(iz),beta))
             scond(iw)=scond(iw) + AA*FF*fmesh*wt(ik)*vk%x*vk%y
             ! endif
          enddo
       enddo
       scond(iw)=scond(iw)*pi/ts/ome(iw)
       call eta(iw,L/2)
    enddo
    call stop_timer
    call splot("OC_realw.ipt",ome,scond)
  end subroutine get_optical_conductivity

end program hmipt_2dsquare
