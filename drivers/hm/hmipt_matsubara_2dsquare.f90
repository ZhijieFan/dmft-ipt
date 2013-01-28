
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
  implicit none

  logical                :: converged
  real(8)                :: n,z
  integer                :: i,Lk
  complex(8)             :: zeta
  type(matsubara_gf)     :: fg,sigma
  complex(8),allocatable :: fg0(:)
  real(8),allocatable    :: wm(:),tau(:),wt(:),epsik(:),nk(:)

  call read_input("inputIPT.in")
  !allocate functions:
  allocate(wm(L),tau(0:L))
  call allocate_gf(fg,L)
  call allocate_gf(sigma,L)
  allocate(fg0(L))

  !build freq. array
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  !dmft loop:
  sigma%iw=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta = xi*wm(i) - sigma%iw(i)
        fg%iw(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo
     call fftgf_iw2tau(fg%iw,fg%tau,beta)
     n   = -fg%tau(L)
     fg0 = one/(one/fg%iw + sigma%iw)
     sigma%iw= solve_ipt_matsubara(fg0)
     converged=check_convergence(sigma%iw,eps_error,nsuccess,nloop)
     z=1.d0 - dimag(sigma%iw(1))/wm(1);z=1.d0/z
     call splot("nVSiloop.ipt",iloop,n,append=TT)
     call splot("zetaVSiloop.ipt",iloop,z,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  call close_file("zetaVSiloop.ipt")
  call splot("G_iw.ipt",wm,fg%iw,append=printf)
  call splot("G0_iw.ipt",wm,fg0,append=printf)
  call splot("Sigma_iw.ipt",wm,sigma%iw,append=printf)
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

end program hmipt_matsuara_2dsquare
