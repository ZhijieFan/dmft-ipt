program hmipt_2dsquare
  !########################################################
  !     Program  : HMIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model using DMFT-IPT
  !     AUTHORS  : Adriano Amaricci
  !########################################################
  !LOCAL:
  USE DMFT_IPT
  USE SQUARE_LATTICE
  USE IOTOOLS
  implicit none

  logical                :: converged
  real(8)                :: n
  integer                :: i,Lk
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)
  real(8),allocatable    :: wr(:),wt(:),epsik(:),nk(:)

  call read_input("inputIPT.in")
  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(wr(L))


  !build freq. array
  wr = linspace(-wmax,wmax,L,mesh=fmesh)

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  !dmft loop:
  sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i)
        fg(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo
     n   = sum(aimag(fg)*fermi(wr,beta))/sum(aimag(fg))
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

contains

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

end program hmipt_2dsquare
