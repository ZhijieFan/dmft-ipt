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
  USE IOTOOLS
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

end program hmipt
