
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
  USE ERROR
  implicit none

  logical                :: converged
  real(8)                :: n,z,epot,doble
  integer                :: i,Lk,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: fg(:),fg0(:),sigma(:)
  real(8),allocatable    :: wm(:),wt(:),epsik(:),nk(:),sigma_tau(:),SxG(:),fg_tau(:)

  call read_input("inputIPT.in")
  !allocate functions:
  allocate(fg(L),sigma(L),fg0(L))

  !build freq. array
  allocate(wm(L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx,Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  !get or read first sigma 
  call  get_inital_sigma(Sigma,"Sigma.restart")

  !dmft loop:
  iloop=0 ; converged=.false.
  do while (.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !SELF-CONSISTENCY:
     do i=1,L
        zeta = xi*wm(i) - sigma(i)
        fg(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo
     n   = get_local_density(fg,beta)
     fg0 = one/(one/fg + sigma)
     !
     !IMPURITY SOLVER
     sigma=solve_ipt_matsubara(fg0)
     sigma=xi*dimag(sigma)
     converged=check_convergence(fg0,dmft_error,nsuccess,nloop)
     !GET OBSERVABLES
     z=1.d0 - dimag(sigma(1))/wm(1);z=1.d0/z
     call splot("observables_all.ipt",dble(iloop),u,beta,n,z,append=.true.)
  enddo
  call splot("G_iw.ipt",wm,fg)
  call splot("G0_iw.ipt",wm,fg0)
  call splot("Sigma_iw.ipt",wm,sigma)
  call splot("observables.ipt",u,beta,n,z)


  allocate(sigma_tau(0:L),fg_tau(0:L))
  open(100,file="Sigma_tau.ipt")
  do i=0,L
     read(100,*)z,sigma_tau(i)
  enddo
  close(100)

  allocate(SxG(0:L))
  call fftgf_iw2tau(fg,fg_tau(0:),beta)
  do i=0,L
     SxG(i) = sigma_tau(L-i)*fg(i)*beta/dble(L)
  enddo

  Epot=0.d0
  do i=1,L
     Epot=Epot+dreal(Sigma(i)*fg(i))
  enddo
  Epot=2.d0*Epot/beta
  doble=0.5d0*n - 0.25d0
  if(u > 0.01d0)doble=Epot/U + 0.5d0*n - 0.25d0
  print*,doble,Epot

  nk = square_lattice_momentum_distribution(Lk)
  call splot("nkVSepsk.ipt",epsik,nk)
contains

  function square_lattice_momentum_distribution(Lk) result(nk)
    integer                :: Lk
    integer                :: ik,i
    real(8)                :: nk(Lk)
    do ik=1,Lk
       fg=one/(xi*wm - epsik(ik) - sigma)
       nk(ik)= get_local_density(fg,beta)
    enddo
  end function square_lattice_momentum_distribution


  subroutine get_inital_sigma(self,file)
    complex(8),dimension(:) :: self
    real(8),dimension(size(self)) :: wm
    character(len=*)        :: file
    logical                 :: check
    inquire(file=file,exist=check)
    if(check)then
       print*,'Reading sigma'
       call sread(file,wm,self)
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       self=zero !U*(n-1/2)
    endif
  end subroutine get_inital_sigma

end program hmipt_matsuara_2dsquare
