program hmipt_matsuara
  USE DMFT_IPT
  USE IOTOOLS
  USE ERROR
  implicit none
  logical                :: converged,check
  real(8)                :: n,z,de,e
  integer                :: i,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: fg(:),fg0(:),sigma(:)
  real(8),allocatable    :: wm(:),epsi(:),dos(:)

  call read_input("inputIPT.in")

  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))

  !build freq. array
  allocate(wm(L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)

  allocate(dos(Nx),epsi(Nx))
  D=sqrt(2.d0)*ts
  epsi = linspace(-5.d0*D,5.d0*D,Nx,mesh=de)
  do i=1,Nx
     dos(i) = dens_hyperc(epsi(i),D)
  enddo
  call splot("DOShyperc.ipt",epsi,dos)
  print*,trapz(de,dos)

  !get or read first sigma 
  call  get_inital_sigma(Sigma,"Sigma.restart")

  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !SELF-CONSISTENCY:
     fg=zero
     do i=1,L
        zeta = xi*wm(i) - sigma(i)
        fg(i) = sum_overk_zeta(zeta,epsi,dos)*de
     enddo
     n   = get_local_density(fg,beta)
     fg0 = one/(one/fg + sigma)
     !
     !IMPURITY SOLVER
     sigma= solve_ipt_matsubara(fg0)
     sigma=xi*dimag(sigma)
     converged=check_convergence(fg0,dmft_error,nsuccess,nloop)
     !GET OBSERVABLES
     z=1.d0 - dimag(sigma(1))/wm(1);z=1.d0/z
     call splot("observables_all.ipt",dble(iloop),u,z,beta,append=.true.)
  enddo
  call splot("G_iw.ipt",wm,fg)
  call splot("G0_iw.ipt",wm,fg0)
  call splot("Sigma_iw.ipt",wm,sigma)
  call splot("observables.ipt",u,beta,n,z)

contains

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

end program hmipt_matsuara
