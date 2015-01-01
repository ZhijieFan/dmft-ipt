program hmipt_matsubara
  USE DMFT_IPT
  USE CONSTANTS
  USE IOTOOLS
  USE ERROR
  USE ARRAYS
  USE FUNCTIONS
  USE DMFT_TOOLS
  implicit none
  logical                :: converged,check
  real(8)                :: wmix,D
  integer                :: i,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: fg(:),fg0(:),sigma(:),GFold(:)
  real(8),allocatable    :: wm(:)
  real(8) :: n,docc,z,energy(3)

  call parse_input_variable(D,"wband","inputIPT.in",default=1d0)
  call parse_input_variable(wmix,"WMIX","inputIPT.in",default=0.75d0)
  call read_input("inputIPT.in")

  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(GFold(L))

  !build freq. array
  allocate(wm(L))
  wm(:)  = pi/beta*(2*arange(1,L)-1)

  !get or read first sigma 
  call  get_inital_sigma(Sigma,"Sigma_iw.ipt")

  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !SELF-CONSISTENCY:
     do i=1,L
        zeta = xi*wm(i) - sigma(i)
        fg(i) = gfbethe(wm(i),zeta,D)
     enddo
     n   = fft_get_density(fg,beta)
     GFold=fg0
     fg0 = one/(one/fg + sigma)
     if(iloop>1)fg0 = wmix*fg0 + (1.d0-wmix)*GFold
     !
     !IMPURITY SOLVER
     sigma= ipt_solve_matsubara(fg0)

     !GET OBSERVABLES
     n    = ipt_measure_dens_matsubara(sigma,fg0)
     z    = ipt_measure_zeta_matsubara(sigma,fg0)
     docc = ipt_measure_docc_matsubara(sigma,fg0)
     write(*,"(3F15.9,1x)",advance="no")n,docc,z
     call splot("observables_all.ipt",n,docc,z)
     converged=check_convergence(fg0,dmft_error,nsuccess,nloop)
  enddo
  call splot("G_iw.ipt",wm,fg)
  call splot("G0_iw.ipt",wm,fg0)
  call splot("Sigma_iw.ipt",wm,sigma)
  energy = ipt_measure_energy_matsubara(Sigma,fg0,100,D)
  call splot("observables_last.ipt",n,z,docc,energy(1),energy(2),energy(3))
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

end program hmipt_matsubara
