program hmipt_matsubara
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none
  logical                :: converged,check
  real(8)                :: wmix,D
  integer                :: i,iloop,L,unit
  complex(8)             :: zeta
  complex(8),allocatable :: fg(:),fg0(:),sigma(:),GFold(:),delta(:)
  real(8),allocatable    :: wm(:),gtau(:)
  real(8)                :: n,docc,z,energy(3)
  character(len=24)      :: finput

  call parse_cmd_variable(finput,"finput",default="inputIPT.conf")
  call parse_input_variable(D,"wband",finput,default=1d0)
  call parse_input_variable(L,"L",finput,default=4096)
  call parse_input_variable(wmix,"WMIX",finput,default=0.75d0)
  call read_input(finput)

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
     n    = ipt_measure_dens_matsubara(fg)
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
     open(free_unit(unit),file="observables_all.ipt")
     write(unit,*)n,docc,z
     close(unit)
     converged=check_convergence(fg0,dmft_error,nsuccess,nloop)
  enddo

  allocate(delta(L))
  delta = xi*wm+xmu-one/fg0
  call splot("G_iw.ipt",wm,fg)
  call splot("G0_iw.ipt",wm,fg0)
  call splot("Delta_iw.ipt",wm,delta)
  call splot("Sigma_iw.ipt",wm,sigma)


  energy = ipt_measure_energy_matsubara(Sigma,fg0,100,D)
  open(free_unit(unit),file="observables_last.ipt")
  write(unit,"(6F21.12)")n,z,docc,energy(1),energy(2),energy(3)
  close(unit)

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
