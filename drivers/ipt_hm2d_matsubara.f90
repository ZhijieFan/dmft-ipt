program hmipt_matsubara
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none
  logical                             :: converged,check
  real(8)                             :: wmix,ts
  integer                             :: i,iloop,L,Lk,Ne,unit
  complex(8)                          :: zeta
  complex(8),dimension(:),allocatable :: Gmats,Weiss,Smats,Weiss_prev
  real(8),allocatable                 :: wm(:)
  real(8)                             :: n,docc,z,energy(3),depsi
  character(len=24)                   :: finput
  real(8),dimension(:),allocatable    :: epsi,Dos2d

  call parse_cmd_variable(finput,"finput",default="inputIPT.conf")
  call parse_input_variable(ts,"ts",finput,default=0.5d0)
  call parse_input_variable(L,"L",finput,default=4096)
  call parse_input_variable(Ne,"Ne",finput,default=1000)
  call parse_input_variable(wmix,"WMIX",finput,default=0.75d0)
  call read_input(finput)

  !allocate functions:
  allocate(Gmats(L))
  allocate(Smats(L))
  allocate(Weiss(L))
  allocate(Weiss_prev(L))

  !build freq. array
  allocate(wm(L))
  wm(:)  = pi/beta*(2*arange(1,L)-1)

  write(*,*) "Using N_energies="//txtfy(Ne)
  !
  allocate(epsi(Ne))
  allocate(dos2d(Ne))
  epsi = linspace(-wmax,wmax,Ne,mesh=depsi)
  do i=1,Ne
     Dos2d(i)=dens_2dsquare(epsi(i),ts)*depsi
  enddo


  !get or read first sigma
  Smats = zero
  inquire(file="Smats_iw.ipt",exist=check)
  if(check)then
     print*,'Reading sigma'
     call sread("Smats_iw.ipt",wm,Smats)
  endif

  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !
     !SELF-CONSISTENCY:
     do i=1,L
        zeta  = xi*wm(i) - Smats(i)
        Gmats(i) = sum_overk_zeta(zeta,epsi,dos2d)
     enddo
     !
     !< Self-consistency
     Weiss  = one/(one/Gmats + Smats)
     !
     call splot("Gmats_iw.ipt",wm,Gmats)
     call splot("Smats_iw.ipt",wm,Smats)
     call splot("Weiss_iw.ipt",wm,Weiss)
     !
     if(iloop>1)Weiss = wmix*Weiss + (1.d0-wmix)*Weiss_prev
     Weiss_prev=Weiss
     !
     !IMPURITY SOLVER
     Smats= ipt_solve_matsubara(Weiss)

     !GET OBSERVABLES
     n    = ipt_measure_dens_matsubara(Smats,Gmats)
     z    = ipt_measure_zeta_matsubara(Smats,Gmats)
     docc = ipt_measure_docc_matsubara(Smats,Gmats)
     write(*,"(3F15.9,1x)",advance="no")n,docc,z
     !
     converged=check_convergence(Weiss,dmft_error,nsuccess,nloop)
  enddo


  call save_array("g0.restart",Weiss)

  energy = ipt_measure_energy_matsubara(Smats,Gmats,epsi,dos2d)
  open(free_unit(unit),file="observables_last.ipt")
  write(unit,*)n,z,docc,energy(1),energy(2),energy(3)
  close(unit)


end program hmipt_matsubara
