program hmipt_matsubara
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none
  logical                          :: converged,check
  real(8)                          :: wmix,ts
  integer                          :: i,iloop,L,Lk,Nx
  complex(8)                       :: zeta
  complex(8),allocatable           :: fg(:),fg0(:),sigma(:),GFold(:)
  real(8),allocatable              :: wm(:)
  real(8)                          :: n,docc,z,energy(3),h1,h2,hdc
  character(len=24)                :: finput
  real(8),dimension(:),allocatable :: kxgrid,kygrid,Wtk,Epsik

  call parse_cmd_variable(finput,"finput",default="inputIPT.conf")
  call parse_input_variable(ts,"ts",finput,default=0.5d0)
  call parse_input_variable(L,"L",finput,default=4096)
  call parse_input_variable(Nx,"NX",finput,default=21)
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

  !BUILD THE LATTICE STRUCTURE (use tight_binding):
  Lk = Nx*Nx
  allocate(Epsik(Lk),Wtk(Lk))
  allocate(kxgrid(Nx),kygrid(Nx))
  write(*,*) "Using Nk_total="//txtfy(Lk)
  kxgrid = kgrid(Nx)
  kygrid = kgrid(Nx)
  Epsik  = build_hk_model(hk_model,kxgrid,kygrid,[0d0])
  Wtk    = 1d0/Lk
  call write_hk_w90("Hk2d.dat",1,1,0,1,dcmplx(Epsik,0d0),kxgrid,kygrid,[0d0])
  call get_free_dos(Epsik,Wtk)


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
        fg(i) = sum_overk_zeta(zeta,epsik,wtk)
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
  energy = ipt_measure_energy_matsubara(Sigma,fg0,Epsik,Wtk)
  call splot("observables_last.ipt",n,z,docc,energy(1),energy(2),energy(3))

  h1 = sum(Epsik)/Lk
  h2 = sum(Epsik**2)/Lk
  hdc= 0d0
  open(free_unit(i),file="g0_header.ipt")
  write(i,"(4F18.9)")xmu,h1,h2,hdc
  close(i)

  open(free_unit(i),file="sigma_header.ipt")
  write(i,"(2F18.9)")uloc(1),n
  close(i)

contains


  function hk_model(kpoint) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky
    real(8)              :: hk
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = -2d0*ts*(cos(kx)+cos(ky))
  end function hk_model


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
