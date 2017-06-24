!     PURPOSE  : Solve the Hubbard model using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
program hmmpt_matsubara
  USE DMFT_IPT
  USE DMFT_TOOLS
  USE SCIFOR
  implicit none
  real(8)                             :: x(1),D,wmix,energy(3)
  logical                             :: check
  complex(8)                          :: zeta
  logical                             :: converged
  integer                             :: i,iloop,L
  real(8),dimension(:),allocatable    :: wm
  real(8)                             :: xmu0,n,n0,z,docc
  complex(8),dimension(:),allocatable :: sigma,fg,fg0,gamma,sold

  call parse_input_variable(L,"L","inputIPT.conf",default=4096)
  call parse_input_variable(D,"wband","inputIPT.conf",default=1d0,comment="half-bandwidth, energy unit")
  call parse_input_variable(wmix,"wmix","inputIPT.conf",default=0.5d0,comment="mixing parameter")
  call read_input("inputIPT.conf")

  allocate(fg(L))
  allocate(fg0(L))
  allocate(sigma(L))
  allocate(gamma(L))
  allocate(sold(L))
  allocate(wm(L))

  wm = pi/beta*(2*arange(1,L)-1)

  call  get_initial_sigma(Sigma,"Sigma_iw.restart")

  xmu0=xmu 



  iloop=0 ; converged=.false.
  do while (.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)")"DMFT-loop",iloop
     do i=1,L
        zeta    = xi*wm(i) + xmu - sigma(i)
        fg(i) = gfbethe(wm(i),zeta,D)
     enddo

     !Get the regular Weiss Field: \Gamma(w+xmu)
     gamma = one/(one/fg + sigma)
     n     = ipt_measure_dens_matsubara(fg)

     !Fix the xmu0 w/ condition n0=n
     x(1)=xmu0
     call broydn(mpt_condition_n_eq_n0,x,check,tolf=1d-3,tolmin=1d-4)
     xmu0=x(1)
     sigma = mpt_solve_matsubara(fg0,n,n0,xmu0)
     if(iloop>1)sigma = wmix*sigma + (1.d0-wmix)*sold
     sold = sigma
     z    = ipt_measure_zeta_matsubara(sigma,fg0)
     docc = ipt_measure_docc_matsubara(sigma,fg0)
     n    = ipt_measure_dens_matsubara(fg)
     call splot("obserbables_all.ipt",n,docc,z,append=.true.)
     converged=check_convergence(sigma,dmft_error,Nsuccess,Nloop)
  enddo
  call splot("G_iw.ipt",wm,fg)
  call splot("G0_iw.ipt",wm,fg0)
  call splot("Sigma_iw.ipt",wm,sigma)
  call splot("obserbables_last.ipt",n,docc,z,energy(1),energy(2),energy(3))

contains


  function mpt_condition_n_eq_n0(x) result(funcv)
    implicit none
    real(8),dimension(:),intent(in)  ::  x
    real(8),dimension(size(x))       ::  funcv
    real(8)                          ::  zn0
    xmu0=x(1)
    !Hartree corrected WF is: \tilde{\calG0}^-1 = \calG0^-1 +xmu -xmu0 -Uloc*n
    fg0 = one/(one/gamma +xmu-xmu0-uloc(1)*(n-0.5d0))
    n0  = ipt_measure_dens_matsubara(fg0)
    funcv(1)=n-n0
    write(*,"(3(f13.9))")n,n0,xmu0
  end function mpt_condition_n_eq_n0


  !Get the initial Sigma: either by reading from file or assuming a HF form
  subroutine get_initial_sigma(self,file)
    complex(8),dimension(:) :: self
    real(8),dimension(size(self)) :: wm
    character(len=*)        :: file
    logical                 :: check
    inquire(file=file,exist=check)
    if(check)then
       print*,'Reading sigma'
       open(100,file=file)
       do i=1,size(self)
          read(100,*)wm(i),self(i)
       enddo
       close(100)
    else
       print*,"Using Hartree-Fock self-energy U*(n-1/2) [n=1/2]"
       print*,"===================================="
       n=0.5d0
       self= Uloc(1)*(n-0.5d0)
    endif
  end subroutine get_initial_sigma



end program hmmpt_matsubara






