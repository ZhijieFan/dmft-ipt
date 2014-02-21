program hmipt_matsuara_dmft_g0
  USE DMFT_IPT
  USE IOTOOLS
  USE ERROR
  implicit none
  logical                :: converged,check
  real(8)                :: n,z
  integer                :: i,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: fg(:),fg0(:),sigma(:),GFold(:)
  real(8),allocatable    :: wm(:),sigt(:)
  real(8) :: C0,C1,n0

  call read_input("inputIPT.in")

  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(sigt(0:L))

  !build freq. array
  allocate(wm(L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)

  !set the energy scale:
  D=2.d0*ts

  !get or read first sigma 
  call  get_initial_weiss_field(fg0,"g0.restart")
  call splot("guessG0_iw.ipt",wm,fg0)

  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !IMPURITY SOLVER
     sigma= solve_ipt_matsubara(fg0)
     !
     sigma=xi*dimag(sigma)!ACTHUNG
     !
     !SELF-CONSISTENCY:
     do i=1,L
        zeta  = xi*wm(i) - sigma(i)
        fg(i) = gfbethe(wm(i),zeta,D)
     enddo
     n   = get_local_density(fg,beta)
     !update weiss-field:
     fg0 = one/(one/fg + sigma)
     !if(iloop>1)fg0 = weight*fg0 + (1.d0-weight)*GFold
     !
     converged=check_convergence(fg0,dmft_error,nsuccess,nloop)
     !GET OBSERVABLES
     z=1.d0 - dimag(sigma(1))/wm(1);z=1.d0/z
     call splot("observables_all.ipt",iloop,beta,xmu,u,z,n,append=.true.)
  enddo
  call splot("G_iw.ipt",wm,fg)
  call splot("G0_iw.ipt",wm,fg0)
  call splot("Sigma_iw.ipt",wm,sigma)
  call splot("observables.ipt",iloop,u,beta,xmu,z,n)

  n0=n/2.d0
  C0=U*(n0-0.5d0)
  C1=U**2*n0*(1.d0-n0)
  print*,n0,C0,C1
  sigma=sigma- C0 - C1/(xi*wm)
  call fftgf_iw2tau(sigma,sigt(0:),beta,notail=.true.)
  sigt=sigt-C1*0.5d0
  open(100,file="fft_sigma_iw.ipt")
  do i=0,L
     write(100,*)i*beta/dble(L),sigt(i)
  enddo
  close(100)

contains

  subroutine get_initial_weiss_field(self,file)
    complex(8),dimension(:)       :: self
    real(8),dimension(size(self)) :: wm_
    character(len=*)              :: file
    logical                       :: check
    inquire(file=file,exist=check)
    if(check)then
       print*,'Reading G0'
       call sread(file,wm_,self)
    else
       print*,"Using non-interacting GF as guess:"
       call bethe_guess_g0(self,D,beta,0.d0)
    endif
  end subroutine get_initial_weiss_field

end program hmipt_matsuara_dmft_g0
