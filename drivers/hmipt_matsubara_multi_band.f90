program hmipt_matsuara_multi_band
  Use DMFT_IPT
  USE IOTOOLS
  USE ERROR
  implicit none
  logical                :: converged,check

  integer                :: i,iloop
  integer                :: Nband
  complex(8),allocatable :: zeta(:),fg(:,:),fg0(:,:),sigma(:,:)
  real(8),allocatable    :: wm(:),n(:),z(:)

  call read_input("inputIPT.in")
  call parse_cmd_variable(nband,"NBAND",default=2)

  !allocate functions:
  allocate(fg(Nband,L))
  allocate(sigma(Nband,L))
  allocate(fg0(Nband,L))
  allocate(zeta(2),n(2),z(2))

  !build freq. array
  allocate(wm(L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)

  !get or read first sigma 
  sigma=zero !U_a*(n_a-1/2)

  !dmft loop:
  D=2.d0*ts ;  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !SELF-CONSISTENCY:
     do i=1,L
        zeta(1) = xi*wm(i) - sigma(1,i)
        zeta(2) = xi*wm(i) - sigma(2,i)
        fg(1,i) = gfbethe(wm(i),zeta(1),D)
        fg(2,i) = gfbethe(wm(i),zeta(2),D)
     enddo
     n(1)   = get_local_density(fg(1,:),beta)
     n(2)   = get_local_density(fg(2,:),beta)
     fg0(1,:) = one/(one/fg(1,:) + sigma(1,:))
     fg0(2,:) = one/(one/fg(2,:) + sigma(2,:))
     !
     !IMPURITY SOLVER
     sigma = solve_ipt_matsubara(fg0)
     converged=check_convergence(fg0(1,:)+fg0(2,:),dmft_error,nsuccess,nloop)
     !GET OBSERVABLES
     z(1)=1.d0 - dimag(sigma(1,1))/wm(1);z(1)=1.d0/z(1)
     z(2)=1.d0 - dimag(sigma(2,1))/wm(1);z(2)=1.d0/z(2)
     call splot("observables_all.ipt",dble(iloop),u,n(1),n(2),z(1),z(2),append=.true.)
  enddo
  call splot("G_iw.ipt",wm,fg(1,:),fg(2,:))
  call splot("G0_iw.ipt",wm,fg0(1,:),fg0(2,:))
  call splot("Sigma_iw.ipt",wm,sigma(1,:),sigma(2,:))
  call splot("observables.ipt",u,n(1),n(2),z(1),z(2))


  ! call fftgf_iw2tau(sigma,sigt(0:),beta,notail=.true.)
  ! open(100,file="fft_sigma_iw.ipt")
  ! do i=0,L
  !    write(100,*)i*beta/dble(L),sigt(i)
  ! enddo
  ! close(100)

contains

  subroutine get_inital_sigma(self,file)
    complex(8),dimension(:,:)       :: self
    real(8),dimension(size(self,2)) :: wm
    character(len=*)                :: file
    logical                         :: check
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

end program hmipt_matsuara_multi_band
