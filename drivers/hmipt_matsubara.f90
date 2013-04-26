program hmipt_matsuara
  USE DMFT_IPT
  USE IOTOOLS
  implicit none
  logical                :: converged
  real(8)                :: n,z
  integer                :: i,iloop
  complex(8)             :: zeta
  type(matsubara_gf)     :: fg
  complex(8),allocatable :: fg0(:),sigma(:),GFold(:)
  real(8),allocatable    :: wm(:)

  call read_input("inputIPT.in")
  !allocate functions:
  allocate(wm(L))
  call allocate_gf(fg,L)
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(GFold(L))

  !build freq. array
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)

  !dmft loop:
  D=2.d0*ts 
  sigma=zero!one/(xi*wm)
  inquire(file="Sigma.restart",exist=converged)
  if(converged)then
     print*,'reading sigma'
     call sread("Sigma.restart",sigma)
  end if
  GFold=sigma
  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta = xi*wm(i) - sigma(i)
        fg%iw(i) = gfbethe(wm(i),zeta,D)
     enddo
     call fftgf_iw2tau(fg%iw,fg%tau,beta)
     n   = -fg%tau(L)
     GFold=fg0
     fg0 = one/(one/fg%iw + sigma)
     if(iloop>1)fg0 = weight*fg0 + (1.d0-weight)*GFold
     sigma= solve_ipt_matsubara(fg0)
     converged=check_convergence(fg0,eps_error,nsuccess,nloop)
     z=1.d0 - dimag(sigma(1))/wm(1);z=1.d0/z
     call splot("observables.ipt",dble(iloop),u,z,beta,xmu,append=.true.)
  enddo
  call splot("G_iw.ipt",wm,fg%iw,append=printf)
  call splot("G0_iw.ipt",wm,fg0,append=printf)
  call splot("Sigma_iw.ipt",wm,sigma,append=printf)
  call splot("Sigma.restart",sigma)
  call splot("observables_last.ipt",u,z,beta,xmu,append=printf)
end program hmipt_matsuara
