program hmipt
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  logical                :: converged
  real(8)                :: n,dw!,A
  integer                :: i,Lk,iloop
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)
  real(8),allocatable    :: wr(:)

  call read_input("inputIPT.in")
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(wr(L))

  !build freq. array
  wr  = linspace(-wmax,wmax,L,mesh=dw)
  sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !GET GLOC:
     fg=zero
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i)
        fg(i) = gfbether(wr(i),zeta,2.d0*ts)
     enddo
     fg0     = one/(one/fg + sigma)
     sigma   = ipt_solve_keldysh(fg0,wmax)
     converged=check_convergence(sigma,dmft_error,nsuccess,nloop)
  enddo
  call splot("G0_realw.ipt",wr,-dimag(fg0)/pi,dreal(fg0))
  call splot("G_realw.ipt",wr,-dimag(fg)/pi,dreal(fg))
  call splot("Sigma_realw.ipt",wr,sigma)

end program hmipt
