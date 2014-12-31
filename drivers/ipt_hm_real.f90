program hmipt
  USE DMFT_IPT
  USE ERROR
  USE ARRAYS
  USE FUNCTIONS
  USE IOTOOLS
  implicit none
  integer                :: i,iloop
  logical                :: converged
  complex(8)             :: zeta
  real(8)                :: n
  real(8),allocatable    :: wr(:)
  complex(8),allocatable :: sigma(:),fg(:),fg0(:),sold(:)

  call read_input("inputIPT.in")
  allocate(fg(L),sigma(L),fg0(L),wr(L),sold(L))

  wr=linspace(-wmax,wmax,L)

  sigma=zero ; iloop=0 ; converged=.false.       
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta  = cmplx(wr(i),eps) - sigma(i)
        fg(i) = gfbether(wr(i),zeta,2.d0*ts)
     enddo
     fg0 = one/(one/fg + sigma)
     sold = sigma
     sigma= ipt_solve_real(fg0,wr)
     sigma = wmix*sigma + (1.d0-wmix)*sold
     converged=check_convergence(sigma,dmft_error,nsuccess,nloop)
  enddo
  call splot("G_realw.ipt",wr,-dimag(fg)/pi,dreal(fg))
  call splot("G0_realw.ipt",wr,fg0)
  call splot("Sigma_realw.ipt",wr,sigma)

end program hmipt

