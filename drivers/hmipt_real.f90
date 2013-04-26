program hmipt
  USE DMFT_IPT
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
     sigma= solve_ipt_sopt(fg0,wr)
     sigma = weight*sigma + (1.d0-weight)*sold
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
  enddo
  call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
  call splot("G_realw.ipt",wr,fg,append=printf)
  call splot("G0_realw.ipt",wr,fg0,append=printf)
  call splot("Sigma_realw.ipt",wr,sigma,append=printf)

end program hmipt

