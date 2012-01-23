program hmipt
  !########################################################
  !     Program  : HMIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model using DMFT-IPT
  !     AUTHORS  : Adriano Amaricci
  !########################################################
  !LOCAL:
  USE DMFT_IPT
  implicit none

  integer    :: i
  logical    :: converged
  complex(8) :: zeta
  real(8)    :: n
  real(8),allocatable :: wr(:)
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)

  call read_input("inputIPT.in")
  allocate(fg(L),sigma(L),fg0(L),wr(L))

  wr=linspace(-wmax,wmax,L)

  sigma=zero ; iloop=0 ; converged=.false.             
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta  = cmplx(wr(i),eps) - sigma(i)
        fg(i) = gfbether(wr(i),zeta,1.d0)
     enddo
     fg0 = one/(one/fg + sigma)
     sigma= solve_ipt_sopt(fg0,wr)
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
  enddo
  call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
  call splot("G_realw.ipt",wr,fg,append=printf)
  call splot("G0_realw.ipt",wr,fg0,append=printf)
  call splot("Sigma_realw.ipt",wr,sigma,append=printf)

end program hmipt

