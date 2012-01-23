!########################################################
!     Program  : AHMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the attractive Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
module COMMON
  USE DMFT_IPT
  USE LATTICE
  implicit none
  complex(8),allocatable :: sigma(:,:),fg(:,:),fg0(:,:),calG(:,:),sold(:,:)
  real(8)                :: n,delta,n0,delta0,ntarget

contains

  subroutine dmft_solve_search_mu(nread,niter,nerror)    
    real(8),save :: ndelta=0.1d0
    real(8)      :: nread,nerror
    real(8) :: ndelta1,ntmp
    integer :: niter,nindex,nindex1,iteration
    interface 
       function get_dens(x)
         real(8),intent(in) :: x
         real(8)            :: get_dens
       end function get_dens
    end interface

    nindex=0    
    print*,"Search mu within:",nerror,niter
    iteration=0
    do while(iteration < niter)
       iteration=iteration+1
       !====================
       ntmp= get_dens(xmu)
       !====================
       nindex1=nindex
       ndelta1=ndelta
       if((ntmp.ge.nread+nerror))then
          nindex=-1
       elseif(ntmp.le.nread-nerror)then
          nindex=1
       else
          nindex=0
          niter=1
       endif
       if(nindex1+nindex.eq.0)then !avoid loop forth and back
          ndelta=ndelta1/2.d0 !decreasing the step
          xmu=xmu+dble(nindex)*ndelta
       else
          ndelta=ndelta1
          xmu=xmu+dble(nindex)*ndelta
       endif
       write(*,"(A,f15.12,A,f15.12,A,f15.12)")" n=",ntmp,"/",nread,"| ",ndelta
       call splot("muVSiter.qmc",iteration,xmu,append=.true.)
    enddo
    print*,'----------------------'
    print*,'#iteration=',iteration
    print*,'        mu=',xmu
    print*,'         n=',ntmp,'/',nread
    print*,'     error=',abs(ntmp-nread)
    print*,''
    return
  end subroutine dmft_solve_search_mu
end module COMMON



program ahmmpt_fixn
  USE COMMON
  implicit none
  real(8)    :: x
  integer                :: ik
  complex(8)             :: zeta1,zeta2,det
  call read_input("inputIPT.in")
  ntarget=xmu
  xmu=0.d0
  allocate(fg(2,-L:L),sigma(2,-L:L),fg0(2,-L:L),calG(2,-L:L))
  allocate(sold(2,-L:L))

  call  build_BetheLattice(Nx,D,Nk=Lk)
  allocate(epsik(Lk)) ; call get_epsik_bethe(epsik,D)

  n=0.5d0 ;  delta=deltasc
  sigma=zero ; sigma(2,:)=-delta ; sold=sigma

  call dmft_solve_search_mu(ntarget,nloop,eps_error)

end program ahmmpt_fixn





function get_dens(x)
  USE COMMON
  implicit none
  real(8),intent(in) :: x
  real(8)            :: get_dens
  integer                :: ik
  complex(8)             :: zeta1,zeta2,det

  xmu=x

  fg=zero
  do i=-L,L
     zeta1 = cmplx(wr(i),eps,8) + xmu - sigma(1,i)
     zeta2 = cmplx(wr(i),eps,8) - xmu + sigma(1,-i)
     do ik=1,Lk
        det     = (zeta1-epsik(ik))*(zeta2+epsik(ik)) - sigma(2,i)*sigma(2,i)
        fg(1,i) =fg(1,i) + wt(ik)*(zeta2+epsik(ik))/det
        fg(2,i) =fg(2,i) + wt(ik)*sigma(2,i)/det
     enddo
  enddo
  n    =  sum(aimag(fg(1,:))*fermi(wr,beta))/sum(dimag(fg(1,:))) !*fmesh/pi
  delta = -u*sum(aimag(fg(2,:))*fermi(wr,beta))*fmesh/pi

  do i=-L,L
     det      = fg(1,i)*fg(1,-i) - fg(2,i)*(fg(2,i))
     calG(1,i)= fg(1,-i)/det + sigma(1,i) - u*(n-0.5d0)
     calG(2,i)= fg(2,i)/det  + sigma(2,i) - delta
  end do
  do i=-L,L
     det     =  calG(1,i)*calG(1,-i) - calG(2,i)*calG(2,i)
     fg0(1,i)=  calG(1,-i)/det
     fg0(2,i)=  calG(2,i)/det
  end do
  n0    =  sum(aimag(fg0(1,:))*fermi(wr,beta))/sum(dimag(fg0(1,:))) !*fmesh/pi
  delta0= -u*sum(aimag(fg0(2,:))*fermi(wr,beta))*fmesh/pi
  sigma =  solve_mpt_sc_sopt(fg0,fg,n,n0,delta,delta0)
  !sigma = weigth*sigma + (1.d0-weigth)*sold ; sold=sigma

  get_dens=n
end function get_dens
