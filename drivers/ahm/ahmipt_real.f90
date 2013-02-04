!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt
  USE DMFT_IPT 
  USE IOTOOLS
  implicit none
  integer    :: i,ik,Lk,iloop
  logical    :: converged
  complex(8) :: det,zeta1,zeta2
  real(8)    :: delta,n
  !
  complex(8),allocatable          :: sigma(:,:),fg(:,:),wf0(:,:),calG(:,:),sold(:,:)
  complex(8),allocatable          :: zeta(:)
  type(matsubara_gf),dimension(2) :: gf,sf
  !
  real(8),allocatable :: wt(:),epsik(:),wr(:),wm(:),tau(:)

  include "revision.inc"
  call version(revision)
  call read_input("inputIPT.in")

  allocate(fg(2,L),sigma(2,L))
  allocate(wf0(2,L),calG(2,L))
  allocate(zeta(L),sold(2,L))

  allocate(wr(L))
  wr = linspace(-wmax,wmax,L,mesh=fmesh)

  Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,2.d0*ts)

  call get_initial_sigma


  iloop=0    ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !
     ! G(iw) F(iw)      --> G(iw) F(iw)
     ! F*(-iw) -G(-iw)  --> F(iw) -G*(iw)
     !
     !analytic continuation: iw--> w+ie
     ! G(w+ie) F(w+ie)     --> G(w) F(w)
     ! F*(-w-ie) -G(-w-ie) --> F(-w) -G*(-w)
     !
     ! w+ie-Sigma(w)-e(k); -S(w)        --> zeta(w)-e(k);-S(w)
     ! -S(-w); -(-w-ie-Sigma*(-w)-e(k)) --> -S(-w); -(zeta*(-w)-e(k))
     !
     ! G = (zeta*(-w)-e(k))/{(zeta(w)-e(k))(zeta*(-w)-e(k))+S(w)S(-w)}
     ! F = -S(-w)/{(zeta(w)-e(k))(zeta*(-w)+e(k))+S(w)S(-w)}
     fg=zero
     zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
     do i=1,L
        zeta1 = zeta(i)
        zeta2 = conjg(zeta(L+1-i))
        do ik=1,Lk
           det = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + conjg(sigma(2,L+1-i))*sigma(2,i)
           fg(1,i)=fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det
           fg(2,i)=fg(2,i) - wt(ik)*conjg(sigma(2,L+1-i))/det
        enddo
     enddo
     delta=-u*sum(dimag(fg(2,:))*fermi(wr,beta))*fmesh/pi
     n    =-sum(dimag(fg(1,:))*fermi(wr,beta))*fmesh/pi


     !get the Weiss Field \calG0^-1(w)
     !
     ! wf(1)=calG0^-1 = G*(-w)/(G(w)G*(-w) + F(w)F(-w)) + Sigma(w)
     ! wf(2)=calF0^-1 =-F(-w) /(G(w)G*(-w) + F(w)F(-w)) + S(w)
     do i=1,L
        det     = fg(1,i)*conjg(fg(1,L+1-i)) + conjg(fg(2,L+1-i))*fg(2,i)
        wf0(1,i)= conjg(fg(1,L+1-i))/det  + sigma(1,i)   +  u*(n-0.5d0)
        wf0(2,i)= fg(2,i)/det       + conjg(sigma(2,L+1-i))   + delta
     end do
     do i=1,L
        det      =  wf0(1,i)*conjg(wf0(1,L+1-i)) + conjg(wf0(2,L+1-i))*wf0(2,i)
        calG(1,i)=  conjg(wf0(1,L+1-i))/det
        calG(2,i)=  conjg(wf0(2,L+1-i))/det
     end do

     write(*,"(2f14.9)",advance="no")2.d0*n,delta
     sold=sigma
     sigma =  solve_ipt_sc_sopt(calG(1:2,:),wr,delta,L)
     sigma = weight*sigma + (1.d0-weight)*sold
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=eps_error,N1=Nsuccess,N2=Nloop)

     call splot("DOS.ipt",wr,-dimag(fg(1,:))/pi,append=printf)
     call splot("G_realw.ipt",wr,fg(1,:),append=printf)
     call splot("F_realw.ipt",wr,fg(2,:),append=printf)
     call splot("Sigma_realw.ipt",wr,sigma(1,:),append=printf)
     call splot("Self_realw.ipt",wr,sigma(2,:),append=printf)
     call splot("calG0_realw.ipt",wr,calG(1,:),append=printf)
     call splot("calF0_realw.ipt",wr,calG(2,:),append=printf)
     call splot("observables.ipt",xmu,u,n,delta,beta,dble(iloop),append=printf)
  enddo

  call splot("observables.last",xmu,u,n,delta,beta,dble(iloop),append=printf)
  call splot("DOS.last",wr,-dimag(fg(1,:))/pi)
  call splot("G_realw.last",wr,fg(1,:))
  call splot("F_realw.last",wr,fg(2,:))
  call splot("Sigma_realw.last",wr,sigma(1,:))
  call splot("Self_realw.last",wr,sigma(2,:))
  call splot("calG0_realw.last",wr,calG(1,:))
  call splot("calF0_realw.last",wr,calG(2,:))

contains 

  subroutine search_mu(convergence)
    integer, save         ::nindex
    integer               ::nindex1
    real(8)               :: naverage,ndelta1
    logical,intent(inout) :: convergence
    naverage=n
    nindex1=nindex
    ndelta1=ndelta
    if((naverage >= nread+nerror))then
       nindex=-1
    elseif(naverage <= nread-nerror)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0)then !avoid loop forth and back
       ndelta=real(ndelta1/2.d0,8) !decreasing the step
       xmu=xmu+real(nindex,8)*ndelta
    else
       ndelta=ndelta1
       xmu=xmu+real(nindex,8)*ndelta
    endif
    write(*,"(A,2f15.12,A,f15.12,A,I3,f15.12)")"mu,n=",xmu,naverage,"/",nread,"| ",nindex,ndelta
    if(abs(naverage-nread)>nerror)convergence=.false.
    call splot("muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
  end subroutine search_mu


  subroutine get_initial_sigma()
    logical                :: check1,check2,check
    inquire(file="Sigma_realw.last",exist=check1)
    inquire(file="Self_realw.last",exist=check2)
    check=check1*check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_realw.last",wr,sigma(1,:))
       call sread("Self_realw.last",wr,sigma(2,:))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 ;  delta=deltasc 
       sigma(2,:)=-delta ; sigma(1,:)=zero ; sold=sigma
    endif
  end subroutine get_initial_sigma

end program hmipt
