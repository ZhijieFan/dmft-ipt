!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt
  USE DMFT_IPT 
  USE DMFT_TOOLS
  USE SCIFOR
  implicit none
  integer    :: i,ik,Lk,iloop
  logical    :: converged
  complex(8) :: det,zeta1,zeta2,x1,x2
  real(8)    :: delta,n,D,wmix,fmesh
  !
  complex(8),allocatable          :: sigma(:,:),fg(:,:),wf0(:,:),calG(:,:),sold(:,:)
  complex(8),allocatable          :: zeta(:)
  !
  real(8),allocatable :: wt(:),epsik(:),wr(:),wm(:),tau(:)

  call parse_input_variable(D,"wband","inputIPT.in",default=1d0)
  call parse_input_variable(wmix,"wmix","inputIPT.in",default=0.75d0)
  call parse_input_variable(Lk,"Lk","inputIPT.in",default=1000)
  call read_input("inputIPT.in")

  allocate(wr(L))
  wr = linspace(-wmax,wmax,L,mesh=fmesh)
  print*,fmesh
  allocate(fg(2,L),sigma(2,L))
  allocate(wf0(2,L),calG(2,L))
  allocate(zeta(L),sold(2,L))

  allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D)

  call get_initial_sigma

  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     fg=zero
     zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
     do i=1,L
        zeta1 = zeta(i)
        zeta2 = conjg(zeta(L+1-i))
        x1 = 0.5d0*((zeta1+zeta2) + sqrt((zeta1-zeta2)**2 - 4.d0*conjg(sigma(2,L+1-i))*sigma(2,i) ))
        x2 = 0.5d0*((zeta1+zeta2) - sqrt((zeta1-zeta2)**2 - 4.d0*conjg(sigma(2,L+1-i))*sigma(2,i) ))
        fg(1,i) = zeta2/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
        fg(2,i) =-conjg(sigma(2,L+1-i))/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
        ! do ik=1,Lk
        !    det = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + conjg(sigma(2,L+1-i))*sigma(2,i)
        !    fg(1,i)=fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det
        !    fg(2,i)=fg(2,i) - wt(ik)*conjg(sigma(2,L+1-i))/det
        ! enddo
     enddo
     delta=-uloc*trapz(fmesh,dimag(fg(2,:))*fermi(wr,beta))/pi
     n    =-trapz(fmesh,dimag(fg(1,:))*fermi(wr,beta))/pi

     !get the Weiss Field \calG0^-1(w) = Gloc^-1 + Sigma + Sigma_HFB
     do i=1,L
        det     = fg(1,i)*conjg(fg(1,L+1-i)) + conjg(fg(2,L+1-i))*fg(2,i)
        wf0(1,i)= conjg(fg(1,L+1-i))/det  + sigma(1,i)   +  uloc*(n-0.5d0)
        wf0(2,i)= conjg(fg(2,L+1-i))/det  + sigma(2,i)   + delta
     end do
     do i=1,L
        det      =  wf0(1,i)*conjg(wf0(1,L+1-i)) + conjg(wf0(2,L+1-i))*wf0(2,i)
        calG(1,i)=  conjg(wf0(1,L+1-i))/det
        calG(2,i)=  conjg(wf0(2,L+1-i))/det
     end do

     write(*,"(3f14.9)",advance="no")2.d0*n,delta,delta/uloc
     sold=sigma
     sigma =  ipt_solve_real_sc(calG(1:2,:),wr,delta,L)
     sigma = wmix*sigma + (1.d0-wmix)*sold
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=dmft_error,N1=Nsuccess,N2=Nloop)
  enddo

  call splot("observables.ipt",xmu,uloc,dble(iloop),n,delta/uloc,delta)
  call splot("G_realw.ipt",wr,-dimag(fg(1,:))/pi,dreal(fg(1,:)))
  call splot("F_realw.ipt",wr,-dimag(fg(2,:))/pi,dreal(fg(2,:)))
  call splot("Sigma_realw.ipt",wr,sigma(1,:))
  call splot("Self_realw.ipt",wr,sigma(2,:))
  call splot("calG0_realw.ipt",wr,calG(1,:))
  call splot("calF0_realw.ipt",wr,calG(2,:))

contains 


  subroutine get_initial_sigma()
    logical                :: check1,check2,check
    inquire(file="Sigma_realw.ipt",exist=check1)
    inquire(file="Self_realw.ipt",exist=check2)
    check=check1.AND.check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_realw.ipt",wr,sigma(1,:))
       call sread("Self_realw.ipt",wr,sigma(2,:))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 ;  delta=deltasc 
       sigma(2,:)=-delta ; sigma(1,:)=zero ; sold=sigma
    endif
  end subroutine get_initial_sigma

end program hmipt
