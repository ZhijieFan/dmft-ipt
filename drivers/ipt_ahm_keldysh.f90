program ahmk
  USE DMFT_IPT 
  USE DMFT_TOOLS
  USE SCIFOR
  implicit none
  integer                :: i,ik,Lk,iloop,Lm,L,Lf
  logical                :: converged
  complex(8)             :: zdet,zeta1,zeta2,x1,x2
  real(8)                :: delta,n,A,B,nf,dzeta1,dzeta2,D,dt,fmesh
  !
  complex(8),allocatable :: sigma(:,:),fg(:,:),wf0(:,:),calG(:,:),dummy(:,:)
  complex(8),allocatable :: zeta(:)
  !
  real(8),allocatable    :: wr(:),wrx(:)
  !
  call parse_input_variable(Lk,'Lk','inputIPT.in',default=1000)
  call parse_input_variable(L,"L",'inputIPT.in',default=10000)    
  call parse_input_variable(Lf,'Lf','inputIPT.in',default=5000)
  call parse_input_variable(D,'wband','inputIPT.in',default=1d0)

  call read_input("inputIPT.in")


  allocate(fg(2,L))
  allocate(sigma(2,L))
  allocate(wf0(2,L))
  allocate(calG(2,L))
  allocate(zeta(L))

  allocate(wr(L))
  wr =  linspace(-wmax,wmax,L,mesh=fmesh)
  print*,fmesh

  call get_initial_sigma

  iloop=0    ; converged=.false.


  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     !GETGLOC:
     fg=zero
     zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
     do i=1,L
        zeta1 = zeta(i)
        zeta2 = conjg(zeta(L+1-i))
        x1 = 0.5d0*((zeta1+zeta2) + sqrt((zeta1-zeta2)**2 - 4.d0*conjg(sigma(2,L+1-i))*sigma(2,i) ))
        x2 = 0.5d0*((zeta1+zeta2) - sqrt((zeta1-zeta2)**2 - 4.d0*conjg(sigma(2,L+1-i))*sigma(2,i) ))
        fg(1,i) = zeta2/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
        fg(2,i) =-conjg(sigma(2,L+1-i))/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
     enddo
     delta=-uloc(1)*trapz(fmesh,dimag(fg(2,:))*fermi(wr(:),beta))/pi
     n=-trapz(fmesh,dimag(fg(1,:))*fermi(wr,beta))/pi

     !GET THE WEISS FIELD \calG0^-1(w)
     ! wf(1)=calG0^-1 = G*(-w)/(G(w)G*(-w) + F(w)F(-w)) + Sigma(w)
     ! wf(2)=calF0^-1 =-F(-w) /(G(w)G*(-w) + F(w)F(-w)) + S(w)
     do i=1,L
        zdet     = fg(1,i)*conjg(fg(1,L+1-i)) + conjg(fg(2,L+1-i))*fg(2,i)
        wf0(1,i)= conjg(fg(1,L+1-i))/zdet  + sigma(1,i)   +  uloc(1)*(n-0.5d0)
        wf0(2,i)= conjg(fg(2,L+1-i))/zdet  + sigma(2,i)   +  delta
     end do
     do i=1,L
        zdet      =  wf0(1,i)*conjg(wf0(1,L+1-i)) + conjg(wf0(2,L+1-i))*wf0(2,i)
        calG(1,i)=  conjg(wf0(1,L+1-i))/zdet
        calG(2,i)=  conjg(wf0(2,L+1-i))/zdet
     end do

     sigma = ipt_solve_keldysh_sc(calG,delta,wmax)

     write(*,"(2f14.9)",advance="no")2.d0*n,delta
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=dmft_error,N1=Nsuccess,N2=Nloop)
     !if(printf)call splot("observables.ipt",delta,xmu,u,n,beta,dble(iloop),append=.true.)
  enddo

  call splot("observables.last",delta,xmu,uloc(1),n,beta,dble(iloop))
  call splot("Sigma_realw.restart",wr,sigma(1,:))
  call splot("Self_realw.restart",wr,sigma(2,:))

  !REDUCE size for printing.
  if(mod(Lf,2)/=0)Lf=Lf+1
  allocate(wrx(Lf),dummy(2,Lf))
  wrx = linspace(-wmax,wmax,Lf)

  call cubic_spline(wr(:),fg(1,:),wrx,dummy(1,:))
  call splot("DOS.last",wrx,-dimag(dummy(1,:))/pi)
  !
  call cubic_spline(wr,fg(1,:),wrx,dummy(1,:))
  call cubic_spline(wr,fg(2,:),wrx,dummy(2,:))
  call splot("G_realw.last",wrx,dummy(1,:))
  call splot("F_realw.last",wrx,dummy(2,:))
  !
  call cubic_spline(wr,calG(1,:),wrx,dummy(1,:))
  call cubic_spline(wr,calG(2,:),wrx,dummy(2,:))
  call splot("calG0_realw.last",wrx,dummy(1,:))
  call splot("calF0_realw.last",wrx,dummy(2,:))
  !
  call cubic_spline(wr,sigma(1,:),wrx,dummy(1,:))
  call cubic_spline(wr,sigma(2,:),wrx,dummy(2,:))
  call splot("Sigma_realw.last",wrx,dummy(1,:))
  call splot("Self_realw.last",wrx,dummy(2,:))
  !
  deallocate(wrx,dummy)







  !##################################################################
  !##################################################################
contains 
  !##################################################################
  !##################################################################


  subroutine get_initial_sigma()
    logical                :: check1,check2,check
    inquire(file="Sigma_realw.restart",exist=check1)
    inquire(file="Self_realw.restart",exist=check2)
    check=check1.AND.check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_realw.restart",wr,sigma(1,:))
       call sread("Self_realw.restart",wr,sigma(2,:))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 ;  delta=deltasc 
       sigma(2,:)=-delta ; sigma(1,:)=zero
    endif
  end subroutine get_initial_sigma





end program ahmk
