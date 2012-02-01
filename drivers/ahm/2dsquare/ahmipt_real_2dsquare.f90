!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt
  USE DMFT_IPT
  USE SQUARE_LATTICE
  implicit none
  integer    :: i,ik,Lk
  logical    :: converged
  complex(8) :: zeta1,zeta2,det
  complex(8),allocatable :: sigma(:,:),fg(:,:),fg0(:,:),calG(:,:)
  real(8),allocatable    :: wr(:),wt(:),epsik(:)

  call read_input("inputIPT.in")
  allocate(fg(2,-L:L),sigma(2,-L:L),fg0(2,-L:L),calG(2,-L:L))

  allocate(wr(-L:L))
  wr = linspace(-wmax,wmax,2*L+1,mesh=fmesh)

  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  sigma=zero ; sigma(2,:)=-deltasc ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     fg=zero
     do i=-L,L
        zeta1 = cmplx(wr(i),eps) - sigma(1,i)
        zeta2 = cmplx(wr(i),eps) + conjg(sigma(1,-i))
        do ik=1,Lk
           det = (zeta1-epsik(ik))*(zeta2+epsik(ik)) - sigma(2,i)**2
           fg(1,i)=fg(1,i) + wt(ik)*(zeta2+epsik(ik))/det
           fg(2,i)=fg(2,i) - wt(ik)*(sigma(2,i))/det
        enddo
        if(dimag(fg(1,i)) > 0.d0)fg(1,i)=real(fg(1,i),8)+ (0.d0,1.d-9)
        if(i<0.d0 .AND. dimag(fg(2,i))<0.d0)fg(2,i)=real(fg(2,i),8)+ (0.d0,1.d-9)
        if(i>0.d0 .AND. dimag(fg(2,i))>0.d0)fg(2,i)=real(fg(2,i),8)+ (0.d0,1.d-9)
     enddo

     deltasc=u*sum(aimag(fg(2,-L:L))*fermi(wr(-L:L),beta))*fmesh/pi

     do i=-L,L
        det      = fg(1,i)*conjg(fg(1,-i)) + fg(2,i)**2
        calG(1,i)= conjg(fg(1,-i))/det + sigma(1,i)
        calG(2,i)= fg(2,i)/det + sigma(2,i) - deltasc
     end do
     do i=-L,L
        det     =  calG(1,i)*conjg(calG(1,-i)) + calG(2,i)**2
        fg0(1,i)=  conjg(calG(1,-i))/det
        fg0(2,i)=  calG(2,i)/det
        if(dimag(fg0(1,i)) > 0.d0)fg0(1,i)=real(fg0(1,i),8)+ (0.d0,1.d-9)
        if(i<0.d0 .AND. dimag(fg0(2,i))<0.d0)fg0(2,i)=real(fg0(2,i),8)+ (0.d0,1.d-9)
        if(i>0.d0 .AND. dimag(fg0(2,i))>0.d0)fg0(2,i)=real(fg0(2,i),8)+ (0.d0,1.d-9)
     end do
     sigma =  solve_ipt_sc_sopt(fg0,wr,deltasc)
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps_error,nsuccess,nloop)
     call splot("Gloc_realw.ipt",wr,fg(1,:),append=printf)
     call splot("Floc_realw.ipt",wr,fg(2,:),append=printf)
  enddo

end program hmipt
