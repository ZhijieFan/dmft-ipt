!########################################################
!     Program  : AHMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the attractive Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci & Antonio Privitera
!########################################################
program ahmmpt
  USE DMFT_IPT
  USE SCIFOR_VERSION
  USE IOTOOLS
  USE FUNCTIONS
  USE INTEGRATE
  !USE SQUARE_LATTICE
  implicit none
  integer                :: i,ik,Lk,iloop
  logical                :: converged
  real(8)                :: adummy,wband,de
  logical                :: check1,check2,check

  complex(8)             :: zeta1,zeta2,det
  real(8)                :: n,delta,n0,delta0
  !
  complex(8),allocatable :: sigma(:,:),fg(:,:)
  complex(8),allocatable :: wf0(:,:),calG(:,:)
  complex(8),allocatable :: sold(:,:),zeta(:)
  !
  real(8),allocatable    :: wt(:),epsik(:),wr(:)


  include "revision.inc"
  call version(revision)
  call read_input("inputIPT.in")



  print*,"we are actually using",L,"frequencies"
  allocate(wr(L))
  wr = linspace(-wmax,wmax,L,mesh=fmesh)
  print*,"|omegamax|",wr(L)

  allocate(fg(2,L),sigma(2,L))
  allocate(wf0(2,L),calG(2,L))
  allocate(sold(2,L),zeta(L))


  !build square lattice structure:
  ! Lk   = square_lattice_dimension(Nx)
  ! allocate(wt(Lk),epsik(Lk))
  ! wt   = square_lattice_structure(Lk,Nx)
  ! epsik= square_lattice_dispersion_array(Lk,ts)
  allocate(wt(L),epsik(L))
  wband=4.d0*ts+0.1d0
  epsik = linspace(-wband,wband,L,mesh=de)
  do i=1,L
     wt(i)=dens_2dsquare(epsik(i),ts)
  enddo
  wt=wt/trapz(de,wt)
  call splot("DOS2d.ipt",epsik,wt)
  wt = wt*de



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
        do ik=1,L
           det     = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + conjg(sigma(2,L+1-i))*sigma(2,i)
           fg(1,i) =fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det
           fg(2,i) =fg(2,i) - wt(ik)*conjg(sigma(2,L+1-i))/det
        enddo
     enddo

     n    = -trapz(fmesh,dimag(fg(1,:))*fermi(wr,beta))/pi!sum(dimag(fg(1,:))*fermi(wr,beta))*fmesh/pi ! densita' per spin
     delta= -u*trapz(fmesh,dimag(fg(2,:))*fermi(wr,beta))/pi!(u*sum(dimag(fg(2,:))*fermi(wr,beta))*fmesh/pi) 


     !Hartree corrected WF is: xmu=xmu0
     !\tilde{\calG0} = [G^-1 + Sigma - \Sigma_HFB]^-1
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

     n0    =  -trapz(fmesh,dimag(calG(1,:))*fermi(wr,beta))/pi!-sum(imag(calG(1,:))*fermi(wr,beta))*fmesh/pi
     delta0=  -u*trapz(fmesh,dimag(calG(2,:))*fermi(wr,beta))/pi!(u*sum(dimag(calG(2,:))*fermi(wr,beta))*fmesh/pi)

     write(*,"(4(f16.12))",advance="no"),n,n0,delta,delta0     

     sigma =  solve_mpt_sc_sopt(calG,wr,n,n0,delta,delta0,L)
     sigma = weight*sigma + (1.d0-weight)*sold ; sold=sigma
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=eps_error,N1=Nsuccess,N2=Nloop)

     if(nread/=0.d0)call search_mu(converged)

     call splot("nVSiloop.ipt",iloop,n,append=.true.)
     call splot("deltaVSiloop.ipt",iloop,delta,append=.true.)
     call splot("DOS.ipt",wr,-aimag(fg(1,:))/pi,append=printf)
     call splot("G_realw.ipt",wr,fg(1,:),append=printf)
     call splot("F_realw.ipt",wr,fg(2,:),append=printf)
     call splot("Sigma_realw.ipt",wr,sigma(1,:),append=printf)
     call splot("Self_realw.ipt",wr,sigma(2,:),append=printf)
     call splot("calG0_realw.ipt",wr,calG(1,:),append=printf)
     call splot("calF0_realw.ipt",wr,calG(2,:),append=printf)
     call splot("observables.ipt",iloop,beta,xmu,u,n,n0,delta,delta0,append=printf)
  enddo

  !ndelta=0.01d0  !  serve per evitare che il deltan piccolissimo in uscita
  !  venga resettato da file . Usa chi_trial e' meglio

  open(10,file="used.inputIPT.in")
  write(10,nml=variables) ! in questo modo il potenziale chimico e' updatato e non riparto sempre dallo stesso
  close(10)
  ! problema: lo step nel potenziale chimico cosi' e' salvato e nel ciclo gli step diventano piccolissimi che e' un casino  densita' e' 

  call close_file("nVSiloop.ipt")
  call close_file("deltaVSiloop.ipt")
  call splot("DOS.last",wr,-aimag(fg(1,:))/pi,append=.false.)
  call splot("G_realw.last",wr,fg(1,:),append=.false.)
  call splot("F_realw.last",wr,fg(2,:),append=.false.)
  call splot("Sigma_realw.last",wr,sigma(1,:),append=.false.)
  call splot("Self_realw.last",wr,sigma(2,:),append=.false.)
  call splot("calG0_realw.last",wr,calG(1,:),append=.false.)
  call splot("calF0_realw.last",wr,calG(2,:),append=.false.)
  call splot("observables.last",iloop,beta,xmu,u,n,n0,delta,delta0,append=.false.)


contains 

  subroutine search_mu(convergence)
    integer, save         ::nindex
    integer               ::nindex1
    real(8)               :: naverage,ndelta1
    logical,intent(inout) :: convergence
    naverage=2.d0*n
    nindex1=nindex
    ndelta1=ndelta
    if((naverage >= nread+nerror))then
       nindex=-1
    elseif(naverage <= nread-nerror)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
       ndelta=ndelta1/2.d0
    else
       ndelta=ndelta1
    endif
    xmu=xmu+real(nindex,8)*ndelta
    write(*,"(A,2f15.12,A,f15.12,A,I3,f15.12)")"mu,n=",xmu,naverage,"/",nread,"| ",nindex,ndelta
    if(abs(naverage-nread)>nerror)convergence=.false.
    call splot("muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
  end subroutine search_mu


  subroutine get_initial_sigma()
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

end program ahmmpt
