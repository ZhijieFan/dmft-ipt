!########################################################
!     Program  : AHMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the attractive Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci & Antonio Privitera
!########################################################
program ahmmpt
  ! USE IPT_VARS_GLOBAL
  ! USE MPT_SC_FUNX_SOPT
  USE DMFT_IPT
  USE IOTOOLS

  implicit none
  integer                :: i,ik,Lk
  logical                :: converged
  real(8)                :: adummy
  logical                :: check1,check2,check

  complex(8)             :: zeta1,zeta2,det
  real(8)                :: n,delta,n0,delta0
  real(8)                :: summag,halfsummag,summag0,halfsummag0    ! test variables for debugging
  !
  complex(8),allocatable :: sigma(:,:),fg(:,:)
  complex(8),allocatable :: wf0(:,:),calG(:,:)
  complex(8),allocatable :: sold(:,:),zeta(:)
  !
  real(8),allocatable    :: wt(:),epsik(:),wr(:)

  call read_input("inputIPT.in")

  print*,"we are actually using",L,"frequencies"
  allocate(wr(-L:L))
  wr = linspace(-wmax,wmax,2*L+1,mesh=fmesh)
  print*,"|omegamax|",wr(L)
  ! if (nread==0)then 
  !    print*,"FIXED CHEMICAL POTENTIAL, xmu = ",xmu
  ! else
  !    print*,"FIXED density, n = ",nread
  ! endif

  allocate(fg(2,-L:L),sigma(2,-L:L))
  allocate(wf0(2,-L:L),calG(2,-L:L))
  allocate(sold(2,-L:L),zeta(-L:L))

  Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D=1.d0)!D_=D,eps_=eps)


  ! it is very important to read the self-energy from previous run if available
  ! this is even more important on the real axis where the self-energy is extremely structured


  call get_initial_sigma


  !  sigma(2,:)=-delta ; sigma(1,:)=-u*n; sold=sigma

  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     fg=zero
     zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
     do i=-L,L
        zeta1 = zeta(i)
        zeta2 = conjg(zeta(-i))
        do ik=1,Lk
           !           det     = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + sigma(2,i)*sigma(2,i)             !!  versione correntemente in uso
           det     = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + conjg(sigma(2,-i))*sigma(2,i)              ! versione di BAUER
           !           det     = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + sigma(2,i)*sigma(2,-i)  !! versione vecchia probably errata

           fg(1,i) =fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det

           !          fg(2,i) =fg(2,i) - wt(ik)*sigma(2,i)/det                                          !! versione correntemente in uso
           fg(2,i) =fg(2,i) - wt(ik)*conjg(sigma(2,-i))/det                                 ! versione di BAUER
           !          fg(2,i) =fg(2,i) - wt(ik)*sigma(2,-i)/det  ! versione vecchia probably errata

        enddo
     enddo

     !    la parte immaginaria della Gloc deve essere correttamente normalizzata..
     !    che succede se non lo e' ?


     !      do i=-L,L
     !         fg(1,i)=fg(1,i)/summa
     !      enddo

     n    = -sum(aimag(fg(1,:))*fermi(wr,beta))*fmesh/pi ! densita' per spin
     delta= -(u*sum(aimag(fg(2,:))*fermi(wr,beta))*fmesh/pi) 


     !     write(*,"(2(f16.12))",advance="no")n,delta
     !     stop 

     !     alla prima iterazione la simmetria sbagliata non cambia nulla sulle G e F 
     !     perche la self-energia BCS ha le simmetrie giuste.. Tuttavia anche in BCS 
     !     G e F non sono banali ecco perche le G0 sono sbagliate 

     !     fg=zero
     !     zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
     !     do i=-L,L
     !        zeta1 = zeta(i)
     !        zeta2 = conjg(zeta(-i))
     !        do ik=1,Lk
     !           det     = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + sigma(2,i)*sigma(2,-i)  
     !           fg(1,i) =fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det
     !           fg(2,i) =fg(2,i) - wt(ik)*sigma(2,-i)/det
     !        enddo
     !     enddo
     !     n    = -sum(aimag(fg(1,:))*fermi(wr,beta))*fmesh/pi
     !     delta= -(u*sum(aimag(fg(2,:))*fermi(wr,beta))*fmesh/pi)
     !

     !Hartree corrected WF is: xmu=xmu0
     !\tilde{\calG0} = [G^-1 + Sigma - \Sigma_HFB]^-1

     do i=-L,L
        det     = fg(1,i)*conjg(fg(1,-i)) + conjg(fg(2,-i))*fg(2,i)          !! versione di BAUER
        !       det     = fg(1,i)*conjg(fg(1,-i)) + fg(2,i)*fg(2,i)         !! versione correntemente in uso
        !        det     = fg(1,i)*conjg(fg(1,-i)) + fg(2,i)*fg(2,-i)       !! versione vecchia probably errata

        wf0(1,i)= conjg(fg(1,-i))/det  + sigma(1,i)   +  u*(n-0.5d0) ! versione vecchia particle-hole doppio segno sbagliato [ora corretto]
        !        wf0(1,i)= conjg(fg(1,-i))/det  + sigma(1,i)   +  u*n 

        !        ! stampo i valori di sigma HFB che sottraggo 
        !        write(49,*)i,real(sigma(1,i)),aimag(sigma(1,i)),  u*n 

        wf0(2,i)= fg(2,i)/det          + conjg(sigma(2,-i))   + delta     !! versione di BAUER
        !        wf0(2,i)= fg(2,i)/det          + sigma(2,i)   + delta       !! versione correntemente in uso
        !        wf0(2,i)= fg(2,i)/det          + sigma(2,-i)   + delta ! versione vecchia probably errata

        !        ! stampo i valori di sigma HFB che sottraggo 
        !         write(50,*)i,real(sigma(2,i)),aimag(sigma(2,i)),delta

     end do

     do i=-L,L
        det      =  wf0(1,i)*conjg(wf0(1,-i)) + conjg(wf0(2,-i))*wf0(2,i)   !! versione di BAUER
        !        det      =  wf0(1,i)*conjg(wf0(1,-i)) + wf0(2,i)*wf0(2,i)   !! versione correntemente in uso
        !        det      =  wf0(1,i)*conjg(wf0(1,-i)) + wf0(2,i)*wf0(2,-i) !! versione vecchia probably errata
        calG(1,i)=  conjg(wf0(1,-i))/det

        calG(2,i)=  conjg(wf0(2,-i))/det     !! versione di BAUER
        !        calG(2,i)=  wf0(2,i)/det     !! versione correntemente in uso
        !        calG(2,i)=  wf0(2,-i)/det   !! versione vecchia probably errata
     end do

     n0    =  -sum(aimag(calG(1,:))*fermi(wr,beta))*fmesh/pi
     delta0=  -(u*sum(aimag(calG(2,:))*fermi(wr,beta))*fmesh/pi)

     !    Normalization check

     summag=0.d0; halfsummag=0.d0; summag0=0.d0; halfsummag0=0.d0

     halfsummag=sum(-1.d0/pi*dimag(fg(1,-L:0))*fmesh)
     summag=sum(-1.d0/pi*dimag(fg(1,:))*fmesh)

     halfsummag0=sum(-1.d0/pi*dimag(calG(1,-L:0))*fmesh)
     summag0=sum(-1.d0/pi*dimag(calG(1,:))*fmesh)

     write(*,"(a10,2f10.6,a10,2f10.6)")"fg",summag,halfsummag,"calg0",summag0,halfsummag0


     !      do i=-L,L
     !          calG(1,i)=calG(1,i)/summa
     !      enddo

     write(*,"(4(f16.12))",advance="no"),n,n0,delta,delta0     

     !
     !     do i=-L,L
     !        det     = fg(1,i)*conjg(fg(1,-i)) + fg(2,i)*fg(2,-i)
     !        wf0(1,i)= conjg(fg(1,-i))/det  + sigma(1,i)   - u*(n-0.5d0)
     !        wf0(2,i)= fg(2,-i)/det          + sigma(2,i)   + delta
     !     end do
     !     do i=-L,L
     !        det      =  wf0(1,i)*conjg(wf0(1,-i)) + wf0(2,i)*wf0(2,-i)
     !        calG(1,i)=  conjg(wf0(1,-i))/det
     !        calG(2,i)=  wf0(2,-i)/det
     !     end do
     !     n0    =  -sum(aimag(calG(1,:))*fermi(wr,beta))*fmesh/pi
     !     delta0=  -(u*sum(aimag(calG(2,:))*fermi(wr,beta))*fmesh/pi)

     !sigma =  solve_mpt_sc_sopt(calG,fg,n,n0,delta,delta0)
     sigma =  solve_mpt_sc_sopt(calG,wr,n,n0,delta,delta0,L)
     sigma = weight*sigma + (1.d0-weight)*sold ; sold=sigma
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=eps_error,N1=Nsuccess,N2=Nloop)

     ! if(nread/=0.d0)call search_mu(converged)

     if (converged) then
        open(22,file="summary")
        write(22,"(6(f16.12),i5)") xmu,n,n0,delta,delta0,U,iloop
        close(22)
     endif

     call splot("DOS.ipt",wr,-aimag(fg(1,:))/pi,append=printf)
     call splot("G_realw.ipt",wr,fg(1,:),append=printf)
     call splot("F_realw.ipt",wr,fg(2,:),append=printf)
     call splot("Sigma_realw.ipt",wr,sigma(1,:),append=printf)
     call splot("Self_realw.ipt",wr,sigma(2,:),append=printf)
     call splot("calG0_realw.ipt",wr,calG(1,:),append=printf)
     call splot("calF0_realw.ipt",wr,calG(2,:),append=printf)
     call splot("n.delta_realw.ipt",n,delta,append=TT)

  enddo

  !ndelta=0.01d0  !  serve per evitare che il deltan piccolissimo in uscita
  !  venga resettato da file . Usa chi_trial e' meglio

  open(10,file="used.inputIPT.in")
  write(10,nml=variables) ! in questo modo il potenziale chimico e' updatato e non riparto sempre dallo stesso
  close(10)

  ! problema: lo step nel potenziale chimico cosi' e' salvato e nel ciclo gli step diventano piccolissimi che e' un casino  densita' e' 

  call splot("DOS.last",wr,-aimag(fg(1,:))/pi,append=FF)
  call splot("G_realw.last",wr,fg(1,:),append=FF)
  call splot("F_realw.last",wr,fg(2,:),append=FF)
  call splot("Sigma_realw.last",wr,sigma(1,:),append=FF)
  call splot("Self_realw.last",wr,sigma(2,:),append=FF)
  call splot("calG0_realw.last",wr,calG(1,:),append=FF)
  call splot("calF0_realw.last",wr,calG(2,:),append=FF)
  call splot("n.delta_realw.last",n,delta,append=TT)

contains 

  ! subroutine search_mu(convergence)
  !   integer, save         ::nindex
  !   integer               ::nindex1
  !   real(8)               :: naverage,ndelta1
  !   logical,intent(inout) :: convergence
  !   naverage=n
  !   nindex1=nindex
  !   ndelta1=ndelta
  !   if((naverage >= nread+nerror))then
  !      nindex=-1
  !   elseif(naverage <= nread-nerror)then
  !      nindex=1
  !   else
  !      nindex=0
  !   endif
  !   if(nindex1+nindex==0)then !avoid loop forth and back
  !      ndelta=real(ndelta1/2.d0,8) !decreasing the step
  !      xmu=xmu+real(nindex,8)*ndelta
  !   else
  !      ndelta=ndelta1
  !      xmu=xmu+real(nindex,8)*ndelta
  !   endif
  !   write(*,"(A,2f15.12,A,f15.12,A,I3,f15.12)")"mu,n=",xmu,naverage,"/",nread,"| ",nindex,ndelta
  !   if(abs(naverage-nread)>nerror)convergence=.false.
  !   call splot("muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
  ! end subroutine search_mu


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
