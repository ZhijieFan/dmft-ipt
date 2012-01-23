!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt_matsubara
  USE DMFT_IPT
  implicit none
  integer                :: ik,Lk
  logical                :: converged
  complex(8)             :: zeta,cdet
  real(8)                :: n,delta
  !
  complex(8),allocatable :: fg(:,:),fg0(:,:),sigma(:,:),calG(:,:),det(:)
  real(8),allocatable    :: fgt(:,:)
  !
  real(8),allocatable    :: wt(:),epsik(:)

  call read_input("inputIPT.in")
  allocate(fg(2,L),fgt(2,0:L))
  allocate(fg0(2,L),sigma(2,L),calG(2,L),det(L))


  Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D_=D,eps_=eps)


  delta=deltasc
  sigma=zero ; sigma(2,:)=-delta 
  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     fg=zero
     do i=1,L
        zeta =  xi*wm(i) + xmu - sigma(1,i)
        do ik=1,Lk
           cdet = abs(zeta-epsik(ik))**2 + (sigma(2,i))**2
           fg(1,i)=fg(1,i) + wt(ik)*(conjg(zeta)-epsik(ik))/cdet
           fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/cdet
        enddo
     enddo
     call fftgf_iw2tau(fg(1,:),fgt(1,0:L),beta)
     call fft_iw2tau(fg(2,:),fgt(2,0:L),beta,L) !; fgt(2,:)=-fgt(2,:)
     n=-real(fgt(1,L),8) ; delta= -u*fgt(2,0)


     det       = abs(fg(1,:))**2 + (fg(2,:))**2
     fg0(1,:) =  conjg(fg(1,:))/det + sigma(1,:)
     fg0(2,:) =  fg(2,:)/det        + sigma(2,:) + delta

     det       =  abs(fg0(1,:))**2 + (fg0(2,:))**2
     calG(1,:) =  conjg(fg0(1,:))/det
     calG(2,:) =  fg0(2,:)/det


     write(*,"(4(f16.12))",advance="no")delta
     sigma =  solve_ipt_sc_matsubara(calG,fg,delta)
     converged = check_convergence(sigma(1,:)+sigma(2,:))
  enddo

  call splot("Sigma_iw.last",wm,sigma(1,:),append=TT)
  call splot("Self_iw.last",wm,sigma(2,:),append=TT)
  call splot("G_iw.last",wm,fg(1,:),append=TT)
  call splot("F_iw.last",wm,fg(2,:),append=TT)
  call splot("calG_iw.last",wm,calG(1,:),append=TT)
  call splot("calF_iw.last",wm,calG(2,:),append=TT)
  call splot("n.delta_iw.last",n,delta,append=TT)
end program hmipt_matsubara
