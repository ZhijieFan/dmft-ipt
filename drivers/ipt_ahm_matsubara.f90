program hmipt_matsubara
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                :: i,ik,Lk,iloop,L
  logical                :: converged,bcsflag
  complex(8)             :: zeta,zeta1,zeta2,x1,x2,zsqrt
  real(8)                :: n,delta,D,cdet,delta0,Ds(2),alpha
  !
  complex(8),allocatable :: fg(:,:),fg0(:,:),sigma(:,:),calG(:,:),Gdet(:),calG_prev(:,:)
  !
  real(8),allocatable    :: wt(:),epsik(:),wm(:)

  call parse_input_variable(L,"L","inputIPT.conf",default=4096)
  call parse_input_variable(D,"wband","inputIPT.conf",default=1d0)
  call parse_input_variable(Lk,"Lk","inputIPT.conf",default=1000)
  call parse_input_variable(alpha,"alpha","inputIPT.conf",default=0.5d0)
  call parse_input_variable(bcsflag,"BCSFLAG","inputIPT.conf",default=.false.)
  call read_input("inputIPT.conf")

  allocate(wm(L))
  wm  = pi/beta*(2*arange(1,L)-1)

  allocate(fg(2,L))
  allocate(fg0(2,L),sigma(2,L),calG(2,L),Gdet(L),calG_prev(2,L))

  allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D)

  call get_initial_sigma

  iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop

     fg=zero
     do i=1,L
        zeta =  xi*wm(i) + xmu - sigma(1,i)
        fg(:,i)=zero
        do ik=1,Lk
           cdet = abs(zeta-epsik(ik))**2 + (sigma(2,i))**2
           fg(1,i)=fg(1,i) + wt(ik)*(conjg(zeta)-epsik(ik))/cdet
           fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/cdet
        enddo
     enddo


     n = ipt_measure_dens_matsubara(fg(1,:))
     delta = Uloc(1)*ipt_measure_phi_matsubara(fg(2,:))

     ! Ds(1) = get_ds1()
     ! Ds(2) = get_ds2()

     Gdet     = abs(fg(1,:))**2 + (fg(2,:))**2
     fg0(1,:) =  conjg(fg(1,:))/Gdet + sigma(1,:)
     fg0(2,:) =  fg(2,:)/Gdet        + sigma(2,:) + delta


     Gdet      =  abs(fg0(1,:))**2 + (fg0(2,:))**2
     calG(1,:) =  conjg(fg0(1,:))/Gdet
     calG(2,:) =  fg0(2,:)/Gdet

     write(*,"(5f14.9)",advance="no")2*n,delta!,Ds(1),Ds(2)
     if(iloop>1)calG= alpha*calG + (1d0-alpha)*calG_prev
     calG_prev =  calG

     if(bcsflag)then
        sigma(2,:)=-delta ; sigma(1,:)=zero
     else
        sigma =  ipt_solve_matsubara(calG,delta)
     endif
     converged = check_convergence(calG(1,:)+calG(2,:),eps=dmft_error,N1=Nsuccess,N2=Nloop)

  enddo

  call splot("Sigma_iw.ipt",wm,sigma(1,:),append=.false.)
  call splot("Self_iw.ipt",wm,sigma(2,:),append=.false.)
  call splot("G_iw.ipt",wm,fg(1,:),append=.false.)
  call splot("F_iw.ipt",wm,fg(2,:),append=.false.)
  call splot("calG_iw.ipt",wm,calG(1,:),append=.false.)
  call splot("calF_iw.ipt",wm,calG(2,:),append=.false.)
  ! call splot("observables.ipt",uloc(1),n,delta,append=.false.)
  !Ds(1),Ds(2)

  ! call get_sc_internal_energy

contains


  subroutine get_initial_sigma()
    logical :: check1,check2,check
    inquire(file="Sigma_iw.restart",exist=check1)
    inquire(file="Self_iw.restart",exist=check2)
    check=check1.AND.check2
    if(check)then
       write(*,*)"Reading Sigma in input:"
       call sread("Sigma_iw.restart",wm,sigma(1,:))
       call sread("Self_iw.restart",wm,sigma(2,:))
    else
       print*,"Using Hartree-Fock self-energy"
       print*,"===================================="
       n=0.5d0 ; delta=deltasc
       sigma(2,:)=-delta ; sigma(1,:)=zero
    endif
  end subroutine get_initial_sigma


  ! function get_ds1 result(Ds)
  !   integer    :: i,ik
  !   real(8)    :: Ds,nk,vk,Dssum,matssum,dens,chisum,Chij
  !   real(8)    :: Ekin
  !   complex(8) :: iw,gkw,fkw,g0kw,f0kw
  !   real(8)    :: Sigma_infty,S_infty,det,det_infty,csi,Ei,free
  !   !Get asymptotic self-energies
  !   Sigma_infty =   dreal(sigma(1,L))
  !   S_infty     =   dreal(sigma(2,L))
  !   !
  !   !1. Ds = 4/beta\sum_mats \int_ek N(ek)V(ek)F(ek,iw)F(ek,iw)
  !   !1.:
  !   !N(ek) = wtk(ik) (DOS)
  !   !V(ek) = (D**2 - ek**2)/3.d0
  !   Ds=0.d0                       ! superfluid stiffness (Bethe)
  !   do ik=1,Lk
  !      Dssum   = 0.d0
  !      nk = wt(ik)
  !      vk = (D**2-epsik(ik)**2)/3.d0
  !      do i=1,L
  !         iw       = xi*wm(i)
  !         det      = abs(iw+xmu-epsik(ik)-sigma(1,i))**2 + dreal(sigma(2,i))**2
  !         gkw = (-iw+xmu - epsik(ik) - conjg(sigma(1,i)) )/det
  !         fkw = -sigma(2,i)/det
  !         Dssum   = Dssum    +  fkw*fkw
  !      enddo
  !      Ds=Ds + 2d0*4d0/beta*nk*vk*Dssum
  !   enddo
  ! end function get_ds1





  ! function get_ds2 result(Ds)
  !   integer    :: i,ik
  !   real(8)    :: Ds,nk,vk,Dssum,matssum,dens,chisum,Chij
  !   real(8)    :: Ekin
  !   complex(8) :: zeta,g11kw,g22kw,fkw,g0kw,f0kw
  !   real(8)    :: Sigma_infty,S_infty,det,det_infty,csi,Ei,free
  !   !Get asymptotic self-energies
  !   Sigma_infty =   dreal(sigma(1,L))
  !   S_infty     =   dreal(sigma(2,L))
  !   !
  !   !2. Ds = -Ekin - \chi_jj(q=0,\omega=0)
  !   !
  !   !2. 
  !   ! - get Ekin:
  !   Ekin=0d0
  !   Chij=0d0
  !   do ik=1,Lk
  !      vk    = (D**2-epsik(ik)**2)/3.d0
  !      csi   = epsik(ik)-(xmu-Sigma_infty)
  !      Ei    = dsqrt(csi**2 + S_infty**2)
  !      free  = 0.5d0*(1d0 - csi/Ei)*dtanh(0.5d0*beta*Ei)
  !      matssum=0d0
  !      chisum =0d0
  !      do i=1,L
  !         zeta      = xi*wm(i)+xmu-sigma(1,i)
  !         det       = abs(zeta-epsik(ik))**2 + dreal(sigma(2,i))**2
  !         g11kw     = conjg(zeta - epsik(ik))/det
  !         fkw       = -sigma(2,i)/det
  !         !
  !         det_infty = wm(i)**2 + (epsik(ik)-(xmu-Sigma_infty))**2 + S_infty**2
  !         g0kw      = (-xi*wm(i) - (epsik(ik)-(xmu-Sigma_infty)))/det_infty
  !         f0kw      =-S_infty/det_infty
  !         !
  !         matssum   =  matssum +  dreal(g11kw)-dreal(g0kw)
  !         chisum    =  chisum  + abs(g11kw)**2 + dreal(fkw)**2
  !      enddo
  !      nk    = 4d0/beta*matssum + 2d0*free
  !      chisum=-2d0/beta*chisum
  !      Ekin  = Ekin   + wt(ik)*nk*epsik(ik)
  !      Chij  = Chij   + wt(ik)*vk*chisum
  !   enddo
  !   Ds = -Ekin-Chij
  ! end function get_ds2




  ! subroutine get_sc_internal_energy
  !   integer    :: i,ik
  !   real(8)    :: matssum,fmatssum,checkP,checkdens,vertex,Dssum
  !   complex(8) :: iw,gkw,fkw,g0kw,f0kw
  !   real(8)    :: Epot,Etot,Eint,kin,kinsim,Ds,docc
  !   real(8)    :: Sigma_infty,S_infty,det,det_infty,csi,Ei,thermal_factor
  !   real(8)    :: free(Lk),Ffree(Lk),n_k(Lk)

  !   !Get asymptotic self-energies
  !   Sigma_infty =   real(sigma(1,L),8)
  !   S_infty     =   real(sigma(2,L),8)

  !   checkP=0.d0 ; checkdens=0.d0 ;          ! test variables

  !   kin=0.d0                      ! kinetic energy (generic)
  !   Ds=0.d0                       ! superfluid stiffness (Bethe)

  !   do ik=1,Lk

  !      csi            = epsik(ik)-(xmu-Sigma_infty)
  !      Ei             = dsqrt(csi**2 + S_infty**2)
  !      thermal_factor = dtanh(0.5d0*beta*Ei)
  !      free(ik)        = 0.5d0*(1.d0 - csi/Ei)*thermal_factor
  !      Ffree(ik)       =-(0.5d0*S_infty)/Ei*thermal_factor

  !      fmatssum= 0.d0
  !      matssum = 0.d0
  !      Dssum   = 0.d0

  !      vertex=(D**2-epsik(ik)**2)/3.d0

  !      do i=1,L
  !         iw       = xi*wm(i)
  !         det      = abs(iw+xmu-epsik(ik)-sigma(1,i))**2 + real(sigma(2,i),8)**2
  !         det_infty= wm(i)**2 + (epsik(ik)-(xmu-Sigma_infty))**2 + S_infty**2

  !         gkw = (-iw+xmu - epsik(ik) - conjg(sigma(1,i)) )/det
  !         fkw = -sigma(2,i)/det

  !         g0kw= (-iw - (epsik(ik)-(xmu-Sigma_infty)))/det_infty
  !         f0kw=-S_infty/det_infty

  !         matssum =  matssum +  real(gkw,8)-real(g0kw,8)
  !         !        matssum =matssum + real(gkw,8)        ! without tails corrections

  !         fmatssum= fmatssum +  real(fkw,8)-real(f0kw,8)
  !         Dssum   = Dssum    +  fkw*fkw

  !      enddo

  !      n_k(ik)   = 4.d0/beta*matssum + 2.d0*free(ik)
  !      !    n_k(ik)   = 4.d0/beta*matssum + 1.d0     ! without tails corrections
  !      checkP    = checkP    - wt(ik)*(2.d0/Beta*fmatssum+Ffree(ik))
  !      !    print*,checkP,Ffree(ik)
  !      checkdens = checkdens + wt(ik)*n_k(ik)
  !      kin    = kin    + wt(ik)*n_k(ik)*epsik(ik)
  !      Ds=Ds + 8.d0/beta* wt(ik)*vertex*Dssum
  !   enddo

  !   !  call splot("nk0_distribution.last",epsik,2.d0*free)
  !   !  call splot("fnk0_distribution.last",epsik,Ffree)

  !   kinsim=0
  !   kinsim = sum(fg(1,:)*fg(1,:)+conjg(fg(1,:))*conjg(fg(1,:))-2.d0*fg(2,:)*fg(2,:))*D**2/beta
  !   ! kinsim = kinsim - 1.d0/(2*pi*wm(L)) !check out where this term comes from??!!

  !   Epot=zero
  !   Epot = sum(fg(1,:)*sigma(1,:) + fg(2,:)*sigma(2,:))/beta*2.d0

  !   docc = n**2
  !   if(uloc(1) > 0.01d0)docc=-Epot/uloc(1) + n - 0.25d0


  !   Eint=kin+Epot

  !   ! Ds=zero
  !   ! Ds = sum(fg(2,:)*fg(2,:))/beta*2.d0

  !   ! kinsim=zero
  !   ! do i=1,L
  !   !    kinsim=kinsim + fg(1,i)*fg(1,i) + conjg(fg(1,i))*conjg(fg(1,i)) - 2.d0*fg(2,i)*fg(2,i)
  !   ! enddo
  !   ! kinsim=0.5d0*kinsim/Beta
  !   ! kinsim=kinsim- 1.d0/(2*pi*om(Iwmax2))

  !   write(*,*)'========================================='      
  !   write(*,*)"Asymptotic Self-Energies",Sigma_infty, S_infty
  !   write(*,*)'========================================='
  !   write(*,*)"n,delta",n,delta
  !   write(*,*)"Dn% ,Ddelta%",(n-0.5d0*checkdens)/n,(delta + uloc*checkP)/delta ! u is positive
  !   write(*,*)'========================================='
  !   write(*,*)"Kinetic energy",kin
  !   write(*,*)'========================================='
  !   write(*,*)"double occupancy   =",docc
  !   write(*,*)'========================================='
  !   write(*,*) 'Kinetic Energy TEST (simple formula)'
  !   write(*,*) '###ACTHUNG: FOR BETHE ONLY####',kinsim
  !   write(*,*) 'Dkin%',(kin-kinsim)/kin
  !   write(*,*)'========================================='
  !   write(*,*) 'Superfluid stiffness',Ds
  !   write(*,*) 'Potential Energy U(n_up-1/2)(n_do-1/2)',Epot
  !   write(*,*) 'Internal Energy',Eint
  !   write(*,*)'========================================='
  !   call splot("nk_distribution.ipt",epsik,n_k,2.d0*free)
  !   call splot("thermodynamics.ipt",L,n,0.5d0*checkdens,kin,kinsim,docc,Ds,append=.true.)

  !   !  call splot("fnk0_distribution.ipt",epsik,Ffree)
  !   !   write(42,143)xmu,totdens,checkdens,dreal(kinsim),Tcheck,double
  !   ! 143 format(f6.4,5(1x,f12.9))
  !   return 
  ! end subroutine get_sc_internal_energy

end program hmipt_matsubara
