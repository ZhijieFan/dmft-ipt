!########################################################
!     Program  : AHMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Attractive Hubbard Model using IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt
  USE DMFT_IPT 
  USE SCIFOR_VERSION
  USE IOTOOLS
  USE INTERPOLATE
  USE INTEGRATE
  USE OPTIMIZE
  USE ERROR
  USE TIMER
  implicit none
  integer                      :: i,ik,Lk,iloop,Lm,p,q,Lf,Lmts
  logical                      :: converged
  complex(8)                   :: det,zeta1,zeta2,x1,x2
  real(8)                      :: delta,n,A,B,nf,nfL,nfR,ekin,dzeta1,dzeta2
  !
  complex(8),allocatable       :: sigma(:,:),fg(:,:),wf0(:,:),calG(:,:),dummy(:,:)
  complex(8),allocatable       :: zeta(:),fgm(:,:),sm(:,:)
  !
  real(8),allocatable          :: wt(:),epsik(:),wr(:),t(:),wrx(:),dos(:,:),en(:),nk(:),ddet(:),ipoles(:)
  !
  type(keldysh_equilibrium_gf) :: fg0k(2),sk(2),calG11,calG22,calF12,calF21
  !
  real(8)                      :: vbias,wxmax
  logical                      :: type,thermo,poles,octype
  type(finter_type)            :: det_finter


  include "revision.inc"
  call version(revision)

  call read_input("inputIPT.in")
  call parse_cmd_variable(vbias,'VBIAS',default=0.d0)
  call parse_cmd_variable(poles,"POLES",default=.false.)
  call parse_cmd_variable(octype,"OCTYPE",default=.false.)
  call parse_cmd_variable(thermo,"THERMO",default=.false.)
  call parse_cmd_variable(wxmax,"WXMAX",default=12.d0)
  call parse_cmd_variable(Lf,"LF",default=4096)
  call parse_cmd_variable(Lmts,"Lmts",default=4096)

  allocate(fg(2,L))
  allocate(sigma(2,L))
  allocate(wf0(2,L))
  allocate(calG(2,L))
  allocate(zeta(L))
  Lm=L/2
  call allocate_gf(fg0k(1),Lm);  call allocate_gf(fg0k(2),Lm)
  call allocate_gf(sk(1),Lm)  ;  call allocate_gf(sk(2),Lm)
  call allocate_gf(calG11,Lm)
  call allocate_gf(calG22,Lm)
  call allocate_gf(calF12,Lm)
  call allocate_gf(calF21,Lm)

  allocate(wr(L),t(-L/2:L/2))
  wr =  linspace(-wmax,wmax,L,mesh=fmesh)
  dt  = pi/wmax
  t   = linspace(-dt*Lm,dt*Lm,L+1,mesh=dt)

  call msg("MESH="//txtfy(fmesh))

  D=2.d0*ts

  call get_initial_sigma

  iloop=0    ; converged=.false.

  if(poles)then
     converged=.true.
     Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk),nk(Lk),ddet(L),ipoles(Lk))
     call bethe_lattice(wt,epsik,Lk,D)
     zeta(:) = wr(:) + xmu - sigma(1,:)
     do ik=1,Lk
        do i=1,L
           dzeta1 = dreal(zeta(i))
           dzeta2 = dreal(zeta(L+1-i))
           ddet(i) = (dzeta1-epsik(ik))*(dzeta2-epsik(ik))+dreal(sigma(2,L+1-i))*dreal(Sigma(2,i))
        enddo
        call init_finter(det_finter,wr,ddet,5)
        ipoles(ik) = fzero_brentq(det_poles,0.d0,wr(L))
     enddo
     call splot("poles.last",epsik,ipoles)
     stop
  endif



  if(thermo)then
     !GET Matsubara functions (not working well: aka only at large enough temperature, 
     !because of the finite broadening you can not resolve infinitely small energies/temp)
     converged=.true.
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
     delta=-u*sum(dimag(fg(2,:))*fermi(wr,beta))*fmesh/pi
     n    =-sum(dimag(fg(1,:))*fermi(wr,beta))*fmesh/pi
     allocate(wrx(Lmts),fgm(2,Lmts),sm(2,Lmts))
     wrx = pi/beta*real(2*arange(1,Lmts)-1,8)
     write(*,*)"Get Matsubara GF:",Lmts
     print*,'G(iw)'
     call get_matsubara_gf_from_dos(wr,fg(1,:),fgm(1,:),beta)
     print*,'Sigma(iw)'
     call get_matsubara_gf_from_dos(wr,sigma(1,:),sm(1,:),beta)
     print*,'F(iw)'
     call get_matsubara_gf_from_dos(wr,fg(2,:),fgm(2,:),beta)
     print*,'S(iw)'
     call get_matsubara_gf_from_dos(wr,sigma(2,:),sm(2,:),beta)
     fgm(2,:)=dreal(fgm(2,:))
     sm(2,:)=dreal(sm(2,:))-delta
     call splot("G_iw.last",wrx,fgm(1,:))
     call splot("F_iw.last",wrx,fgm(2,:))
     call splot("Sigma_iw.last",wrx,sm(1,:))
     call splot("Self_iw.last",wrx,sm(2,:))
     Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk),nk(Lk))
     call bethe_lattice(wt,epsik,Lk,D)
     call get_sc_internal_energy(Lmts,wrx,fgm,sm)
     stop
  endif
  ! if(thermo)then
  !    !GET Matsubara functions (not working well: aka only at large enough temperature, 
  !    !because of the finite broadening you can not resolve infinitely small energies/temp)
  !    allocate(wrx(Lmts),fgm(2,Lmts),sm(2,Lmts))
  !    wrx = pi/beta*real(2*arange(1,Lmts)-1,8)
  !    write(*,*)"Get Matsubara GF:"
  !    print*,'G(iw)'
  !    call get_matsubara_gf_from_dos(wr,fg(1,:),fgm(1,:),beta)
  !    print*,'Sigma(iw)'
  !    call get_matsubara_gf_from_dos(wr,sigma(1,:),sm(1,:),beta)
  !    print*,'F(iw)'
  !    call get_matsubara_gf_from_dos(wr,fg(2,:),fgm(2,:),beta)
  !    print*,'S(iw)'
  !    call get_matsubara_gf_from_dos(wr,sigma(2,:),sm(2,:),beta)
  !    fgm(2,:)=dreal(fgm(2,:))
  !    sm(2,:)=dreal(sm(2,:))-delta
  !    call splot("G_iw.last",wrx,fgm(1,:))
  !    call splot("F_iw.last",wrx,fgm(2,:))
  !    call splot("Sigma_iw.last",wrx,sm(1,:))
  !    call splot("Self_iw.last",wrx,sm(2,:))
  !    Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk),nk(Lk))
  !    call bethe_lattice(wt,epsik,Lk,D)
  !    call get_sc_internal_energy(Lmts,wrx,fgm,sm)
  ! endif




  if(octype)then
     converged=.true.
     if(mod(Lf,2)/=0)Lf=Lf+1
     allocate(wrx(Lf),dummy(2,Lf))
     wrx = linspace(-wxmax,wxmax,Lf)
     call cubic_spline(wr,sigma(1,:),wrx,dummy(1,:))
     call cubic_spline(wr,sigma(2,:),wrx,dummy(2,:))
     call splot("Sigma_realw_oc.last",wrx,dummy(1,:))
     call splot("Self_realw_oc.last",wrx,dummy(2,:))
     Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
     call bethe_lattice(wt,epsik,Lk,D)
     call get_sc_optical_conductivity(Lf,wrx,dummy)  
     deallocate(wt,epsik,wrx,dummy)
     stop
  endif




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
     delta=-u*trapz(fmesh,dimag(fg(2,:))*istep(wr))/pi
     n=-trapz(fmesh,dimag(fg(1,:))*istep(wr))/pi

     !GET THE WEISS FIELD \calG0^-1(w)
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


     fg0k(1)=zero ; fg0k(2)=zero
     do i=1,L
        A = -dimag(calG(1,i))/pi
        B = -dimag(calG(2,i))/pi
        nfL= istep(wr(i)-vbias/2.d0)!fermi(wr(i)-vbias/2.d0,beta)
        nfR= istep(wr(i)+vbias/2.d0)!fermi(wr(i)+vbias/2.d0,beta)
        nf= (nfL+nfR)/2.d0
        fg0k(1)%less%w(i) = pi2*xi*nf*A
        fg0k(1)%gtr%w(i)  = pi2*xi*(nf-1.d0)*A
        fg0k(2)%less%w(i) = pi2*xi*nf*B       
        fg0k(2)%gtr%w(i)  = pi2*xi*(nf-1.d0)*B
     enddo

     calG11%less%w = fg0k(1)%less%w    
     calG11%gtr%w = fg0k(1)%gtr%w    
     forall(i=1:L)calG22%less%w(i) = -(fg0k(1)%gtr%w(L+1-i))
     forall(i=1:L)calG22%gtr%w(i) = -(fg0k(1)%less%w(L+1-i))

     calF12%less%w   = -conjg(fg0k(2)%less%w)
     calF12%gtr%w   =  fg0k(2)%gtr%w

     calF21%less%w   = fg0k(2)%less%w
     calF21%gtr%w   = -conjg(fg0k(2)%gtr%w)

     call fftgf_rw2rt(calG11%less%w,calG11%less%t,Lm) ; calG11%less%t=fmesh/pi2*calG11%less%t
     call fftgf_rw2rt(calG11%gtr%w,calG11%gtr%t,Lm)   ; calG11%gtr%t =fmesh/pi2*calG11%gtr%t
     call fftgf_rw2rt(calG22%less%w,calG22%less%t,Lm) ; calG22%less%t=fmesh/pi2*calG22%less%t
     call fftgf_rw2rt(calG22%gtr%w,calG22%gtr%t,Lm)   ; calG22%gtr%t =fmesh/pi2*calG22%gtr%t
     call fftgf_rw2rt(calF12%less%w,calF12%less%t,Lm) ; calF12%less%t=fmesh/pi2*calF12%less%t
     call fftgf_rw2rt(calF12%gtr%w,calF12%gtr%t,Lm)   ; calF12%gtr%t =fmesh/pi2*calF12%gtr%t
     call fftgf_rw2rt(calF21%less%w,calF21%less%t,Lm) ; calF21%less%t=fmesh/pi2*calF21%less%t
     call fftgf_rw2rt(calF21%gtr%w,calF21%gtr%t,Lm)   ; calF21%gtr%t =fmesh/pi2*calF21%gtr%t

     do i=-Lm,Lm 
        sk(1)%less%t(i) = U**2*(calG11%less%t(i)*calG22%less%t(i) - calF12%less%t(i)*calF21%less%t(i))*calG22%gtr%t(-i)
        sk(1)%gtr%t(i) =  U**2*(calG11%gtr%t(i)*calG22%gtr%t(i) - calF12%gtr%t(i)*calF21%gtr%t(i))*calG22%less%t(-i)
        !
        sk(2)%less%t(i) = U**2*(calF12%less%t(i)*calF21%less%t(i) - calG11%less%t(i)*calG22%less%t(i))*calF12%gtr%t(-i)
        sk(2)%gtr%t(i) =  U**2*(calF12%gtr%t(i)*calF21%gtr%t(i)  - calG11%gtr%t(i)*calG22%gtr%t(i))*calF12%less%t(-i)
        !
        sk(1)%ret%t(i) =heaviside(t(i))*(sk(1)%gtr%t(i)-sk(1)%less%t(i))
        sk(2)%ret%t(i) =heaviside(t(i))*(sk(2)%gtr%t(i)-sk(2)%less%t(i))
     enddo
     if(heaviside(0.d0)==1.d0)sk(1)%ret%t(0)=sk(1)%ret%t(0)/2.d0 
     if(heaviside(0.d0)==1.d0)sk(2)%ret%t(0)=sk(2)%ret%t(0)/2.d0
     call fftgf_rt2rw(sk(1)%ret%t,sk(1)%ret%w,Lm) ;      sk(1)%ret%w=dt*sk(1)%ret%w
     call fftgf_rt2rw(sk(2)%ret%t,sk(2)%ret%w,Lm) ;      sk(2)%ret%w=dt*sk(2)%ret%w

     sigma(1,:) = sk(1)%ret%w
     sigma(2,:) = -delta + sk(2)%ret%w

     write(*,"(2f14.9)",advance="no")2.d0*n,delta
     converged = check_convergence(sigma(1,:)+sigma(2,:),eps=eps_error,N1=Nsuccess,N2=Nloop)
     if(printf)call splot("observables.ipt",vbias,delta,xmu,u,n,beta,dble(iloop),append=.true.)
  enddo

  call splot("observables.last",vbias,delta,xmu,u,n,beta,dble(iloop))
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
  ! Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  ! call bethe_lattice(wt,epsik,Lk,D)
  ! call get_sc_optical_conductivity(Lm,wrx,dummy)  
  ! deallocate(wt,epsik)
  deallocate(wrx,dummy)







  !##################################################################
  !##################################################################
contains 
  !##################################################################
  !##################################################################


  function det_poles(w) result(det)
    real(8) :: w
    real(8) :: det
    det = finter(det_finter,w)
  end function det_poles


  elemental function istep(x) result(out)
    real(8),intent(in) :: x
    real(8)            :: out
    if(x < 0.d0) then
       out = 1.0d0
    elseif(x==0.d0)then
       out = 0.50d0
    else
       out = 0.0d0
    endif
  end function istep


  subroutine get_initial_sigma()
    logical                :: check1,check2,check
    inquire(file="Sigma_realw.restart",exist=check1)
    inquire(file="Self_realw.restart",exist=check2)
    check=check1*check2
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


  subroutine get_sc_optical_conductivity(L,wr,sigma)
    integer                :: L
    real(8)                :: wr(L),D,dw
    complex(8)             :: sigma(2,L)
    integer                :: i,ik,iv,iw,Nw
    real(8)                :: vel2,dos,Dfermi,A2,B2,ock
    complex(8)             :: det,zeta1,zeta2,fg(2)
    real(8),allocatable    :: oc(:), Ak(:,:,:),cDOS(:,:),ocw(:)
    complex(8),allocatable :: zeta(:)

    print*,"Get OC with:",L,"freq."
    Nw=L/2
    allocate(Ak(2,Lk,L),zeta(L))

    dw = abs(wr(2)-wr(1))
    allocate(cDOS(2,L))
    cDOS=0.d0
    zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
    do i=1,L
       zeta1 = zeta(i)
       zeta2 = conjg(zeta(L+1-i))
       do ik=1,Lk
          det = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + conjg(sigma(2,L+1-i))*sigma(2,i)
          fg(1)=(zeta2-epsik(ik))/det
          fg(2)=conjg(sigma(2,L+1-i))/det
          Ak(1,ik,i)=-dimag(fg(1))/pi
          Ak(2,ik,i)=-dimag(fg(2))/pi
          cDOS(1,i)=cDOS(1,i)+Ak(1,ik,i)*wt(ik)
          cDOS(2,i)=cDOS(2,i)+Ak(2,ik,i)*wt(ik)
       enddo
    enddo
    call splot("ocDOS.last",wr,cDOS(1,:),cDOS(2,:))
    deallocate(cDOS)
    ! do ik=1,Lk
    !    do i=1,L
    !       write(100,*)wr(i),Ak(1,ik,i)
    !       write(200,*)wr(i),Ak(2,ik,i)
    !    enddo
    !    write(100,*)""
    !    write(200,*)""
    ! enddo

    !Changing the loop order does not affect the calculation.
    D=2.d0*ts
    allocate(oc(Nw),ocw(L))
    oc=0.d0
    call start_progress
    do iv=1,Nw
       ocw   = 0.d0
       do iw=1,L-iv
          Dfermi  = istep(wr(iw)) - istep(wr(iw+iv))
          ock=0.d0
          do ik=1,Lk
             A2 = Ak(1,ik,iw)*Ak(1,ik,iw+iv)
             B2 = Ak(2,ik,iw)*Ak(2,ik,iw+iv)
             vel2= (D**2-epsik(ik)**2)/3.d0
             ock = ock + vel2*(A2-B2)*wt(ik)
          enddo
          ocw(iw) = Dfermi*ock
       enddo
       oc(iv)=trapz(dw,ocw(:L-iv))/wr(Nw+iv)
       call progress(iv,Nw)
    enddo
    call stop_progress
    call splot("OC_realw.ipt",wr(Nw+1:L),oc(:))
    call splot("OC_integral.ipt",vbias,trapz(dw,oc))
  end subroutine get_sc_optical_conductivity





  subroutine get_sc_internal_energy(L,wm,fg,sigma)
    integer    :: L
    real(8)    :: wm(L)
    complex(8) :: fg(2,L),sigma(2,L)
    integer    :: i,ik
    real(8)    :: matssum,fmatssum,checkP,checkdens,vertex,Dssum
    complex(8) :: iw,gkw,fkw,g0kw,f0kw
    real(8)    :: Epot,Etot,Eint,kin,kinsim,Ds,docc
    real(8)    :: Sigma_infty,S_infty,det,det_infty,csi,Ei,thermal_factor
    real(8)    :: free(Lk),Ffree(Lk),n_k(Lk)

    !Get asymptotic self-energies
    Sigma_infty =   dreal(sigma(1,L))
    S_infty     =   dreal(sigma(2,L))

    checkP=0.d0 ; checkdens=0.d0 ;          ! test variables
    call start_progress
    kin=0.d0                      ! kinetic energy (generic)
    Ds=0.d0                       ! superfluid stiffness (Bethe)

    do ik=1,Lk

       csi            = epsik(ik)-(xmu-Sigma_infty)
       Ei             = dsqrt(csi**2 + S_infty**2)
       thermal_factor = dtanh(0.5d0*beta*Ei)
       free(ik)        = 0.5d0*(1.d0 - csi/Ei)*thermal_factor
       Ffree(ik)       =-(0.5d0*S_infty)/Ei*thermal_factor

       fmatssum= 0.d0
       matssum = 0.d0
       Dssum   = 0.d0

       vertex=(4.d0*ts**2-epsik(ik)**2)/3.d0

       do i=1,L
          iw       = xi*wm(i)
          det      = abs(iw+xmu-epsik(ik)-sigma(1,i))**2 + real(sigma(2,i),8)**2
          det_infty= wm(i)**2 + (epsik(ik)-(xmu-Sigma_infty))**2 + S_infty**2

          gkw = (-iw+xmu - epsik(ik) - conjg(sigma(1,i)) )/det
          fkw = -sigma(2,i)/det

          g0kw= (-iw - (epsik(ik)-(xmu-Sigma_infty)))/det_infty
          f0kw=-S_infty/det_infty

          matssum =  matssum +  real(gkw,8)-real(g0kw,8)
          !        matssum =matssum + real(gkw,8)        ! without tails corrections

          fmatssum= fmatssum +  real(fkw,8)-real(f0kw,8)
          Dssum   = Dssum    +  fkw*fkw

       enddo

       n_k(ik)   = 4.d0/beta*matssum + 2.d0*free(ik)
       !    n_k(ik)   = 4.d0/beta*matssum + 1.d0     ! without tails corrections
       checkP    = checkP    - wt(ik)*(2.d0/Beta*fmatssum+Ffree(ik))
       !    print*,checkP,Ffree(ik)
       checkdens = checkdens + wt(ik)*n_k(ik)
       kin    = kin    + wt(ik)*n_k(ik)*epsik(ik)
       Ds=Ds + 8.d0/beta* wt(ik)*vertex*Dssum
       call progress(ik,Lk)
    enddo
    call stop_progress

    kinsim=0
    kinsim = sum(fg(1,:)*fg(1,:)+conjg(fg(1,:))*conjg(fg(1,:))-2.d0*fg(2,:)*fg(2,:))*2.d0*ts**2/beta

    Epot=zero
    Epot = sum(fg(1,:)*sigma(1,:) + fg(2,:)*sigma(2,:))/beta*2.d0

    docc = 0.5d0*n**2
    if(u > 0.01d0)docc=-Epot/u + n - 0.25d0


    Eint=kin+Epot

    Ds=zero
    Ds = sum(fg(2,:)*fg(2,:))/beta*2.d0


    write(*,*)"Asymptotic Self-Energies",Sigma_infty, S_infty
    write(*,*)"n,delta",n,delta
    write(*,*)"Dn% ,Ddelta%",(n-0.5d0*checkdens)/n,(delta + u*checkP)/delta ! u is positive
    write(*,*)'========================================='
    write(*,*)"Kinetic energy",kin
    write(*,*)'========================================='
    write(*,*)"double occupancy   =",docc
    write(*,*)'========================================='
    write(*,*) 'Kinetic Energy TEST (simple formula)'
    write(*,*) '###ACTHUNG: FOR BETHE ONLY####',kinsim
    write(*,*) 'Dkin%',(kin-kinsim)/kin
    write(*,*)'========================================='
    write(*,*) 'Superfluid stiffness',Ds
    write(*,*) 'Potential Energy U(n_up-1/2)(n_do-1/2)',Epot
    write(*,*) 'Internal Energy',Eint
    write(*,*)'========================================='
    call splot("nk_distribution.ipt",epsik,n_k/2.d0,free)
    open(100,file="columns.ipt")
    write(100,"(11A21)")"1vbias","2u","3beta","4n","5kin","6docc","7Ds","8Epot","9Eint"
    close(100)
    open(200,file="thermodynamics.ipt")
    write(200,"(11F21.12)")vbias,u,beta,n,kinsim,docc,Ds,Epot,Eint
    close(200)
    return 
  end subroutine get_sc_internal_energy




end program hmipt
