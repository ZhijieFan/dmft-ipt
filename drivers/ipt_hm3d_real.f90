program hmipt_matsubara
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none
  logical                :: converged,check
  real(8)                :: wmix,ts
  integer                :: i,ik,ie,iloop,L,Lk,Le,unit,Loc,Nx
  complex(8)             :: zeta,gf
  complex(8),allocatable :: fg(:),fg0(:),sigma(:),fg0_prev(:)
  real(8),allocatable    :: wreal(:),wtk(:),edos(:),Ndos(:)
  real(8)                :: wmin,fmesh,wband,de
  character(len=24)      :: finput
  complex(8),allocatable :: Ek(:,:,:)
  !


  call parse_cmd_variable(finput,"finput",default="inputIPT.conf")
  call parse_input_variable(ts,"ts",finput,default=1d0/6d0)
  call parse_input_variable(L,"L",finput,default=2000)
  call parse_input_variable(Nx,"Nx",finput,default=10)
  call parse_input_variable(Le,"Le",finput,default=1000)
  call parse_input_variable(wmix,"WMIX",finput,default=0.5d0)
  call read_input(finput)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(1,"norb")
  call add_ctrl_var(1,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  wmin = -wmax
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(0.01d0,"eps")


  Lk = Nx*Nx*Nx
  call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
  allocate(Ek(1,1,Lk))
  allocate(Wtk(Lk))
  call TB_build_model(Ek,hk_model,1,[Nx,Nx,Nx])
  Wtk = 1d0/Lk
  ! call get_free_dos(dreal(Ek(1,1,:)),Wtk)
  allocate(Edos(Le))
  allocate(Ndos(Le))
  wband = 6d0*ts
  Edos = linspace(-wband,wband,Le,mesh=de)
  Ndos = 0d0
  do ie=1,Le
     zeta = dcmplx(Edos(ie),eps)
     gf   = sum(wtk/( zeta-dreal(Ek(1,1,:)) ))
     Ndos(ie) = -dimag(gf)/pi
  enddo
  Ndos = Ndos/trapz(Ndos,de)
  call splot("DOS3d.ipt",Edos,Ndos)

  !build freq. array
  allocate(wreal(L))
  wreal = linspace(wmin,wmax,L,mesh=fmesh)

  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(fg0_prev(L))

  !get or read first sigma
  sigma=zero
  inquire(file="Sigma.restart",exist=check)
  if(check)call read_array("Sigma.restart",sigma)


  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !
     !SELF-CONSISTENCY:
     fg=zero
     do i=1,L
        zeta = dcmplx(wreal(i),eps) - sigma(i)
        fg(i) = simps(Ndos(:)/(zeta-Edos(:)),de)
     enddo
     !
     fg0_prev=fg0
     fg0 = one/(one/fg + sigma)
     call broyden_mix(fg0_prev,fg0,0.5d0,5,iloop)
     !
     !IMPURITY SOLVER
     sigma = ipt_solve_real(fg0,wreal)
     !
     converged=check_convergence(fg0,dmft_error,nsuccess,nloop)
     !
  enddo


  call splot("Gloc_wreal.ipt",wreal,fg)
  call splot("G0_wreal.ipt",wreal,fg0)
  call splot("Sigma_wreal.ipt",wreal,sigma)
  call save_array("Sigma.restart",sigma)

  call get_optical_conductivity()


contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky,kz
    complex(8)           :: hk(N,N)
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    Hk = -2d0*ts*(cos(kx)+cos(ky)+cos(kz))
  end function hk_model



  subroutine get_optical_conductivity()
    integer                       :: Nw
    real(8),dimension(Le)         :: Phi_xx
    !
    real(8),dimension(-L/2+1:L/2) :: ReOc,Omega,ImOc
    complex(8),dimension(L)       :: Oc,Edie
    real(8),dimension(L)          :: Ref
    real(8)                       :: Fw,Fw_p_v,DeltaF
    real(8)                       :: OCtmp
    real(8)                       :: sumAe
    complex(8)                    :: Zw,Zw_p_v
    real(8)                       :: om,nu,om_p_nu
    real(8)                       :: Ak,Aw,Aw_p_v
    integer                       :: i,iw,iv,ie
    !
    Nw=L/2
    Omega(-Nw+1:Nw:1)=wreal
    !
    Phi_xx= 0d0
    do ie=1,Le
       Phi_xx(ie) = -1d0/3d0*trapz( Edos(:ie)*Ndos(:ie), de)       
    enddo
    !
    call start_timer
    ReOc=0d0
    do iv=1,Nw
       nu = Omega(iv)
       !
       do iw=1,L-iv
          om      = wreal(iw)
          om_p_nu = om+nu
          !
          Fw     = fermi(om,beta)
          Fw_p_v = fermi(om_p_nu,beta)
          DeltaF = (Fw - Fw_p_v)/nu
          if(abs(nu)<1d0/beta)DeltaF=-dfermi(om,beta)
          if(abs(DeltaF)<1d-12)cycle
          !
          zw      = dcmplx(om,eps)-Sigma(iw)
          zw_p_v  = dcmplx(om_p_nu,eps)-Sigma(iw+iv)
          !
          sumAe=0d0
          do ie=1,Le
             Aw      = -dimag(1d0/( Zw     - Edos(ie)) )/pi
             Aw_p_v  = -dimag(1d0/( Zw_p_v - Edos(ie)) )/pi
             sumAe   =  sumAe + Phi_xx(ie)*Aw*Aw_p_v*de
          enddo
          !
          ReOC(iv) = ReOC(iv) + DeltaF*sumAe*fmesh*pi2
       enddo
       call eta(iv,Nw)
    enddo
    call stop_timer
    ReOc(-Nw+1:0) = ReOc(Nw:1:-1)
    ImOc = kronig(reOc,omega,size(reOc))
    call splot("reOC_realw.ipt",Omega,ReOC)
    call splot("imOC_realw.ipt",Omega,ImOC)
    !
    Oc = ReOc + xi*ImOc
    call splot("Oc_realw.ipt",wreal,Oc)
    !
    Edie = 1d0 + xi*4*pi*Oc/wreal
    call splot("Ed_realw.ipt",wreal,Edie)
    call splot("Ed_mod_realw.ipt",wreal,abs(Edie))
    !
    Ref = abs( (sqrt(Edie)-1d0)/(sqrt(Edie)+1d0))**2
    call splot("Ref_realw.ipt",wreal,Ref)
    !
    call splot("Oc_over_Edie2_realw.ipt",wreal,Oc/abs(Edie)**2)
  end subroutine get_optical_conductivity





  !+-------------------------------------------------------------------+
  !PURPOSE  : calculate the Fermi-Dirac distribution
  !+-------------------------------------------------------------------+
  elemental function dfermi(x,beta,limit)
    real(8),intent(in)          :: x, beta
    real(8),optional,intent(in) :: limit
    real(8)                     :: fe,arg,limit_,dfermi
    limit_ = 200d0 ; if(present(limit))limit_=abs(limit)
    arg    = x*beta
    fe     = fermi(x,beta,limit_)
    dfermi = -exp(arg)*fe**2
  end function dfermi



end program hmipt_matsubara
