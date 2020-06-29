program hmipt_matsubara
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none
  logical                :: converged,check
  real(8)                :: wmix,ts
  integer                :: i,ik,iloop,L,Lk,Le,unit,Loc
  complex(8)             :: zeta
  complex(8),allocatable :: fg(:),fg0(:),sigma(:),fg0_prev(:)
  real(8),allocatable    :: wreal(:),edos(:),Ndos(:)
  real(8)                :: wmin,fmesh,wband,de,x
  character(len=24)      :: finput
  type(finter_type)      :: fSigma
  !


  call parse_cmd_variable(finput,"finput",default="inputIPT.conf")
  call parse_input_variable(ts,"ts",finput,default=0.25d0)
  call parse_input_variable(L,"L",finput,default=2000)
  call parse_input_variable(Le,"Le",finput,default=1000)
  call parse_input_variable(wmix,"WMIX",finput,default=0.5d0)
  call read_input(finput)
  call save_input(finput)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(1,"norb")
  call add_ctrl_var(1,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  wmin = -wmax
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(0.01d0,"eps")


  !BUILD THE LATTICE STRUCTURE (use tight_binding):
  allocate(Edos(Le))
  allocate(Ndos(Le))
  wband = 4d0*ts
  Edos = linspace(-wband,wband,Le,mesh=de)
  do ik=1,Le
     Ndos(ik) = dens_2dsquare(edos(ik),ts)
  enddo
  call splot("DOS2d.ipt",Edos,Ndos)


  allocate(wreal(L))
  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(fg0_prev(L))


  !build freq. array
  wreal = linspace(wmin,wmax,L,mesh=fmesh)



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

  call init_finter(fSigma,wreal,Sigma,5)
  call get_optical_conductivity()


contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky
    complex(8)              :: hk(N,N)
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = -one*2d0*ts*(cos(kx)+cos(ky))
  end function hk_model

  function finter_sigma(x) result(sigma)
    real(8)    :: x
    complex(8) :: sigma
    sigma = cinter(FSigma,x)
    if(uloc(1)==0d0)sigma=zero
  end function finter_sigma

  subroutine get_optical_conductivity()
    real(8),dimension(Le)   :: Phi_xx
    real(8),dimension(L/2)    :: ReOc,ImOc!,Omega
    complex(8),dimension(L/2) :: Oc,Edie,ncgs
    real(8),dimension(L/2)    :: Ref
    real(8),dimension(L/2)    :: Omega
    real(8)                 :: Fw,Fw_p_v,DeltaF
    real(8)                 :: OCtmp
    real(8)                 :: sumAe
    complex(8)              :: Zw,Zw_p_v
    real(8)                 :: om,nu,om_p_nu
    real(8)                 :: Ak,Aw,Aw_p_v
    integer                 :: i,iw,iv,ie

    !
    Phi_xx= 0d0
    do ie=1,Le
       Phi_xx(ie) = -1d0/2d0*trapz( Edos(:ie)*Ndos(:ie), de)
    enddo
    !
    Omega = wreal(L/2+1:L)
    !
    call start_timer
    ReOc=0d0
    do iv=1,L/2
       nu = Omega(iv)
       !
       do iw=1,L
          om      = wreal(iw)
          om_p_nu = om+nu
          if(om_p_nu>wmax)cycle
          !
          Fw     = fermi(om,beta)
          Fw_p_v = fermi(om_p_nu,beta)
          DeltaF = (Fw - Fw_p_v)/nu
          if(abs(nu) < 1d0/beta )DeltaF=-dfermi(om,beta)
          !
          if(abs(DeltaF)<1d-9)cycle
          !
          zw      = dcmplx(om,eps)-finter_sigma(om)
          zw_p_v  = dcmplx(om_p_nu,eps)-finter_sigma(om_p_nu)
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
       call progress(iv,size(Omega))
    enddo
    call stop_timer

    ImOc = -one/pi*kronig(reOc,Omega,size(reOc))
    call splot("reOC_realw.ipt",Omega,ReOC)
    call splot("imOC_realw.ipt",Omega,ImOC)

    Oc = ReOc + xi*ImOc
    call splot("Oc_realw.ipt",Omega,Oc)


    Edie = 1d0 + xi*4*pi*Oc/Omega
    call splot("Ed_realw.ipt",Omega,Edie)
    call splot("invEd_realw.ipt",Omega,one/Edie)
    call splot("Ed_mod_realw.ipt",Omega,abs(Edie))
    !
    ncgs= zsqrt(Edie)
    Ref = abs( (one-ncgs)/(ncgs+one) )**2
    call splot("Ref_realw.ipt",Omega,Ref,append=.true.)
    !
    call splot("ReOc_over_Edie2_realw.ipt",Omega,-Omega*dimag(one/Edie))


    

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
