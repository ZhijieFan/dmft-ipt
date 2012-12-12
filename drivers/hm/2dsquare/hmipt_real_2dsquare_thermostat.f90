program hmipt_2dsquare
  !########################################################
  !     Program  : HMIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Hubbard model using DMFT-IPT
  !     AUTHORS  : Adriano Amaricci
  !########################################################
  !LOCAL:
  USE DMFT_IPT
  USE SQUARE_LATTICE
  USE IOTOOLS
  implicit none

  logical                :: converged
  real(8)                :: n,Wbath,Vbath
  integer                :: i,Lk
  complex(8)             :: zeta
  complex(8),allocatable :: sigma(:),fg(:),fg0(:)
  real(8),allocatable    :: wr(:),wt(:),epsik(:),nk(:)
  !
  real(8),allocatable    :: bath_dens(:)
  complex(8),allocatable :: bath_sigma(:)

  call read_input("inputIPT.in")

  wbath=20.d0
  Vbath=0.d0
  do i=1,command_argument_count()
     nml_var=get_cmd_variable(i)
     if(nml_var%name=="VBATH")read(nml_var%value,*)vbath
     if(nml_var%name=="WBATH")read(nml_var%value,*)wbath
  end do


  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(wr(L))
  allocate(bath_dens(L),bath_sigma(L))

  !build freq. array
  wr = linspace(-wmax,wmax,L,mesh=fmesh)



  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  !build bath
  call get_bath

  !dmft loop:
  sigma=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta = cmplx(wr(i),eps) - sigma(i) - bath_sigma(i)
        fg(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo
     n   = sum(aimag(fg)*fermi(wr,beta))/sum(aimag(fg))
     fg0 = one/(one/fg + sigma)
     sigma= solve_ipt_sopt(fg0,wr)
     converged=check_convergence(sigma,eps_error,nsuccess,nloop)
     call splot("nVSiloop.ipt",iloop,n,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  call splot("DOS.ipt",wr,-aimag(fg)/pi,append=printf)
  call splot("G_realw.ipt",wr,fg,append=printf)
  call splot("G0_realw.ipt",wr,fg0,append=printf)
  call splot("Sigma_realw.ipt",wr,sigma,append=printf)
  nk = square_lattice_momentum_distribution(Lk)
  call splot("nkVSepsk.ipt",epsik,nk,append=printf)
contains

  function square_lattice_momentum_distribution(Lk) result(nk)
    integer            :: Lk
    integer,parameter  :: M=2048
    integer            :: ik,i
    type(matsubara_gf) :: gm,sm
    real(8)            :: nk(Lk),wm(M),w,Epot,docc,nd
    call allocate_gf(gm,M)
    call allocate_gf(sm,M)
    wm   = pi/beta*real(2*arange(1,M)-1,8)
    call get_matsubara_gf_from_DOS(wr,sigma,sm%iw,beta)
    do ik=1,Lk
       gm%iw=one/(xi*wm - epsik(ik) - sm%iw)
       call fftgf_iw2tau(gm%iw,gm%tau,beta)
       nk(ik)=-gm%tau(M)
    enddo
    call get_matsubara_gf_from_DOS(wr,fg,gm%iw,beta)
    call splot("Sigma_iw.ipt",wm,sm%iw)
    call splot("G_iw.ipt",wm,gm%iw)

    Epot=0.d0
    do i=1,M
       Epot=Epot+sm%iw(i)*gm%iw(i)
       !Epot=dot_product(conjg(sm%iw(:)),gm%iw(:))*2.d0/beta
    enddo
    Epot=Epot/beta*2.d0

    docc=0.25d0
    if(u>=1.d-3)docc = Epot/u + 0.25d0
    call splot("docc.ipt",temp,u,vbath**2/2.d0/wbath,epot,docc,append=printf)
  end function square_lattice_momentum_distribution

  subroutine get_bath
    complex(8) :: peso,bath_gf
    real(8)    :: w,wini,wfin,dw,de,nless,ngtr,en,epsimu(Lk),wfreq(L)
    integer    :: i,iw,im

    ! dw=2.d0*wbath/dble(Lk)
    ! do i=1,Lk
    !    epsimu(i)=-wbath+dble(i)*dw
    !    print*,epsimu(i)
    ! enddo
    epsimu = linspace(-2.d0*wbath,2.d0*wbath,Lk,mesh=de)
    wfreq  = linspace(-2.d0*wbath,2.d0*wbath,L,mesh=dw)

    bath_sigma=zero
    do i=1,L
       bath_gf=zero
       do im=1,Lk
          bath_gf=bath_gf+1.d0/(wr(i)+xi*eps-epsimu(im))/dble(Lk)
       enddo
       bath_sigma=vbath**2*bath_gf
       bath_dens=-aimag(bath_gf)/pi
    enddo

    ! do i=1,L
    !    w=wr(i)
    !    bath_dens(i) = heaviside(Wbath-abs(w))/(2.d0*Wbath)
    !    bath_gf      = log(abs((Wbath+w)/(Wbath-w)))-xi*pi*bath_dens(i)
    !    bath_sigma(i)= Vbath**2*bath_gf
    ! enddo
    call splot("bathDens.ipt",wfreq,bath_dens)
    call splot("bathSigma.ipt",wfreq,bath_sigma)
  end subroutine get_bath

end program hmipt_2dsquare
