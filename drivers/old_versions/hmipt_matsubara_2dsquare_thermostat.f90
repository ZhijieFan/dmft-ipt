
!########################################################
!     Program  : HMIPT
!     TYPE     : Main program
!     PURPOSE  : Solve the Hubbard model using DMFT-IPT
!     AUTHORS  : Adriano Amaricci
!########################################################
program hmipt_matsuara_2dsquare
  USE DMFT_IPT
  USE SQUARE_LATTICE
  USE IOTOOLS
  implicit none

  logical                :: converged
  real(8)                :: n,z,Wbath,Vbath
  integer                :: i,Lk
  complex(8)             :: zeta
  type(matsubara_gf)     :: fg,sigma
  complex(8),allocatable :: fg0(:)
  real(8),allocatable    :: wm(:),tau(:),wt(:),epsik(:),nk(:)
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
  allocate(wm(L),tau(0:L))
  call allocate_gf(fg,L)
  call allocate_gf(sigma,L)
  allocate(bath_sigma(L))
  allocate(fg0(L))

  !build freq. array
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)

  !build square lattice structure:
  Lk   = square_lattice_dimension(Nx)
  allocate(wt(Lk),epsik(Lk),nk(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)

  !build bath
  call get_bath

  !dmft loop:
  sigma%iw=zero ; iloop=0 ; converged=.false.
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     do i=1,L
        zeta = xi*wm(i) - sigma%iw(i) -bath_sigma(i)
        fg%iw(i) = sum_overk_zeta(zeta,epsik,wt)
     enddo
     call fftgf_iw2tau(fg%iw,fg%tau,beta)
     n   = -fg%tau(L)
     fg0 = one/(one/fg%iw + sigma%iw)
     sigma%iw= solve_ipt_matsubara(fg0)
     converged=check_convergence(sigma%iw,eps_error,nsuccess,nloop)
     z=1.d0 - dimag(sigma%iw(1))/wm(1);z=1.d0/z
     call splot("nVSiloop.ipt",iloop,n,append=TT)
     call splot("zetaVSiloop.ipt",iloop,z,append=TT)
  enddo
  call close_file("nVSiloop.ipt")
  call close_file("zetaVSiloop.ipt")
  call splot("G_iw.ipt",wm,fg%iw,append=printf)
  call splot("G0_iw.ipt",wm,fg0,append=printf)
  call splot("Sigma_iw.ipt",wm,sigma%iw,append=printf)
  nk = square_lattice_momentum_distribution(Lk)
  call splot("nkVSepsk.ipt",epsik,nk,append=printf)

contains

  function square_lattice_momentum_distribution(Lk) result(nk)
    integer            :: Lk
    integer            :: ik,i
    real(8)            :: nk(Lk),docc,Epot

    Epot=dot_product(conjg(sigma%iw(:)),fg%iw(:))*2.d0/beta
    docc=0.25d0
    if(u>=1.d-3)docc = Epot/u + 0.25d0
    call splot("docc.ipt",beta,u,vbath**2/2.d0/wbath,epot,docc,append=printf)

    do ik=1,Lk
       fg%iw=one/(xi*wm - epsik(ik) - sigma%iw)
       call fftgf_iw2tau(fg%iw,fg%tau,beta)
       nk(ik)=-fg%tau(L)
    enddo
  end function square_lattice_momentum_distribution


  subroutine get_bath
    complex(8) :: peso,bath_gf
    real(8)    :: w,wini,wfin,dw,de,nless,ngtr,en,epsimu(Lk)
    integer    :: i,iw,im
    epsimu = linspace(-2.d0*wbath,2.d0*wbath,Lk,mesh=de)

    bath_sigma=zero
    do i=1,L
       bath_gf=zero
       do im=1,Lk
          bath_sigma(i)=bath_sigma(i)+vbath**2/(xi*wm(i)-epsimu(im))/dble(Lk)
       enddo
    enddo
    call splot("bathSigma.ipt",wm,bath_sigma)
  end subroutine get_bath

end program hmipt_matsuara_2dsquare
