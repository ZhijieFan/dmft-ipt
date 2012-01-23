!########################################################
!     Program  : PAMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the PAM using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
!     AUTHORS  : Adriano Amaricci
!########################################################
program pammpt
  USE VARS_GLOBAL
  USE MPT_FUNX_SOPT

  implicit none
  real(8)                             :: x(2)
  logical                             :: check
  complex(8)                          :: zeta,alpha
  real(8)                             :: npp,gzero,gmu,ntot,shift
  complex(8),dimension(:),allocatable :: fgp,fgp0,sigmap

  call read_input("inputIPT.in")
  call alloc_memory()
  allocate(fgp(-L:L),fgp0(-L:L),sigmap(-L:L))

  gmu=xmu  ; gzero=0.d0
  if((ed0-ep0) > 0.d0)gzero=0.5*(ep0+ed0+sqrt((ep0-ed0)**2 + 4*Vpd**2))
  if((ed0-ep0) < 0.d0)gzero=0.5*(ep0+ed0-sqrt((ep0-ed0)**2 + 4*Vpd**2))
  if((ed0-ep0) /=0.d0)xmu=gmu+gzero !true ED chemical potential
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  sigma=zero ; sigma=U*(n-0.5d0) ; xmu0=xmu ; n0=nobj ; n=nobj
  do iloop=1,nloop
     write(*,"(A,i5,A,i5)")"DMFT-loop",iloop,"/",nloop

     do i=-L,L
        alpha  = cmplx(wr(i),eps)  +xmu -ed0 - sigma(i)
        sigmap(i) = Vpd**2/alpha
        zeta  = cmplx(wr(i),eps)  +xmu -ep0 - sigmap(i)
        fgp(i) = gfbether(wr(i),zeta,D)
        fg(i) = one/alpha + Vpd**2/alpha**2*fgp(i)
     enddo
     n=sum(fermi(wr,beta)*aimag(fg))/sum(aimag(fg))

     npp = 2.d0*sum(aimag(fgp(-L:L))*fermi(wr(-L:L),beta))/sum(aimag(fgp(-L:L)))
     ntot = npp+2.d0*n

     !Get the hybridization functions: \Gamma(w+xmu)
     gamma= cmplx(wr,eps) + xmu - ed0 - sigma - one/fg

     !Fix the xmu0 w/ condition n0=n
     x(1)=xmu
     x(2)=xmu0
     call broydn(x,check)
     ! xmu=x(1)
     ! xmu0=x(2)

     call solve_spt_finiteT(iloop)
     call splot("npVSiloop.ipt",iloop,npp,append=TT)
     call splot("ntotVSiloop.ipt",iloop,ntot,append=TT)
  enddo

  call splot("npVS"//extension,xout,npp,append=TT)
  call splot("ntotVS"//extension,xout,ntot,append=TT)
  call splot("DOSp.ipt",wr,-aimag(fgp)/pi,append=TT)
  call splot("Sigmap_realw.ipt",wr,sigmap)
  call splot("Gp_realw.ipt",wr,fgp)
  Lk=(Nx+1)**2
  call pam_getenergy_spt(Lk)

contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  subroutine pam_getenergy_spt(Lk)    
    integer   :: Lk,Lm
    real(8)   :: nd,np,docc,Etot,Ekin,Epot,Ehyb,Eepsi,Emu,nk(Lk),npk(Lk)
    real(8)   :: epsi,de
    complex(8),allocatable :: gf(:),sf(:),gff(:)
    complex(8),allocatable :: gp(:),sp(:),gpp(:)
    real(8),allocatable    :: gtau(:),gptau(:)
    complex(8) :: gamma,alpha

    Lm=int(L*beta/pi);if(beta<1)Lm=L;if(Lm<L)Lm=L
    Lm=2**14
    allocate(gf(Lm),sf(Lm),gtau(0:L),gff(L))
    allocate(gp(Lm),sp(Lm),gptau(0:L),gpp(L))

    call getGmats(wr,fg,gf,beta)
    call getGmats(wr,fgp,gp,beta)
    call getGmats(wr,sigma,sf,beta)
    call getGmats(wr,sigmap,sp,beta)
    call splot("Giw.ipt",(/(pi/beta*dble(2*i-1),i=1,Lm)/),gf)
    call splot("Gpiw.ipt",(/(pi/beta*dble(2*i-1),i=1,Lm)/),gp)

    call cfft_iw2tau(gf,gtau,beta)
    call cfft_iw2tau(gp,gptau,beta)
    call splot("Gtau.ipt",(/(dble(i)*beta/dble(L),i=0,L)/),gtau(0:L))
    call splot("Gptau.ipt",(/(dble(i)*beta/dble(L),i=0,L)/),gptau(0:L))

    nd=-2.d0*gtau(L)
    np=-2.d0*gptau(L)

    !Energy && k-dispersed quantities
    Ekin=0.d0
    de=2.d0*D/dble(Lk)
    do ik=1,Lk
       epsi=-D + de*dble(ik)
       do i=1,L
          w=pi/beta*dble(2*i-1)
          gamma=xi*w+xmu-ed0-sf(i)
          alpha=xi*w+xmu-ep0-epsi
          gff(i)=one/(gamma - vpd**2/alpha)
          gpp(i)=one/(alpha-sp(i))
       enddo
       call cfft_iw2tau(gff,gtau,beta)
       call cfft_iw2tau(gpp,gptau,beta)
       nk(ik)=-gtau(L)
       npk(ik)=-gptau(L)
       Ekin=Ekin + nk(ik)*epsi/dble(Lk)
    enddo

    Epot=dot_product(conjg(sf(:)),gf(:))
    Ehyb=4.d0*dot_product(conjg(sp(:)),gp(:))


    Ehyb=2.d0*Ehyb/beta
    Epot=2.d0*Epot/beta

    docc=0.5*nd  - 0.25d0
    if(u>=1.d-2)docc = Epot/U + 0.5*nd - 0.25d0

    Eepsi=ed0*(nd-1.d0) + ep0*(np-1.d0)
    Emu=-xmu*(ntot-2.d0)

    Etot=Ekin+Epot+Ehyb+Eepsi+Emu

    write(*,"(A,f13.9)")"docc   =",docc
    write(*,"(A,f13.9)")"nd     =",nd
    write(*,"(A,f13.9)")"np     =",np
    write(*,"(A,f13.9)")"ntot   =",np+nd
    write(*,"(A,f13.9)")"Etot   =",Etot
    write(*,"(A,f13.9)")"Ekin   =",Ekin
    write(*,"(A,f13.9)")"Epot   =",Epot
    write(*,"(A,f13.9)")"Ehyb   =",Ehyb
    write(*,"(A,f13.9)")"Eepsi  =",Eepsi
    write(*,"(A,f13.9)")"Emu    =",Emu

    call splot("nkVSepsik.ipt",(/(-D + de*dble(ik),ik=1,Lk )/),nk(1:Lk),npk(1:Lk))
    call splot("EtotVS"//trim(extension),xout,Ekin+Epot,append=TT)
    call splot("EkinVS"//trim(extension),xout,Ekin,append=TT)
    call splot("EpotVS"//trim(extension),xout,Epot,append=TT)
    call splot("EhybVS"//trim(extension),xout,Ehyb,append=TT)
    call splot("EmuVS"//trim(extension),xout,Emu,append=TT)
    call splot("EepsiVS"//trim(extension),xout,Eepsi,append=TT)
    call splot("doccVS"//trim(extension),xout,docc,append=TT)
    deallocate(gf,sf,gtau,gff)
    deallocate(gp,sp,gptau,gpp)
  end subroutine pam_getenergy_spt
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************

end program pammpt



function funcv(x)
  USE VARS_GLOBAL
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  complex(8)                       :: alpha,zeta,sigpp,gpp
  xmu=x(1)
  do i=-L,L
     alpha  = cmplx(wr(i),eps)  +xmu -ed0 - sigma(i)
     sigpp  = Vpd**2/alpha
     zeta   = cmplx(wr(i),eps)  +xmu -ep0 - sigpp
     gpp    = gfbether(wr(i),zeta,D)
     fg(i) = one/alpha + Vpd**2/alpha**2*gpp
  enddo
  n=sum(fermi(wr,beta)*aimag(fg))/sum(aimag(fg))


  xmu0=x(2)
  fg0 = one/(cmplx(wr,eps) + xmu0 - ed0 - gamma -U*(n-0.5d0))
  n0=sum(fermi(wr,beta)*aimag(fg0))/sum(aimag(fg0))

  funcv(1)=nobj-n
  funcv(2)=nobj-n0
  write(*,"(4(f13.9))")n,n0,xmu,xmu0
end function funcv



