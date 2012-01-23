!########################################################
!     Program  : PAMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the PAM using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
!     AUTHORS  : Adriano Amaricci
!########################################################
program pammpt
  USE VARS_GLOBAL
  USE MPT_FUNX_MATS

  implicit none
  real(8)    :: x(2)
  logical    :: check
  complex(8) :: zeta,alpha
  real(8)    :: npp,gzero,gmu,ntot
  type(matsubara_gf):: fgmp,fg0mp,sigmamp

  call read_input("inputIPT.in")
  call alloc_memory()
  call allocate_gf(fgmp,L)
  call allocate_gf(fg0mp,L)
  call allocate_gf(sigmamp,L)

  xmu=xmu+abs(nobj-0.5d0);print*,xmu
  gmu=xmu ; gzero=0.d0
  if((ed0-ep0) > 0.d0)gzero=0.5*(ep0+ed0+sqrt((ep0-ed0)**2 + 4*Vpd**2))
  if((ed0-ep0) < 0.d0)gzero=0.5*(ep0+ed0-sqrt((ep0-ed0)**2 + 4*Vpd**2))
  if((ed0-ep0) /=0.d0)xmu=gmu+gzero !true ED chemical potential
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  sigmam=zero ; sigmam%iw=U*(n-0.5d0) ; xmu0=xmu ; n0=nobj ; n=nobj
  do iloop=1,nloop
     write(*,"(A,i5,A,i5)")"DMFT-loop",iloop,"/",nloop

     do i=1,L
        alpha  = xi*wm(i)  +xmu -ed0 - sigmam%iw(i)
        sigmamp%iw(i) = Vpd**2/alpha
        zeta  = xi*wm(i)  +xmu -ep0 - sigmamp%iw(i)
        fgmp%iw(i) = gfbether(wm(i),zeta,D)
        fgm%iw(i) = one/alpha + Vpd**2/alpha**2*fgmp%iw(i)
     enddo
     call cfft_iw2tau(fgm%iw,fgm%tau,beta)  ; 
     call cfft_iw2tau(fgmp%iw,fgmp%tau,beta); 
     n=-real(fgm%tau(L))
     npp=-2.d0*real(fgmp%tau(L))
     ntot = npp+2.d0*n

     !Get the hybridization functions: \Gamma(w+xmu)
     gammam%iw= xi*wm + xmu - ed0 - sigmam%iw - one/fgm%iw

     !Fix the xmu0 w/ condition n0=n
     x(1)=xmu
     x(2)=xmu0
     call broydn(x,check)

     call solve_mpt_matsubara(iloop)
     if(printf)then
        call splot("npVSiloop.ipt",iloop,npp,append=TT)
        call splot("ntotVSiloop.ipt",iloop,ntot,append=TT)
     endif
  enddo
  if(printf)then
     call splot("npVS"//extension,xout,npp,append=TT)
     call splot("ntotVS"//extension,xout,ntot,append=TT)
     call splot("Sigmap_iw.ipt",wm,sigmamp%iw)
     call splot("Gp_iw.ipt",wm,fgmp%iw)
     call splot("Gp_tau.ipt",tau,fgmp%tau)
  endif
  if(getenergy)call pam_getener()

contains

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine pam_getener()
    implicit none
    integer :: i,Nf
    real(8) :: depsi,e,free
    real(8) :: w,doble,temp
    real(8) :: znd,znp,zntot,znpd
    real(8) :: Etot,Ekin,Epot,Ehyb,Eepsi,Emu
    real(8) :: kag,kag0,coda1
    complex(8) :: iw,greenp0,greenp,selfp,greend0,greend,selfd,zita,zita0
    complex(8) :: ptail,dtail,gamma

    write(*,"(A)")"Getting the energy"
    temp=1.d0/beta
    free=0.d0
    depsi=2.d0*d/dble(L)
    do i=1,L
       e=-d+dble(i)*depsi
       free=free+depsi*e*dens_bethe(e)*fermi(e,beta)
    enddo

    Nf=L
    write(*,"(A,I10)")"Get En with=",Nf

    Ekin=zero
    Ehyb=zero
    Epot=0.d0
    znd=0.d0
    znp=0.d0
    do i=1,Nf
       w=pi/beta*dble(2*i-1) ; iw=xi*w
       selfp=sigmamp%iw(i)
       selfd=sigmam%iw(i)
       zita0  = iw
       zita   = iw+xmu-ep0-selfp
       greenp0= gfbethe(w,zita0,d)
       greenp = gfbethe(w,zita,d)
       gamma  = iw+xmu-ed0-selfd
       greend = one/gamma + Vpd**2/gamma**2*greenp

       !Kinetic energy     
       Ekin=Ekin+2.d0*(zita*greenp-zita0*greenp0)
       !Hybridization energy
       Ehyb=Ehyb+4.d0*(selfp*greenp)
       ! !Potential Energy second type
       Epot=Epot+(selfd*greend)

       znd=znd+2.d0*temp*real(greend)
       znp=znp+2.d0*temp*real(greenp)
    enddo

    znd=2.d0*znd;znp=2.d0*znp
    zntot=znd+znp

    Ekin=2.d0*Ekin
    Ehyb=2.d0*Ehyb
    Epot=2.d0*Epot

    ! kag0=-real(greenp0)*w
    ! kag=-real(greenp)*w
    ! kag0=abs(kag0)
    ! kag=abs(kag)
    ! coda1=(kag0-kag)*pi+2.d0*(kag*atan(Nf/kag)-kag0*atan(Nf/kag0))
    ! coda1=temp*coda1

    Ekin=temp*Ekin+2.d0*free !+2.d0*coda1
    Epot=temp*Epot
    Ehyb=temp*Ehyb 
    Eepsi=ed0*znd + ep0*znp
    Emu=-xmu*zntot!-2.d0*temp*xmu*zntot

    Etot=Ekin+Epot+Ehyb+Eepsi+Emu

    doble=0.5d0*(znd+1.d0) - 0.25d0
    if(u > 0.01d0)doble=real(Epot)/u + 0.5d0*(znd+1.d0) - 0.25d0

    print*,""
    write(*,"(A,f13.9)")"docc   =",doble
    write(*,"(A,f13.9,f13.9)")"nd     =",znd+1.d0,2.d0*n
    write(*,"(A,f13.9,f13.9)")"np     =",znp+1.d0,npp
    write(*,"(A,f13.9,f13.9)")"ntot   =",znp+znd+2.d0,ntot
    !write(*,"(A,f13.9)")"code   =",coda1

    call splot("EfreeVS"//extension,temp,free,append=TT)
    call splot("EtotVS"//extension,xout,Etot,append=TT)
    call splot("EkinVS"//extension,xout,Ekin,append=TT)
    call splot("EpotVS"//extension,xout,Epot,append=TT)
    call splot("EhybVS"//extension,xout,Ehyb,append=TT)
    call splot("EmuVS"//extension,xout,Emu,append=TT)
    call splot("EepsiVS"//extension,xout,Eepsi,append=TT)
    call splot("EndVS"//extension,xout,znd+1.d0,append=TT)
    call splot("EnpVS"//extension,xout,znp+1.d0,append=TT)
    call splot("EntotVS"//extension,xout,zntot+2.d0,append=TT)
    call splot("EdoccVS"//extension,xout,doble,append=TT)
    return
  end subroutine pam_getener
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************

end program pammpt



function funcv(x)
  USE VARS_GLOBAL
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  complex(8)                       :: alpha,zeta,sigpp,gpp

  xmu=x(1)
  do i=1,L
     alpha  = xi*wm(i)  +xmu -ed0 - sigmam%iw(i)
     sigpp  = Vpd**2/alpha
     zeta   = xi*wm(i)  +xmu -ep0 - sigpp
     gpp    = gfbether(wm(i),zeta,D)
     fgm%iw(i) = one/alpha + Vpd**2/alpha**2*gpp
  enddo
  call cfft_iw2tau(fgm%iw,fgm%tau,beta)  
  n=-real(fgm%tau(L))

  xmu0=x(2)
  fg0m%iw = one/(xi*wm + xmu0 - ed0 - gammam%iw -U*(n-0.5d0))
  call cfft_iw2tau(fg0m%iw,fg0m%tau,beta) 
  n0=-real(fg0m%tau(L)) 

  funcv(1)=nobj-n
  funcv(2)=nobj-n0
  write(*,"(4(f13.9))")n,n0,xmu,xmu0
end function funcv



