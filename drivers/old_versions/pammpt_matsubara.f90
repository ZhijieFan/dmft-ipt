!########################################################
!     Program  : PAMMPT
!     TYPE     : Main program
!     PURPOSE  : Solve the PAM using DMFT away from half-filling
! with modified itertive perturbation scheme (MPT)  
!     AUTHORS  : Adriano Amaricci
!########################################################
module COMMON
  USE BROYDEN
  USE DMFT_IPT
  implicit none
  !Put here vars in common with the BROYDN function
  real(8)                :: xmu0,n,n0,zd,zp
  type(matsubara_gf)     :: fg,fg0,fgp
  complex(8),allocatable :: sigma(:),sigmap(:)
  complex(8),allocatable :: gamma(:)
  real(8),allocatable    :: wm(:),tau(:)
end module COMMON

function funcv(x)
  USE COMMON
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0=x(1)
  fg0%iw = one/(xi*wm + xmu0 - ed0 - gamma - u*(n-0.5d0))
  call fftgf_iw2tau(fg0%iw,fg0%tau,beta)
  n0=-real(fg0%tau(L))
  funcv(1)=n-n0
  write(*,"(3(f13.9))")n,n0,xmu0
end function funcv

program pammpt
  USE COMMON
  USE DMFT_IPT
  USE IOTOOLS
  implicit none
  integer                :: i,iloop
  real(8)                :: x(1)
  logical                :: check,converged
  complex(8)             :: zeta,alpha
  real(8)                :: np,gzero,gmu,ntot,dtau
  complex(8),allocatable :: sold(:)

  call read_input("inputIPT.in")

  call allocate_gf(fg,L)
  call allocate_gf(fg0,L)
  call allocate_gf(fgp,L)
  allocate(sigma(L),sigmap(L),sold(L))
  allocate(gamma(L))
  allocate(wm(L),tau(0:L))
  wm  = pi/beta*real(2*arange(1,L)-1,8)
  tau = linspace(0.d0,beta,L+1,mesh=dtau)
  gmu=xmu  ; gzero=0.d0
  if((ed0-ep0) > 0.d0)gzero=0.5*(ep0+ed0+sqrt((ep0-ed0)**2 + 4*Vpd**2))
  if((ed0-ep0) < 0.d0)gzero=0.5*(ep0+ed0-sqrt((ep0-ed0)**2 + 4*Vpd**2))
  if((ed0-ep0) /=0.d0)xmu=gmu+gzero !true ED chemical potential
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  D=2.d0*ts
  n=0.5d0
  sigma=zero
  xmu0=xmu
  iloop=0
  converged=.false. ; sold=sigma
  do while (.not.converged)
     iloop=iloop+1
     write(*,"(A,i5)")"DMFT-loop",iloop

     do i=1,L
        alpha     = xi*wm(i)  +xmu -ed0 - sigma(i)
        sigmap(i) = Vpd**2/alpha
        zeta      = xi*wm(i)  +xmu -ep0 - sigmap(i)
        fgp%iw(i) = gfbethe(wm(i),zeta,D)
        fg%iw(i)  = one/alpha + Vpd**2/alpha**2*fgp%iw(i)
     enddo
     call fftgf_iw2tau(fg%iw,fg%tau,beta)
     call fftgf_iw2tau(fgp%iw,fgp%tau,beta)
     n   =-real(fg%tau(L)) ; np  =-2.d0*real(fgp%tau(L))
     ntot= np+2.d0*n
     print*,2.d0*n,np,ntot


     !Get the hybridization functions: \Gamma(w+xmu)
     gamma= xi*wm + xmu - ed0 - sigma - one/fg%iw

     !Fix the xmu0 w/ condition n0=n
     x(1)=xmu
     call broydn(x,check)
     xmu0=x(1)

     sigma= solve_mpt_matsubara(fg0%iw,n,n0,xmu0)
     sigma=weight*sigma + (1.d0-weight)*sold
     sold=sigma
     converged=check_convergence(sigma,eps_error,Nsuccess,Nloop)
     zd=1.d0/(1.d0+abs(dimag(Sigma(1))/wm(1)))
     zp=1.d0/(1.d0+abs(dimag(Sigmap(1))/wm(1)))
     call splot("observables.ipt",xmu,u,beta,dble(iloop),n,np,ntot,zd,zp,append=TT)
  enddo
  call splot("observables_last.ipt",xmu,u,vpd,beta,n,np,ntot,zd,zp,append=printf)
  call splot("nd.np.ntot.ipt",n,np,ntot,append=TT)
  call splot("Sigma_iw.ipt",wm,sigma,append=printf)
  call splot("Sigmap_iw.ipt",wm,sigmap,append=printf)
  call splot("G_iw.ipt",wm,fg%iw,append=printf)
  call splot("Gp_iw.ipt",wm,fgp%iw,append=printf)
  call splot("G_tau.ipt",tau,fg%tau,append=printf)
  call splot("Gp_tau.ipt",tau,fgp%tau,append=printf)
  !call pam_getener()


  ! contains

  !   !+------------------------------------------------------------------+
  !   !PROGRAM  : 
  !   !TYPE     : subroutine
  !   !PURPOSE  : 
  !   !+------------------------------------------------------------------+
  !   subroutine pam_getener()
  !     implicit none
  !     integer :: i,Nf
  !     real(8) :: depsi,e,free
  !     real(8) :: w,doble,temp
  !     real(8) :: znd,znp,zntot,znpd
  !     real(8) :: Etot,Ekin,Epot,Ehyb,Eepsi,Emu
  !     real(8) :: kag,kag0,coda1
  !     complex(8) :: iw,greenp0,greenp,selfp,greend0,greend,selfd,zita,zita0
  !     complex(8) :: ptail,dtail,gamma

  !     write(*,"(A)")"Getting the energy"
  !     temp=1.d0/beta
  !     free=0.d0
  !     depsi=2.d0*d/dble(L)
  !     do i=1,L
  !        e=-d+dble(i)*depsi
  !        free=free+depsi*e*dens_bethe(e)*fermi(e,beta)
  !     enddo

  !     Nf=L
  !     write(*,"(A,I10)")"Get En with=",Nf

  !     Ekin=zero
  !     Ehyb=zero
  !     Epot=0.d0
  !     znd=0.d0
  !     znp=0.d0
  !     do i=1,Nf
  !        w=pi/beta*dble(2*i-1) ; iw=xi*w
  !        selfp=sigmap(i)
  !        selfd=sigma(i)
  !        zita0  = iw
  !        zita   = iw+xmu-ep0-selfp
  !        greenp0= gfbethe(w,zita0,d)
  !        greenp = gfbethe(w,zita,d)
  !        gamma  = iw+xmu-ed0-selfd
  !        greend = one/gamma + Vpd**2/gamma**2*greenp

  !        !Kinetic energy     
  !        Ekin=Ekin+2.d0*(zita*greenp-zita0*greenp0)
  !        !Hybridization energy
  !        Ehyb=Ehyb+4.d0*(selfp*greenp)
  !        ! !Potential Energy second type
  !        Epot=Epot+(selfd*greend)

  !        znd=znd+2.d0*temp*real(greend)
  !        znp=znp+2.d0*temp*real(greenp)
  !     enddo

  !     znd=2.d0*znd;znp=2.d0*znp
  !     zntot=znd+znp

  !     Ekin=2.d0*Ekin
  !     Ehyb=2.d0*Ehyb
  !     Epot=2.d0*Epot

  !     Ekin=temp*Ekin+2.d0*free
  !     Epot=temp*Epot
  !     Ehyb=temp*Ehyb 
  !     Eepsi=ed0*znd + ep0*znp
  !     Emu=-xmu*zntot!-2.d0*temp*xmu*zntot

  !     Etot=Ekin+Epot+Ehyb+Eepsi+Emu

  !     doble=0.5d0*(znd+1.d0) - 0.25d0
  !     if(u > 0.01d0)doble=real(Epot)/u + 0.5d0*(znd+1.d0) - 0.25d0

  !     print*,""
  !     write(*,"(A,f13.9)")"docc   =",doble
  !     write(*,"(A,f13.9,f13.9)")"nd     =",znd+1.d0,2.d0*n
  !     write(*,"(A,f13.9,f13.9)")"np     =",znp+1.d0,np
  !     write(*,"(A,f13.9,f13.9)")"ntot   =",znp+znd+2.d0,ntot

  !     open(20,file="Ekin.pot.hyb.tot.ipt",access="append")
  !     open(21,file="Efree.mu.epsi.ipt",access="append")
  !     open(22,file="End.np.ntot.docc.ipt",access="append")
  !     write(20,"(4(f13.9))")Ekin,Epot,Ehyb,Etot
  !     write(21,"(3(f13.9))")free,Emu,Eepsi
  !     write(21,"(4(f13.9))")znd+1.d0,znp+1.d0,zntot+2.d0,doble
  !     close(20);close(21);close(22)
  !     return
  !   end subroutine pam_getener
  !   !*********************************************************************
  !   !*********************************************************************
  !   !*********************************************************************

end program pammpt







