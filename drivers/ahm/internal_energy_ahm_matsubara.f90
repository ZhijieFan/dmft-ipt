subroutine get_sc_internal_energy
  integer    :: i,ik
  real(8)    :: matssum,fmatssum,checkP,checkdens,vertex,Dssum
  complex(8) :: iw,gkw,fkw,g0kw,f0kw
  real(8)    :: Epot,Etot,Eint,kin,kinsim,Ds,docc
  real(8)    :: Sigma_infty,S_infty,det,det_infty,csi,Ei,thermal_factor
  real(8)    :: free(Lk),Ffree(Lk),n_k(Lk)

  !Get asymptotic self-energies
  Sigma_infty =   real(sigma(1,L),8)
  S_infty     =   real(sigma(2,L),8)

  checkP=0.d0 ; checkdens=0.d0 ;          ! test variables

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
  enddo

  !  call splot("nk0_distribution.last",epsik,2.d0*free)
  !  call splot("fnk0_distribution.last",epsik,Ffree)

  kinsim=0
  kinsim = sum(fg(1,:)*fg(1,:)+conjg(fg(1,:))*conjg(fg(1,:))-2.d0*fg(2,:)*fg(2,:))*2.d0*ts**2/beta
  ! kinsim = kinsim - 1.d0/(2*pi*wm(L)) !check out where this term comes from??!!

  Epot=zero
  Epot = sum(fg(1,:)*sigma(1,:) + fg(2,:)*sigma(2,:))/beta*2.d0

  docc = 0.5d0*n**2
  if(u > 0.01d0)docc=-Epot/u + n - 0.25d0


  Eint=kin+Epot

  Ds=zero
  Ds = sum(fg(2,:)*fg(2,:))/beta*2.d0

  ! kinsim=zero
  ! do i=1,L
  !    kinsim=kinsim + fg(1,i)*fg(1,i) + conjg(fg(1,i))*conjg(fg(1,i)) - 2.d0*fg(2,i)*fg(2,i)
  ! enddo
  ! kinsim=0.5d0*kinsim/Beta
  ! kinsim=kinsim- 1.d0/(2*pi*om(Iwmax2))

  write(*,*)'========================================='      
  write(*,*)"Asymptotic Self-Energies",Sigma_infty, S_infty
  write(*,*)'========================================='
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
  call splot("nk_distribution.ipt",epsik,n_k,2.d0*free)
  call splot("thermodynamics.ipt",L,n,0.5d0*checkdens,kin,kinsim,docc,Ds,append=TT)

  !  call splot("fnk0_distribution.last",epsik,Ffree)
  !   write(42,143)xmu,totdens,checkdens,dreal(kinsim),Tcheck,double
  ! 143 format(f6.4,5(1x,f12.9))
  return 
end subroutine get_sc_internal_energy
