  !########################################################
  !     Program  : AHMIPT
  !     TYPE     : Main program
  !     PURPOSE  : Solve the Attractive Hubbard Model using IPT
  !     AUTHORS  : Adriano Amaricci
  !########################################################
  include "IPT_VARS_GLOBAL"
  include "IPT_SC_FUNX_SOPT"
  program hmipt
    !USE DMFT_IPT
    USE IPT_VARS_GLOBAL
    USE IPT_SC_FUNX_SOPT
    implicit none
    integer    :: ik,Lk
    logical    :: converged
    complex(8) :: det,zeta1,zeta2
    real(8)    :: delta,n
    !
    complex(8),allocatable          :: sigma(:,:),fg(:,:),wf0(:,:),calG(:,:)
    complex(8),allocatable          :: zeta(:)
    type(matsubara_gf),dimension(2) :: gf,sf
    !
    real(8),allocatable :: wt(:),epsik(:)

    call read_input("inputIPT.in")
    allocate(fg(2,-L:L),sigma(2,-L:L))
    allocate(wf0(2,-L:L),calG(2,-L:L))
    allocate(zeta(-L:L))

    Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
    call bethe_lattice(wt,epsik,Lk,D_=D,eps_=eps)

    delta=deltasc 
    sigma=zero ; sigma(2,:)=-delta 
    iloop=0    ; converged=.false.
    do while (.not.converged)
       iloop=iloop+1
       write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
       !
       ! G(iw) F(iw)      --> G(iw) F(iw)
       ! F*(-iw) -G(-iw)  --> F(iw) -G*(iw)
       !
       !analytic continuation: iw--> w+ie
       ! G(w+ie) F(w+ie)     --> G(w) F(w)
       ! F*(-w-ie) -G(-w-ie) --> F(-w) -G*(-w)
       !
       ! w+ie-Sigma(w)-e(k); -S(w)        --> zeta(w)-e(k);-S(w)
       ! -S(-w); -(-w-ie-Sigma*(-w)-e(k)) --> -S(-w); -(zeta*(-w)-e(k))
       !
       ! G = (zeta*(-w)-e(k))/{(zeta(w)-e(k))(zeta*(-w)-e(k))+S(w)S(-w)}
       ! F = -S(-w)/{(zeta(w)-e(k))(zeta*(-w)+e(k))+S(w)S(-w)}
       fg=zero
       zeta(:) = cmplx(wr(:),eps,8) + xmu - sigma(1,:)
       do i=-L,L
          zeta1 = zeta(i)
          zeta2 = conjg(zeta(-i))!cmplx(wr(i),eps,8) + sigma(1,-i)
          do ik=1,Lk
             det = (zeta1-epsik(ik))*(zeta2-epsik(ik)) + sigma(2,i)*sigma(2,i)
             fg(1,i)=fg(1,i) + wt(ik)*(zeta2-epsik(ik))/det
             fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/det
          enddo
       enddo
       delta=-u*sum(aimag(fg(2,-L:L))*fermi(wr(-L:L),beta))*fmesh/pi
       n=-sum(aimag(fg(1,-L:L))*fermi(wr(-L:L),beta))*fmesh/pi


       !get the Weiss Field \calG0^-1(w)
       !
       ! wf(1)=calG0^-1 = G*(-w)/(G(w)G*(-w) + F(w)F(-w)) + Sigma(w)
       ! wf(2)=calF0^-1 =-F(-w) /(G(w)G*(-w) + F(w)F(-w)) + S(w)
       do i=-L,L
          det      = fg(1,i)*conjg(fg(1,-i)) + fg(2,i)*fg(2,i)
          wf0(1,i) = conjg(fg(1,-i))/det + sigma(1,i)
          wf0(2,i) = fg(2,i)/det         + sigma(2,i) + delta
       end do
       do i=-L,L
          det       =  wf0(1,i)*conjg(wf0(1,-i)) + wf0(2,i)*wf0(2,-i)
          calG(1,i)=  conjg(wf0(1,-i))/det
          calG(2,i)=  wf0(2,-i)/det
       end do

       write(*,"(2f14.9)",advance="no")2.d0*n,delta
       sigma =  solve_ipt_sc_sopt(calG(1:2,:),fg,delta)
       converged = check_convergence(sigma(1,:)+sigma(2,:))
    enddo

    call splot("DOS.last",wr,-aimag(fg(1,:))/pi,append=TT)
    call splot("G_realw.last",wr,fg(1,:),append=TT)
    call splot("F_realw.last",wr,fg(2,:),append=TT)
    call splot("Sigma_realw.last",wr,sigma(1,:),append=TT)
    call splot("Self_realw.last",wr,sigma(2,:),append=TT)
    call splot("calG0_realw.last",wr,calG(1,:),append=TT)
    call splot("calF0_realw.last",wr,calG(2,:),append=TT)
    call splot("n.delta_realw.last",n,delta,append=TT)

    deallocate(wm,tau)
    call init_wmgrid(wm,beta,4*L)
    call init_taugrid(tau,beta/dble(4*L),4*L)
    call allocate_gf(gf,4*L);call allocate_gf(sf,4*L)
    call getGmats(wr,fg(1,:),gf(1)%iw,beta)
    call getGmats(wr,fg(2,:),gf(2)%iw,beta)
    call fftgf_iw2tau(gf(1)%iw,gf(1)%tau,beta)
    call fft_iw2tau(gf(2)%iw,gf(2)%tau,beta,4*L)
    print*,-2.d0*gf(1)%tau(4*L),-u*gf(2)%tau(0)
    call splot('GM_iw.last',wm,gf(1)%iw,append=TT)
    call splot('FM_iw.last',wm,gf(2)%iw,append=TT)
    call splot('GM_tau.last',tau,gf(1)%tau,append=TT)
    call splot('FM_tau.last',tau,gf(2)%tau,append=TT)
  end program hmipt
