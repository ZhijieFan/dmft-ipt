program hmipt_real
  USE DMFT_IPT
  USE SCIFOR
  USE DMFT_TOOLS

  implicit none
  logical                                 :: converged,check
  real(8)                                 :: wmix,D
  integer                                 :: i,ik,iloop,L,Lk,Nx,unit
  complex(8)                              :: zeta
  complex(8),allocatable                  :: fg(:),fg0(:),sigma(:),fg0_prev(:)
  real(8),allocatable                     :: wreal(:)
  real(8)                                 :: n,docc,z,wmin
  character(len=24)                       :: finput
  real(8),dimension(:),allocatable        :: Wtk
  complex(8),dimension(:,:,:),allocatable :: Hk

  call parse_cmd_variable(finput,"finput",default="inputIPT.conf")
  call parse_input_variable(D,"D",finput,default=1d0)
  call parse_input_variable(L,"L",finput,default=2000)
  call parse_input_variable(wmix,"WMIX",finput,default=0.5d0)
  call read_input(finput)


  !allocate functions:
  allocate(fg(L))
  allocate(sigma(L))
  allocate(fg0(L))
  allocate(fg0_prev(L))

  !build freq. array
  allocate(wreal(L))
  wreal = linspace(-1d0*wmax,wmax,L)

  !get or read first sigma 
  inquire(file="Sigma.restart",exist=check)
  sigma=zero
  if(check)then
     call read_array("Sigma.restart",sigma)
  endif

  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5,2x)",advance="no")"DMFT-loop",iloop
     !
     !SELF-CONSISTENCY:
     do i=1,L
        zeta = dcmplx(wreal(i),eps) - sigma(i)
        fg(i) = gfbether(wreal(i),zeta,D)
     enddo
     !
     fg0_prev=fg0
     fg0 = one/(one/fg + sigma)
     ! if(iloop>1)fg0 = wmix*fg0 + (1.d0-wmix)*fg0_prev
     call broyden_mix(fg0_prev,fg0,wmix,5,iloop)
     !
     !IMPURITY SOLVER
     sigma = ipt_solve_real(fg0,wreal)
     !sigma = ph_symmetry(sigma)
     !
     converged=check_convergence(fg0,dmft_error,nsuccess,nloop)
     !
  enddo


  call splot("Gloc_wreal.ipt",wreal,fg)
  call splot("G0_wreal.ipt",wreal,fg0)
  call splot("Sigma_wreal.ipt",wreal,sigma)
  call save_array("Sigma.restart",sigma)


contains

  function ph_symmetry(func) result(sfunc)
    complex(8),dimension(:)          :: func
    complex(8),dimension(size(func)) :: sfunc
    complex(8),dimension(size(func)) :: ifunc
    !Re(F(-w)) = -Re(F(w))
    !Im(F(w))  = Im(F(w))
    L = size(func)
    ifunc =  func(L:1:-1)
    sfunc = (dreal(func) - dreal(ifunc))/2d0 + xi*(dimag(func) + dimag(ifunc))/2d0
  end function ph_symmetry
end program hmipt_real
