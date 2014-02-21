!###############################################################
!     PURPOSE  : Contains global variables
!     AUTHORS  : Adriano Amaricci
!###############################################################
! DESCRIPTION
!   provide access to a set of solvers for the dmft based on IPT
! OPTIONS
!   u=[2]      -- interaction
!   beta=[100] -- inverse temperature
!   xmu=[0]    -- chemical potential
!   nloop=[10] -- max number of iterations
!   L=[2048]   -- number of frequencies
!   ts=[0.5]   -- n.n. hopping parameter
!   tsp=[0]    -- n.n.n. hopping parameter
!   nx=[20]    -- number of points in energy/k-grid
!   wmax=[5]   -- max frequency on real axis
!   eps=[0.01] -- broadening parameter
!   deltasc=[0.1]     -- breaking symmetry parameter
!   nread=[0.0]       -- required density look for mu
!   nerror=[1.d-4]    -- error treshold for mu-loop
!   ndelta=[0.1]      -- mu-step in mu-loop for fixed density
!   eps_error=[1.D-4] -- error treshold
!   success=[2]       -- number of converged events
! ROUTINES
!   normal:
!   solve_ipt_keldysh(fg0_,wr_,t_) result(sigma_) dim(-L:L) 
!   solve_ipt_matsubara(fg0_)      result(sigma_) dim(L)
!   solve_ipt_sopt(fg0_,wr_)       result(sigma_) dim(-L:L)
!   solve_ipt_zeroT(fg0_,dw_,dt_)  result(sigma_) dim(-L:L)
!   solve_mpt_sopt(fg0_,wr_,n_,n0_,xmu0_)  result(sigma_) dim(-L:L)
!   solve_mpt_matsubara(fg0_,n_,n0_,xmu0_) result(sigma_) dim(L)
!  
!   supercond:
!   solve_ipt_sc_sopt(fg0_,wr_,delta_)  result(sigma_) dim(2,-L:L)
!   solve_ipt_sc_matsubara(fg0_,delta_) result(sigma_) dim(2,L)
!   solve_mpt_sc_sopt(fg0_,wr_,n_,n0_,delta_,delta0_)  result(sigma_) dim(2,-L:L)
!   solve_mpt_sc_matsubara(fg0_,n_,n0_,delta_,delta0_) result(sigma_) dim(2,L)
module IPT_VARS_GLOBAL
  USE COMMON_VARS
  USE PARSE_CMD
  USE GREENFUNX
  USE FFTGF
  USE INTEGRATE
  USE TOOLS
  USE ARRAYS
  USE FUNCTIONS
  implicit none

  !Size of the problem: some parameters:
  !=========================================================
  integer :: L              !a large number fix the size of the problem
  integer :: nloop          !dmft loop variables
  real(8) :: d              !bandwidth
  real(8) :: ts             !hopping amplitude
  real(8) :: u              !local  interaction
  real(8) :: xmu            !chemical potential
  real(8) :: dt             !time step
  real(8) :: fmesh          !freq. step
  real(8) :: beta           !inverse temperature
  real(8) :: eps            !broadening
  integer :: Nx
  real(8) :: wmin,wmax
  real(8) :: dmft_error
  integer :: Nsuccess
  real(8) :: weight
  real(8) :: deltasc
  !TO BE REMOVED
  real(8) :: nread,nerror,ndelta

  !Namelists:
  !=========================================================
  namelist/variables/L,beta,U,ts,xmu,wmax,Nx,nloop,eps,weight,&
       dmft_error,Nsuccess,deltasc,nread,nerror,ndelta

contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input(inputFILE)
    character(len=*) :: inputFILE
    integer :: i
    logical :: control
    !variables: default values
    U     = 2.d0
    beta  = 100.d0
    ts    = 0.5d0
    xmu   = 0.d0
    Nx    = 20
    nloop = 10
    eps   = 0.01d0
    wmax  = 5.d0
    L     = 2048
    weight= 0.9d0
    dmft_error= 1.d-4
    Nsuccess = 2
    deltasc  = 0.1d0
    nread=0.d0
    nerror=1.d-4
    ndelta=0.1d0

    inquire(file=adjustl(trim(inputFILE)),exist=control)
    if(control)then
       open(10,file=adjustl(trim(inputFILE)))
       read(10,nml=variables)
       close(10)
    else
       open(10,file="default."//adjustl(trim(inputFILE)))
       write(10,nml=variables)
       close(10)
       call error("can not open INPUT file, dumping a default version in default."//adjustl(trim(inputFILE)))
    endif

    call parse_cmd_variable(u,"U")
    call parse_cmd_variable(beta,"BETA")
    call parse_cmd_variable(ts,"TS")
    call parse_cmd_variable(xmu,"XMU")
    call parse_cmd_variable(nx,"NX")
    call parse_cmd_variable(nloop,"NLOOP")
    call parse_cmd_variable(L,"L")
    call parse_cmd_variable(eps,"EPS")
    call parse_cmd_variable(wmax,"WMAX")
    call parse_cmd_variable(weight,"WEIGHT")
    call parse_cmd_variable(dmft_error,"DMFT_ERROR")
    call parse_cmd_variable(nsuccess,"NSUCCESS")
    call parse_cmd_variable(deltasc,"DELTASC")
    call parse_cmd_variable(nread,"NREAD")
    call parse_cmd_variable(nerror,"NERROR")
    call parse_cmd_variable(ndelta,"NDELTA")

    write(*,nml=variables)
    open(10,file="used."//adjustl(trim(inputFILE)))
    write(10,nml=variables)
    close(10)
    return
  end subroutine read_input


  !******************************************************************
  !******************************************************************
  !******************************************************************


  function get_local_density(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_local_density


end module IPT_VARS_GLOBAL
