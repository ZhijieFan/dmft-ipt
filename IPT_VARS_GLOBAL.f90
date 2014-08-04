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
  USE SCIFOR_VERSION
  USE CONSTANTS
  USE PARSE_INPUT
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
  real(8) :: ust            !local  interaction
  real(8) :: xmu            !chemical potential
  real(8) :: dt             !time step
  real(8) :: fmesh          !freq. step
  real(8) :: beta           !inverse temperature
  real(8) :: eps            !broadening
  real(8) :: wmin,wmax
  real(8) :: dmft_error
  integer :: Nsuccess
  real(8) :: weight
  real(8) :: deltasc
  !TO BE REMOVED
  !real(8) :: nread,nerror,ndelta

  include "revision.inc"


contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input(inputFILE)
    character(len=*) :: inputFILE

    call parse_input_variable(u,"U",inputFILE,default=2.d0,comment="")
    call parse_input_variable(ust,"UST",inputFILE,default=0.d0)
    call parse_input_variable(beta,"BETA",inputFILE,default=100.d0)
    call parse_input_variable(ts,"TS",inputFILE,default=0.5d0)
    call parse_input_variable(xmu,"XMU",inputFILE,default=0.d0)
    call parse_input_variable(nloop,"NLOOP",inputFILE,default=100)
    call parse_input_variable(L,"L",inputFILE,default=2048)
    call parse_input_variable(eps,"EPS",inputFILE,default=0.01d0)
    call parse_input_variable(wmax,"WMAX",inputFILE,default=5.d0)
    call parse_input_variable(dmft_error,"DMFT_ERROR",inputFILE,default=1.d-4)
    call parse_input_variable(nsuccess,"NSUCCESS",inputFILE,default=2)
    call parse_input_variable(deltasc,"DELTASC",inputFILE,default=0.1d0)
    ! call parse_input_variable(nread,"NREAD",inputFILE,default=0.d0)
    ! call parse_input_variable(nerror,"NERROR",inputFILE,default=1.d-4)
    ! call parse_input_variable(ndelta,"NDELTA",inputFILE,default=0.1d0)

    call save_input_file(INPUTFILE)
    call version(revision)

  end subroutine read_input


end module IPT_VARS_GLOBAL
