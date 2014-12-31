module IPT_VARS_GLOBAL
  USE CONSTANTS
  USE PARSE_INPUT
  USE SCIFOR_VERSION
  implicit none

  !Size of the problem: some parameters:
  !=========================================================
  integer :: L              !a large number fix the size of the problem
  integer :: nloop          !dmft loop variables
  real(8) :: uloc           !local  interaction
  real(8) :: xmu            !chemical potential
  real(8) :: beta           !inverse temperature
  real(8) :: eps            !broadening
  real(8) :: wmax
  real(8) :: dmft_error
  integer :: Nsuccess
  real(8) :: deltasc

  include "revision.inc"


contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input(inputFILE)
    character(len=*) :: inputFILE

    call parse_input_variable(uloc,"ULOC",inputFILE,default=2.d0,comment="Local interaction term")
    call parse_input_variable(beta,"BETA",inputFILE,default=100.d0)
    call parse_input_variable(xmu,"XMU",inputFILE,default=0.d0)
    call parse_input_variable(nloop,"NLOOP",inputFILE,default=100)
    call parse_input_variable(L,"L",inputFILE,default=2048,comment="number of frequencies")
    call parse_input_variable(eps,"EPS",inputFILE,default=0.01d0)
    call parse_input_variable(wmax,"WMAX",inputFILE,default=5.d0)
    call parse_input_variable(dmft_error,"DMFT_ERROR",inputFILE,default=1.d-5)
    call parse_input_variable(nsuccess,"NSUCCESS",inputFILE,default=1)
    call parse_input_variable(deltasc,"DELTASC",inputFILE,default=0.1d0)
    !call parse_input_variable(wmix,"WMIX",inputFILE,default=0.75d0)

    call save_input_file(INPUTFILE)
    call sf_version(revision)

  end subroutine read_input


end module IPT_VARS_GLOBAL
