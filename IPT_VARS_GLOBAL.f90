module IPT_VARS_GLOBAL
  USE SF_CONSTANTS
  USE SF_VERSION
  USE SF_PARSE_INPUT

  implicit none

  !Size of the problem: some parameters:
  !=========================================================
  !integer :: Lmats          !a large number fix the size of the problem
  !this is now read from the driver, because the dimensions are passed as input in the solvers
  integer :: Norb           !number of orbitals
  integer :: Nspin          !number of spin channels
  integer :: nloop          !dmft loop variables
  real(8) :: uloc(3)        !local  interaction
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

    call parse_input_variable(uloc,"ULOC",inputFILE,default=[2d0,0d0,0d0],comment="Local interaction term")
    call parse_input_variable(beta,"BETA",inputFILE,default=100.d0,comment="inverse temperature")
    call parse_input_variable(xmu,"XMU",inputFILE,default=0.d0,comment="chemical potential w/ respect to HF shift")
    call parse_input_variable(Norb,"Norb",inputFILE,default=1,comment="number of orbitals")
    call parse_input_variable(Nspin,"Nspin",inputFILE,default=1,comment="number of spin channels")
    call parse_input_variable(nloop,"NLOOP",inputFILE,default=100,comment="maximum number of DMFT iterations")
    call parse_input_variable(eps,"EPS",inputFILE,default=0.01d0,comment="real-axis broadening")
    call parse_input_variable(wmax,"WMAX",inputFILE,default=5.d0,comment="max frequency on the real-axis")
    call parse_input_variable(dmft_error,"DMFT_ERROR",inputFILE,default=1.d-5,comment="DMFT error thereshold")
    call parse_input_variable(nsuccess,"NSUCCESS",inputFILE,default=1,comment="DMFT max number of consecutive successes")
    call parse_input_variable(deltasc,"DELTASC",inputFILE,default=0.1d0,comment="superconducting channel breaking symmetry parameter")

    call save_input_file(INPUTFILE)
    !
    call scifor_version()
    call code_version(revision)

  end subroutine read_input


end module IPT_VARS_GLOBAL
