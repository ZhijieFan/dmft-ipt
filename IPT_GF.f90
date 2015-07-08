!###############################################################
!PURPOSE  : a library for Green's functions type and related 
!operations. This comes really handy when FFT must be used as it
!automatically set the correct dimensions.
!###############################################################
module IPT_GF
  USE SF_CONSTANTS, only: xi,pi
  USE SF_FFT_FFTPACK
  implicit none
  private 

  type,public                        :: matsubara_gf
     complex(8),dimension(:),pointer :: iw
     real(8),dimension(:),pointer    :: tau
     logical                         :: status=.false.
  end type matsubara_gf

  type,public                        ::  keldysh_component
     complex(8),dimension(:),pointer :: w
     complex(8),dimension(:),pointer :: t
  end type keldysh_component

  type,public                        ::  keldysh_gf
     type(keldysh_component)         :: less
     type(keldysh_component)         :: gtr
     type(keldysh_component)         :: ret
     logical                         :: status=.false.
  end type keldysh_gf

  interface allocate_gf
     module procedure &
          allocate_matsubara_gf     ,&
          allocate_keldysh_component,&
          allocate_keldysh_gf
  end interface allocate_gf

  interface deallocate_gf
     module procedure &
          deallocate_matsubara_gf     ,&
          deallocate_keldysh_component,&
          deallocate_keldysh_gf
  end interface deallocate_gf


  public :: allocate_gf
  public :: deallocate_gf
  public :: f_fft_rw2rt
  public :: f_fft_rt2rw

  interface assignment(=)
     module procedure &
          matsubara_gf_identity        ,&
          keldysh_component_gf_identity,&
          keldysh_gf_identity          ,&
                                !
          matsubara_gf_scalar_identity ,&
          keldysh_component_scalar_identity,&
          keldysh_gf_scalar_identity
  end interface assignment(=)
  public :: assignment(=)


  interface operator(+)
     module procedure &
          add_matsubara_gf      ,&
          add_keldysh_component ,&
          add_keldysh_gf
  end interface operator(+)
  public :: operator(+)

contains


  !+----------------------------------------------------------------+
  !PURPOSE  :  Allocate the matsubara and Keldysh GF used in the code.
  ! - note that the time component of the Keldysh GF are allocated 
  ! in the 0:L form. This is because we need an odd (suppose that L 
  ! is even) number of point to go from -T to T thru the zero. 
  ! This is required to perform products of the form f(t)*f(-t)
  ! in the expression of the self-energy and that other wise 
  ! would break the reversal symmetry.
  ! 
  ! The imag. time matsubara GF does not suffer of this problem
  ! because we work in the 0:beta only interval so we don't need 
  ! to go thru zero but include it in the interval 0:beta--> 1:L
  ! 
  ! Last but not the least the 0:L form of the Keldysh GF requires 
  ! to modify the real-time-->real-freq. FFT, which are added here 
  ! below.
  !+----------------------------------------------------------------+
  subroutine allocate_matsubara_gf(fgk,L)
    type(matsubara_gf),intent(inout) :: fgk
    integer,intent(in)               :: L
    allocate(fgk%iw(L),fgk%tau(L))
    fgk%status=.true.
  end subroutine allocate_matsubara_gf

  subroutine allocate_keldysh_component(fgk,L)
    type(keldysh_component),intent(inout) :: fgk
    integer,intent(in)                    :: L
    allocate(fgk%w(L),fgk%t(0:L))
  end subroutine allocate_keldysh_component

  subroutine allocate_keldysh_gf(fgk,L)
    type(keldysh_gf),intent(inout) :: fgk
    integer,intent(in)             :: L
    call allocate_keldysh_component(fgk%less,L)
    call allocate_keldysh_component(fgk%gtr,L)
    call allocate_keldysh_component(fgk%ret,L)
    fgk%status=.true.
  end subroutine allocate_keldysh_gf







  !+----------------------------------------------------------------+
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine deallocate_matsubara_gf(fgk)
    type(matsubara_gf),intent(inout)      :: fgk
    deallocate(fgk%iw,fgk%tau)
    fgk%status=.false.
  end subroutine deallocate_matsubara_gf

  subroutine deallocate_keldysh_component(fgk)
    type(keldysh_component),intent(inout) :: fgk
    deallocate(fgk%w,fgk%t)
  end subroutine deallocate_keldysh_component

  subroutine deallocate_keldysh_gf(fgk)
    type(keldysh_gf),intent(inout)        :: fgk
    call deallocate_keldysh_component(fgk%less)
    call deallocate_keldysh_component(fgk%gtr)
    call deallocate_keldysh_component(fgk%ret)
    fgk%status=.false.
  end subroutine deallocate_keldysh_gf






  !+----------------------------------------------------------------+
  !PURPOSE  : Modified real-axis FFT to deal with the special 0:L 
  ! form of the time-axis Keldysh GF.
  !+----------------------------------------------------------------+
  function f_fft_rw2rt(func_in) result(func_out)
    complex(8),dimension(:)               :: func_in
    complex(8),dimension(0:size(func_in)) :: func_out
    complex(8),dimension(size(func_in))   :: ftmp
    ftmp = func_in
    call fft(ftmp)
    call fftex(ftmp)
    func_out = fftshift(ftmp)*size(ftmp)
  end function f_fft_rw2rt
  !
  function f_fft_rt2rw(func_in) result(func_out)
    complex(8),dimension(:)               :: func_in
    complex(8),dimension(size(func_in)-1) :: func_out
    complex(8),dimension(size(func_in)-1) :: ftmp
    ftmp = func_in
    call ifft(ftmp)
    call fftex(ftmp)
    func_out = ifftshift(ftmp)
  end function f_fft_rt2rw





  !+----------------------------------------------------------------+
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine matsubara_gf_identity(fgk1,fgk2)
    type(matsubara_gf),intent(inout) :: fgk1
    type(matsubara_gf),intent(in)    :: fgk2
    fgk1%iw  = fgk2%iw
    fgk1%tau = fgk2%tau
  end subroutine matsubara_gf_identity

  subroutine keldysh_component_gf_identity(fgk1,fgk2)
    type(keldysh_component),intent(inout) :: fgk1
    type(keldysh_component),intent(in)    :: fgk2
    fgk1%t = fgk2%t
    fgk1%w = fgk2%w
  end subroutine keldysh_component_gf_identity

  subroutine keldysh_gf_identity(fgk1,fgk2)
    type(keldysh_gf),intent(inout) :: fgk1
    type(keldysh_gf),intent(in)    :: fgk2
    fgk1%less%t = fgk2%less%t
    fgk1%gtr%t  = fgk2%gtr%t
    fgk1%ret%t  = fgk2%ret%t
    fgk1%less%w = fgk2%less%w
    fgk1%gtr%w  = fgk2%gtr%w
    fgk1%ret%w  = fgk2%ret%w
  end subroutine keldysh_gf_identity







  !+----------------------------------------------------------------+
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine matsubara_gf_scalar_identity(fgk,C)
    type(matsubara_gf),intent(inout)           :: fgk
    complex(8),intent(in)                      :: C
    fgk%iw = C
    fgk%tau= dreal(C)
  end subroutine matsubara_gf_scalar_identity

  subroutine keldysh_component_scalar_identity(fgk,C)
    type(keldysh_component),intent(inout)      :: fgk
    complex(8),intent(in)                      :: C
    fgk%t = C
    fgk%w = C
  end subroutine keldysh_component_scalar_identity

  subroutine keldysh_gf_scalar_identity(fgk,C)
    type(keldysh_gf),intent(inout) :: fgk
    complex(8),intent(in)          :: C
    fgk%less%t = C
    fgk%gtr%t  = C
    fgk%ret%t  = C
    fgk%less%w = C
    fgk%gtr%w  = C
    fgk%ret%w  = C
  end subroutine keldysh_gf_scalar_identity







  !+----------------------------------------------------------------+
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  elemental function add_matsubara_gf(fgk1,fgk2) result(fgk3)
    type(matsubara_gf),intent(in) :: fgk1,fgk2
    type(matsubara_gf)            :: fgk3
    fgk3%iw  = fgk1%iw  + fgk2%iw
    fgk3%tau = fgk1%tau + fgk2%tau
  end function add_matsubara_gf

  elemental function add_keldysh_component(fgk1,fgk2) result(fgk3)
    type(keldysh_component),intent(in) :: fgk1,fgk2
    type(keldysh_component)            :: fgk3
    fgk3%t = fgk1%t + fgk2%t
    fgk3%w = fgk1%w  + fgk2%w
  end function add_keldysh_component

  elemental function add_keldysh_gf(fgk1,fgk2) result(fgk3)
    type(keldysh_gf),intent(in) :: fgk1,fgk2
    type(keldysh_gf)            :: fgk3
    fgk3%less%t = fgk1%less%t + fgk2%less%t
    fgk3%gtr%t  = fgk1%gtr%t  + fgk2%gtr%t
    fgk3%ret%t  = fgk1%ret%t  + fgk2%ret%t
    fgk3%less%w = fgk1%less%w + fgk2%less%w
    fgk3%gtr%w  = fgk1%gtr%w  + fgk2%gtr%w
    fgk3%ret%w  = fgk1%ret%w  + fgk2%ret%w
  end function add_keldysh_gf




  ! !+----------------------------------------------------------------+
  ! !PURPOSE  : 
  ! !+----------------------------------------------------------------+
  ! function ret_component_t(fgkgtr,fgkless,t) result(fgkret)
  !   integer                                       :: i,L
  !   complex(8),dimension(:),intent(in)            :: fgkgtr
  !   complex(8),dimension(size(fgkgtr)),intent(in) :: fgkless
  !   complex(8),dimension(size(fgkgtr))            :: fgkret
  !   real(8),dimension(size(fgkgtr)),intent(in)   :: t
  !   L=size(fgkgtr)
  !   forall(i=1:L)fgkret(i)=step(t(i))*(fgkgtr(i)-fgkless(i))
  ! contains
  !   pure function step(x)
  !     real(8),intent(in) :: x
  !     real(8)            :: step
  !     if(x < 0.d0) then
  !        step = 0.0d0
  !     elseif(x==0.d0)then
  !        step = 0.50d0
  !     else
  !        step = 1.0d0
  !     endif
  !   end function step
  ! end function ret_component_t



  ! !+----------------------------------------------------------------+
  ! !PURPOSE  : 
  ! !+----------------------------------------------------------------+
  ! function less_component_w(fret,wr,beta) result(fless)
  !   integer                            :: i,L
  !   complex(8),dimension(:),intent(in) :: fret
  !   complex(8),dimension(size(fret))   :: fless
  !   real(8),dimension(size(fret))      :: wr
  !   real(8)                            :: A,beta,w
  !   L=size(fret)
  !   do i=1,L
  !      w       = wr(i)
  !      A       = -dimag(fret(i))/pi
  !      fless(i)= 2d0*pi*xi*fermi(w,beta)*A
  !   enddo
  ! contains
  !   function fermi(x,beta)
  !     real(8) :: fermi, x, beta
  !     if(x*beta > 50d0)then
  !        fermi=0.d0
  !        return
  !     endif
  !     fermi = 1.d0/(1.d0+exp(beta*x))
  !   end function fermi
  ! end function less_component_w



  ! !+----------------------------------------------------------------+
  ! !PURPOSE  : 
  ! !+----------------------------------------------------------------+
  ! function gtr_component_w(fret,wr,beta) result(fgtr)
  !   integer                            :: i,L
  !   complex(8),dimension(:),intent(in) :: fret
  !   complex(8),dimension(size(fret))   :: fgtr
  !   real(8),dimension(size(fret))      :: wr
  !   real(8)                            :: A,beta,w
  !   L=size(fret)
  !   do i=1,L
  !      w      = wr(i)
  !      A      = -dimag(fret(i))/pi
  !      fgtr(i)= 2.d0*pi*xi*(fermi(w,beta)-1.d0)*A
  !   enddo
  ! contains
  !   function fermi(x,beta)
  !     real(8) :: fermi, x, beta
  !     if(x*beta > 50d0)then
  !        fermi=0.d0
  !        return
  !     endif
  !     fermi = 1.d0/(1.d0+exp(beta*x))
  !   end function fermi
  ! end function gtr_component_w

end module IPT_GF
