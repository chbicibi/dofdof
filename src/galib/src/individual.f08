module individual
  use interface
  use util
  implicit none

  private
  public :: TIndiv, TIndivI, TIndivS, TIndivSI, TIndivC, TIndivCI

  type :: TIndiv
    integer :: id = -1
    integer, allocatable :: parents_id(:)

    real(8), allocatable :: dvariables(:)
    integer, allocatable :: ivariables(:)
    real(8), allocatable :: objectives(:)
    real(8), allocatable :: constraints(:)
    logical :: feasible = .true.
    logical :: evaluated = .false.

    contains

    generic :: initialize => initialize_n
    generic :: set_variables => set_ivariables, set_dvariables
    generic :: get_variables => get_ivariables, get_dvariables, get_dvariables_proc
    generic :: obj => obj_1, obj_n
    generic :: print => print_u, print_up
    generic :: print_wv => print_wv_u, print_wv_up
    generic :: print_header => print_header_u
    generic :: print_header_wv => print_header_wv_u

    procedure :: set_id
    procedure :: initialize_n => initialize_d
    procedure :: set_ivariables, set_dvariables
    procedure :: get_ivariables, get_dvariables, get_dvariables_proc
    procedure :: clamp_variables
    procedure :: obj_1, obj_n

    procedure :: is_same => dis_same
    procedure :: print_u => print_d
    procedure :: print_up => print_dp
    procedure :: print_wv_u => print_wv_d
    procedure :: print_wv_up => print_wv_dp
    procedure :: print_header_u => print_header_d
    procedure :: print_header_wv_u => print_header_wv_d

    final :: destroy_instance
  end type TIndiv


  type, extends(TIndiv) :: TIndivI
    contains

    procedure :: initialize_n => initialize_i
    procedure :: is_same => iis_same
    procedure :: print_u => print_i
    procedure :: print_wv_u => print_wv_i
    procedure :: print_header_u => print_header_i
    procedure :: print_header_wv_u => print_header_wv_i
  end type TIndivI


  type, extends(TIndiv) :: TIndivS
    contains
  end type TIndivS


  type, extends(TIndivI) :: TIndivSI
    contains
  end type TIndivSI


  type, extends(TIndiv) :: TIndivC
    contains

    procedure :: print_u => print_dc
    procedure :: print_up => print_dcp
    procedure :: print_wv_u => print_wv_dc
    procedure :: print_wv_up => print_wv_dcp
    procedure :: print_header_u => print_header_dc
    procedure :: print_header_wv_u => print_header_wv_dc
  end type TIndivC


  type, extends(TIndivI) :: TIndivCI
    contains

    procedure :: print_u => print_ic
    procedure :: print_wv_u => print_wv_ic
    procedure :: print_header_u => print_header_ic
    procedure :: print_header_wv_u => print_header_wv_ic
  end type TIndivCI

  integer :: current_id = 0

  contains


  ! ============================================================================
  ! constructor
  ! ============================================================================

  subroutine set_id(this)
    implicit none
    class(TIndiv), intent(inout) :: this

    if (this%id <= 0) then
      current_id = current_id + 1
      this%id = current_id
    end if
  end subroutine set_id

  subroutine initialize_d(this, nvar)
    implicit none
    class(TIndiv), intent(inout) :: this
    integer, intent(in) :: nvar

    call this%set_id
    this%dvariables = random_array(nvar)
  end subroutine initialize_d

  subroutine initialize_i(this, nvar)
    implicit none
    class(TIndivI), intent(inout) :: this
    integer, intent(in) :: nvar

    call this%set_id
    this%ivariables = shuffle(nvar)
  end subroutine initialize_i


  ! ============================================================================
  ! setter
  ! ============================================================================

  subroutine set_ivariables(this, variables)
    implicit none
    class(TIndiv), intent(inout) :: this
    integer, intent(in) :: variables(:)

    call this%set_id
    this%ivariables = variables
  end subroutine set_ivariables

  subroutine set_dvariables(this, variables)
    implicit none
    class(TIndiv), intent(inout) :: this
    real(8), intent(in) :: variables(:)

    call this%set_id
    this%dvariables = variables
  end subroutine set_dvariables


  ! ============================================================================
  ! getter
  ! ============================================================================

  subroutine get_ivariables(this, variables)
    implicit none
    class(TIndiv), intent(inout) :: this
    integer, intent(out) :: variables(:)

    variables = this%ivariables
  end subroutine get_ivariables

  subroutine get_dvariables(this, variables)
    implicit none
    class(TIndiv), intent(inout) :: this
    real(8), intent(out) :: variables(:)

    variables = this%dvariables
  end subroutine get_dvariables

  subroutine get_dvariables_proc(this, variables, proc)
    implicit none
    class(TIndiv), intent(inout) :: this
    real(8), intent(out) :: variables(:)
    procedure(func_1d_1d), pointer :: proc

    if (associated(proc)) then
      variables = proc(this%dvariables)
    else
      variables = this%dvariables
    end if
  end subroutine get_dvariables_proc


  ! ============================================================================
  !
  ! ============================================================================

  subroutine clamp_variables(this, lower, upper)
    implicit none
    class(TIndiv), intent(inout) :: this
    real(8), intent(in) :: lower, upper

    call this%set_variables(min(max(this%dvariables, lower), upper))
  end subroutine clamp_variables

  real(8) function obj_1(this) result(res)
    implicit none
    class(TIndiv), intent(in) :: this

    res = this%obj(1)
  end function obj_1

  real(8) function obj_n(this, n) result(res)
    implicit none
    class(TIndiv), intent(in) :: this
    integer, intent(in) :: n

    res = this%objectives(n)
  end function obj_n


  ! ============================================================================
  !
  ! ============================================================================

  logical function dis_same(this, other) result(res)
    implicit none
    class(TIndiv), intent(in) :: this, other

    res = all(this%dvariables == other%dvariables)
  end function dis_same

  logical function iis_same(this, other) result(res)
    implicit none
    class(TIndivI), intent(in) :: this
    class(TIndiv), intent(in) :: other

    res = all(this%ivariables == other%ivariables)
  end function iis_same


  ! ============================================================================
  !
  ! ============================================================================

  subroutine print_d(this, unit)
    implicit none
    class(TIndiv), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives), "es15.8", ",") // ")"
    write(unit, format, advance='no') this%id, this%objectives
  end subroutine print_d

  subroutine print_dp(this, unit, proc)
    implicit none
    class(TIndiv), intent(in) :: this
    integer, intent(in) :: unit
    procedure(func_1d_1d) :: proc
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives), "es15.8", ",") // ")"
    write(unit, format, advance='no') this%id, this%objectives
  end subroutine print_dp

  subroutine print_i(this, unit)
    implicit none
    class(TIndivI), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives), "i0", ",") // ")"
    write(unit, format, advance='no') this%id, this%objectives
  end subroutine print_i

  subroutine print_dc(this, unit)
    implicit none
    class(TIndivC), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives) + &
                            size(this%constraints), "es15.8','") // "l0)"
    write(unit, format, advance='no') this%id, this%objectives, this%constraints, this%feasible
  end subroutine print_dc

  subroutine print_dcp(this, unit, proc)
    implicit none
    class(TIndivC), intent(in) :: this
    integer, intent(in) :: unit
    procedure(func_1d_1d) :: proc
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives) + &
                            size(this%constraints), "es15.8','") // "l0)"
    write(unit, format, advance='no') this%id, this%objectives, this%constraints, this%feasible
  end subroutine print_dcp

  subroutine print_ic(this, unit)
    implicit none
    class(TIndivCI), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives), "es15.8','") &
                 // nformat(size(this%constraints), "es15.8','") // "l0)"
    write(unit, format, advance='no') this%id, this%objectives, this%constraints, this%feasible
  end subroutine print_ic


  ! ============================================================================
  !
  ! ============================================================================

  subroutine print_wv_d(this, unit)
    implicit none
    class(TIndiv), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives) + &
                            size(this%dvariables), "es15.8", ",") // ")"
    write(unit, format, advance='no') this%id, this%objectives, this%dvariables
  end subroutine print_wv_d

  subroutine print_wv_dp(this, unit, proc)
    implicit none
    class(TIndiv), intent(in) :: this
    integer, intent(in) :: unit
    procedure(func_1d_1d) :: proc
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives) + &
                            size(this%dvariables), "es15.8", ",") // ")"
    write(unit, format, advance='no') this%id, this%objectives, proc(this%dvariables)
  end subroutine print_wv_dp

  subroutine print_wv_i(this, unit)
    implicit none
    class(TIndivI), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives), "es15.8','") &
                 // nformat(size(this%ivariables), "i0", ",") // ")"
    write(unit, format, advance='no') this%id, this%objectives, this%ivariables
  end subroutine print_wv_i

  subroutine print_wv_dc(this, unit)
    implicit none
    class(TIndivC), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives) + &
                            size(this%dvariables) + &
                            size(this%constraints), "es15.8','") // "l0)"
    write(unit, format, advance='no') this%id, this%objectives, this%dvariables, this%constraints, this%feasible
  end subroutine print_wv_dc

  subroutine print_wv_dcp(this, unit, proc)
    implicit none
    class(TIndivC), intent(in) :: this
    integer, intent(in) :: unit
    procedure(func_1d_1d) :: proc
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives) + &
                            size(this%dvariables) + &
                            size(this%constraints), "es15.8','") // "l0)"
    write(unit, format, advance='no') this%id, this%objectives, proc(this%dvariables), this%constraints, this%feasible
  end subroutine print_wv_dcp

  subroutine print_wv_ic(this, unit)
    implicit none
    class(TIndivCI), intent(in) :: this
    integer, intent(in) :: unit
    character(:), allocatable :: format

    format = "(i0','" // nformat(size(this%objectives), "es15.8','") &
                 // nformat(size(this%ivariables), "i0','")     &
                 // nformat(size(this%constraints), "es15.8','") // "l0)"
    write(unit, format, advance='no') this%id, this%objectives, this%ivariables, this%constraints, this%feasible
  end subroutine print_wv_ic


  ! ============================================================================
  !
  ! ============================================================================

  subroutine print_header_d(this, unit)
    implicit none
    class(TIndiv), intent(in) :: this
    integer, intent(in) :: unit
    integer :: l

    l = size(this%objectives)
    write(unit, "(a)", advance='no') "id," // join(serstr("obj", integers(l)), ",")
  end subroutine print_header_d

  subroutine print_header_i(this, unit)
    implicit none
    class(TIndivI), intent(in) :: this
    integer, intent(in) :: unit
    integer :: l

    l = size(this%objectives)
    write(unit, "(a)", advance='no') "id," // join(serstr("obj", integers(l)), ",")
  end subroutine print_header_i

  subroutine print_header_dc(this, unit)
    implicit none
    class(TIndivC), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(2)

    l(1) = size(this%objectives)
    l(2) = size(this%constraints)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
                                serstr("con", integers(l(2)))])
    write(unit, "(a)", advance='no') "id," // join(s, ",") // ",feasible"
  end subroutine print_header_dc

  subroutine print_header_ic(this, unit)
    implicit none
    class(TIndivCI), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(2)

    l(1) = size(this%objectives)
    l(2) = size(this%constraints)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
                                serstr("con", integers(l(2)))])
    write(unit, "(a)", advance='no') "id," // join(s, ",") // ",feasible"
  end subroutine print_header_ic


  ! ============================================================================
  !
  ! ============================================================================

  subroutine print_header_wv_d(this, unit)
    implicit none
    class(TIndiv), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(2)

    l(1) = size(this%objectives)
    l(2) = size(this%dvariables)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
                                serstr("var", integers(l(2)))])
    write(unit, "(a)", advance='no') "id," // join(s, ",")
  end subroutine print_header_wv_d

  subroutine print_header_wv_i(this, unit)
    implicit none
    class(TIndivI), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(2)

    l(1) = size(this%objectives)
    l(2) = size(this%ivariables)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
                                serstr("var", integers(l(2)))])
    write(unit, "(a)", advance='no') "id," // join(s, ",")
  end subroutine print_header_wv_i

  subroutine print_header_wv_dc(this, unit)
    implicit none
    class(TIndivC), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(3)

    l(1) = size(this%objectives)
    l(2) = size(this%dvariables)
    l(3) = size(this%constraints)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
                                serstr("var", integers(l(2))), &
                                serstr("con", integers(l(3)))])
    write(unit, "(a)", advance='no') "id," // join(s, ",") // ",feasible"
  end subroutine print_header_wv_dc

  subroutine print_header_wv_ic(this, unit)
    implicit none
    class(TIndivCI), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(3)

    l(1) = size(this%objectives)
    l(2) = size(this%ivariables)
    l(3) = size(this%constraints)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
                                serstr("var", integers(l(2))), &
                                serstr("con", integers(l(3)))])
    write(unit, "(a)", advance='no') "id," // join(s, ",") // ",feasible"
  end subroutine print_header_wv_ic


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_instance(this)
    implicit none
    type(TIndiv), intent(inout) :: this

    if (allocated(this%dvariables)) deallocate(this%dvariables)
    if (allocated(this%ivariables)) deallocate(this%ivariables)
    if (allocated(this%objectives)) deallocate(this%objectives)
    if (allocated(this%constraints)) deallocate(this%constraints)
  end subroutine destroy_instance
end module individual
