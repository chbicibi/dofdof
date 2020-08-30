module problem
  use interface
  use individual
  implicit none

  private
  public :: TProblem

  type :: TProblem
    character(:), allocatable :: ftype
    logical :: constraint, ready = .false.
    procedure(func_1d_1d), nopass, pointer :: objective_func
    procedure(func_1d_1d), nopass, pointer :: constraint_func
    procedure(func_1d_1d), nopass, pointer :: scaling_func

    contains

    generic :: initialize => initialize_func
    generic :: call => call_indiv,                       &
                       call_func_so, call_func_so_con,   &
                       call_func_mo, call_func_mo_con,   &
                       call_ifunc_so, call_ifunc_so_con, &
                       call_ifunc_mo, call_ifunc_mo_con

    procedure :: initialize_func
    procedure :: prepare
    procedure :: call_indiv
    procedure :: call_func_so, call_func_so_con
    procedure :: call_func_mo, call_func_mo_con
    procedure :: call_ifunc_so, call_ifunc_so_con
    procedure :: call_ifunc_mo, call_ifunc_mo_con
  end type TProblem

  real(8), parameter :: PI = 4d0 * atan(1d0)

  contains

  ! ============================================================================
  ! constructor
  ! ============================================================================

  subroutine initialize_func(this, ftype)
    implicit none
    class(TProblem), intent(inout) :: this
    character(*), intent(in) :: ftype
    procedure(func_1d_1d), pointer :: fobj, fcon, fscl

    nullify(fobj, fcon, fscl)

    select case (ftype)
    case ("ZDT1")
      fobj => testProblem_ZDT1
    case ("ZDT2")
      fobj => testProblem_ZDT2
    case ("ZDT3")
      fobj => testProblem_ZDT3
    case ("ZDT4")
      fobj => testProblem_ZDT4
    case ("ZDT6")
      fobj => testProblem_ZDT6
    case ("DTLZ1")
      fobj => testProblem_DTLZ1
    case ("DTLZ2")
      fobj => testProblem_DTLZ2
    case ("O_K")
      fobj => testProblem_O_K
      fcon => testProblem_O_K_constraints
      fscl => testProblem_O_K_scaling
    case ("TNK")
      fobj => testProblem_TNK
      fcon => testProblem_TNK_constraints
      fscl => testProblem_TNK_scaling
    case ("USER")
    case default
      write(0, "(3a/a)") "Error: unknown keyword '", ftype, "'", &
                         "in subroutine 'TProblem%initialize_func'"
      call exit(1)
    end select

    this%ftype = ftype
    this%objective_func => fobj
    this%constraint_func => fcon
    this%scaling_func => fscl
  end subroutine initialize_func

  subroutine prepare(this)
    implicit none
    class(TProblem), intent(inout) :: this

    if (.not. associated(this%objective_func)) then
      write(0, "(a/a)") "Error: function 'objective_func' is associated", &
                        "in subroutine 'TProblem%prepare'"
      call exit(1)
    end if

    if (.not. associated(this%constraint_func)) this%constraint_func => empty
    if (.not. associated(this%scaling_func)) this%scaling_func => pass

    this%ready = .true.
  end subroutine prepare


  ! ============================================================================
  ! call procedure
  ! ============================================================================

  subroutine call_indiv(this, indiv)
    implicit none
    class(TProblem), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv

    if (indiv%evaluated) return

    select type(indiv)
    class is (TIndiv)
      call this%call(indiv%dvariables, indiv%objectives)
    class is (TIndivI)
      call this%call(indiv%ivariables, indiv%objectives)
    class is (TIndivS)
      allocate(indiv%objectives(1))
      call this%call(indiv%dvariables, indiv%objectives(1), indiv%feasible)
    class is (TIndivSI)
      allocate(indiv%objectives(1))
      call this%call(indiv%ivariables, indiv%objectives(1), indiv%feasible)
    class is (TIndivC)
      call this%call(indiv%dvariables, indiv%objectives, indiv%constraints, indiv%feasible)
    class is (TIndivCI)
      call this%call(indiv%ivariables, indiv%objectives, indiv%constraints, indiv%feasible)
    end select

    indiv%evaluated = .true.
  end subroutine call_indiv

  subroutine call_func_so(this, variables, objective)
    implicit none
    class(TProblem), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out) :: objective

    if (.not. this%ready) call this%prepare
    associate (obj=>this%objective_func(this%scaling_func(variables)))
      objective = obj(1)
    end associate
  end subroutine call_func_so

  subroutine call_func_so_con(this, variables, objective, feasible)
    implicit none
    class(TProblem), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out) :: objective
    logical, intent(out) :: feasible

    if (.not. this%ready) call this%prepare
    associate (x=>this%scaling_func(variables))
      associate (obj=>this%objective_func(x), &
                 con=>this%constraint_func(x))
        objective = obj(1)
        feasible = all(con >= 0)
      end associate
    end associate
  end subroutine call_func_so_con

  subroutine call_func_mo(this, variables, objectives)
    implicit none
    class(TProblem), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out), allocatable :: objectives(:)

    if (.not. this%ready) call this%prepare
    objectives = this%objective_func(this%scaling_func(variables))
  end subroutine call_func_mo

  subroutine call_func_mo_con(this, variables, objectives, constraints, feasible)
    implicit none
    class(TProblem), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out), allocatable :: objectives(:), constraints(:)
    logical, intent(out) :: feasible

    if (.not. this%ready) call this%prepare
    associate (x=>this%scaling_func(variables))
      objectives = this%objective_func(x)
      constraints = this%constraint_func(x)
      feasible = all(constraints >= 0)
    end associate
  end subroutine call_func_mo_con

  subroutine call_ifunc_so(this, variables, objective)
    implicit none
    class(TProblem), intent(inout) :: this
    integer, intent(in) :: variables(:)
    real(8), intent(out) :: objective

    objective = 0d0
  end subroutine call_ifunc_so

  subroutine call_ifunc_so_con(this, variables, objective, feasible)
    implicit none
    class(TProblem), intent(inout) :: this
    integer, intent(in) :: variables(:)
    real(8), intent(out) :: objective
    logical, intent(out) :: feasible

    objective = 0d0
    feasible = .true.
  end subroutine call_ifunc_so_con

  subroutine call_ifunc_mo(this, variables, objectives)
    implicit none
    class(TProblem), intent(inout) :: this
    integer, intent(in) :: variables(:)
    real(8), intent(out), allocatable :: objectives(:)

    allocate(objectives(1), source=0d0)
  end subroutine call_ifunc_mo

  subroutine call_ifunc_mo_con(this, variables, objectives, constraints, feasible)
    implicit none
    class(TProblem), intent(inout) :: this
    integer, intent(in) :: variables(:)
    real(8), intent(out), allocatable :: objectives(:), constraints(:)
    logical, intent(out) :: feasible

    allocate(objectives(1), source=0d0)
    allocate(constraints(1), source=0d0)
    feasible = .true.
  end subroutine call_ifunc_mo_con


  ! ============================================================================
  ! dummy functions
  ! ============================================================================

  function pass(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    res = variables
  end function pass

  function empty(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    allocate(res(0))
  end function empty


  ! ============================================================================
  ! test functions body
  ! ============================================================================

  function testProblem_ZDT1(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: g

    g = 1d0 + 9d0 * sum(variables(2:)) / (size(variables) - 1)
    res = [variables(1), g * (1d0 - sqrt(variables(1) / g))]
  end function testProblem_ZDT1

  ! ******************************************************************

  function testProblem_ZDT2(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: g

    g = 1d0 + 9d0 * sum(variables(2:)) / (size(variables) - 1)
    res = [variables(1), g * (1d0 - (variables(1) / g) ** 2)]
  end function testProblem_ZDT2

  ! ******************************************************************

  function testProblem_ZDT3(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: g, t

    g = 1d0 + 9d0 * sum(variables(2:)) / (size(variables) - 1)
    t = variables(1) / g
    res = [variables(1), g * (1d0 - sqrt(t) - t * sin(10d0 * PI * variables(1)))]
  end function testProblem_ZDT3

  ! ******************************************************************

  function testProblem_ZDT4(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: g

    g = 1d0 + 10d0 * (size(variables) - 1) + &
        sum(variables(2:) ** 2 - 10d0 * cos(4d0 * PI * variables(2:)))
    res = [variables(1), g * (1d0 - sqrt(variables(1) / g))]
  end function testProblem_ZDT4

  ! ******************************************************************

  function testProblem_ZDT6(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: g

    g = 1d0 + 9d0 * (sum(variables(2:)) / (size(variables) - 1)) ** 0.25d0

    res = [1d0 - exp(-4d0 * variables(1)) * sin(6d0 * PI * variables(1)) ** 6, 0d0]
    res(2) = g * (1d0 - (res(1) / g) ** 2)
  end function testProblem_ZDT6

  ! ******************************************************************

  function testProblem_DTLZ1(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: g

    g = 100d0 * ((size(variables) - 2) + &
                 (sum((variables(3:) - 0.5d0) ** 2 - cos(20d0 * PI * (variables(3:) - 0.5d0)))))

    res = (1d0 + g) * [variables(1) * variables(2),         &
                       variables(1) * (1d0 - variables(2)), &
                       (1d0 - variables(1))]
  end function testProblem_DTLZ1

  ! ******************************************************************

  function testProblem_DTLZ2(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: g, pi_h = 0.5d0 * PI

    g = sum((variables(3:) ** 2))

    res = (1d0 + g) * [cos(variables(1) * pi_h) * cos(variables(2) * pi_h), &
                       cos(variables(1) * pi_h) * sin(variables(2) * pi_h), &
                       sin(variables(1) * pi_h)]
  end function testProblem_DTLZ2

  ! ******************************************************************

  ! Osyczka and Kundu function
  function testProblem_O_K(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    res = [-sum([25d0, 1d0, 1d0, 1d0, 1d0] * (variables(:5) - [2d0, 2d0, 1d0, 4d0, 1d0]) ** 2), &
            sum(variables ** 2)]
  end function testProblem_O_K

  function testProblem_O_K_constraints(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    res = [-2d0 + variables(1) + variables(2),              &
            6d0 - variables(1) - variables(2),              &
            2d0 + variables(1) - variables(2),              &
            2d0 - variables(1) + 3d0 * variables(2),        &
            4d0 - (variables(3) - 3d0) ** 2 - variables(4), &
           -4d0 + (variables(5) - 3d0) ** 2 + variables(6)]
  end function testProblem_O_K_constraints

  function testProblem_O_K_scaling(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: low(6), upp(6)

    low = [0d0, 0d0, 1d0, 0d0, 1d0, 0d0]
    upp = [10d0, 10d0, 5d0, 6d0, 5d0, 10d0]
    res = low + (upp - low) * variables
  end function testProblem_O_K_scaling

  ! ******************************************************************

  function testProblem_TNK(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    res = variables
  end function testProblem_TNK

  function testProblem_TNK_constraints(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    associate (x => variables)
      res = [sum(x ** 2) - 1d0 - 0.1d0 * cos(16d0 * atan(x(1) / x(2))), &
             0.5d0 - sum((x - 0.5d0) ** 2)]
    end associate
  end function testProblem_TNK_constraints

  function testProblem_TNK_scaling(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: low(2), upp(2)

    low = [0d0, 0d0]
    upp = [PI, PI]
    res = low + (upp - low) * variables
  end function testProblem_TNK_scaling


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_pronlem(this)
    implicit none
    type(TProblem), intent(inout) :: this

    if (associated(this%objective_func)) nullify(this%objective_func)
    if (associated(this%constraint_func)) nullify(this%constraint_func)
    if (associated(this%scaling_func)) nullify(this%scaling_func)
  end subroutine destroy_pronlem
end module problem
