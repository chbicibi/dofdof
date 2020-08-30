module basic_optimizer
  use interface
  use individual
  use problem
  implicit none

  private
  public :: TOptimizer

  type :: TOptimizer
    integer :: current_step = 0
    integer :: start_time = 0
    integer :: num_var
    integer :: num_obj
    integer :: num_con
    integer :: pop_size

    logical :: show_log = .true.

    character(:), allocatable :: problem_type

    class(TProblem), allocatable :: problem
    class(TIndiv), allocatable :: prototype

    contains

    generic :: set_problem => set_problem1, set_problem2
    generic :: set_prototype => set_prototype1, set_prototype2
    generic :: evaluate => evaluate_single

    procedure :: set_problem1, set_problem2
    procedure :: set_prototype1, set_prototype2
    procedure :: set_objective_function, set_constraint_function, set_scaling_function
    procedure :: init_indiv

    procedure :: evaluate_single

    procedure :: logger
    procedure :: print_indiv
    procedure :: dec

    final :: destroy_instance
  end type TOptimizer

  contains


  subroutine set_problem1(this)
    implicit none
    class(TOptimizer), intent(inout) :: this

    allocate(TProblem::this%problem)
    call this%problem%initialize(this%problem_type)
  end subroutine set_problem1

  subroutine set_problem2(this, problem)
    implicit none
    class(TOptimizer), intent(inout) :: this
    class(TProblem), intent(in) :: problem

    if (allocated(this%problem)) deallocate(this%problem)
    allocate(this%problem, source=problem)
    this%problem_type = "USER"
  end subroutine set_problem2

  subroutine set_prototype1(this)
    implicit none
    class(TOptimizer), intent(inout) :: this

    if (.not. allocated(this%prototype)) allocate(TIndivS::this%prototype)
  end subroutine set_prototype1

  subroutine set_prototype2(this, indiv)
    implicit none
    class(TOptimizer), intent(inout) :: this
    class(TIndiv), intent(in) :: indiv

    if (allocated(this%prototype)) deallocate(this%prototype)
    allocate(this%prototype, source=indiv)
  end subroutine set_prototype2

  subroutine set_objective_function(this, proc)
    implicit none
    class(TOptimizer), intent(inout) :: this
    procedure(func_1d_1d) :: proc

    this%problem%objective_func => proc
    this%problem_type = "USER"
  end subroutine set_objective_function

  subroutine set_constraint_function(this, proc, ncon)
    implicit none
    class(TOptimizer), intent(inout) :: this
    integer, intent(in) :: ncon
    procedure(func_1d_1d) :: proc

    this%problem%constraint_func => proc
    this%num_con = ncon
  end subroutine set_constraint_function

  subroutine set_scaling_function(this, proc)
    implicit none
    class(TOptimizer), intent(inout) :: this
    procedure(func_1d_1d) :: proc

    this%problem%scaling_func => proc
  end subroutine set_scaling_function

  subroutine init_indiv(this, indiv)
    implicit none
    class(TOptimizer), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv

    call indiv%initialize(this%num_var)
  end subroutine init_indiv

  subroutine evaluate_single(this, indiv)
    implicit none
    class(TOptimizer), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv

    if (indiv%evaluated) return
    call this%problem%call(indiv)
    indiv%evaluated = .true.
  end subroutine evaluate_single


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine logger(this, n, total)
    implicit none
    class(TOptimizer), intent(in) :: this
    integer, intent(in) :: n, total

    if (.not. this%show_log) return

    if (n == 0) then
      print "(2a)", "Start: ", this%problem_type
    else if (n == -1) then
      print "(2a/)", "End: ", this%problem_type
    else if (mod(n, 1) == 0) then
      print "(a i0 a)", "Progress: ", n, " steps finished"
    end if
  end subroutine logger

  subroutine print_indiv(this, indiv, unit, with_var)
    implicit none
    class(TOptimizer), intent(in) :: this
    class(TIndiv), intent(in) :: indiv
    integer, intent(in) :: unit
    logical, intent(in), optional :: with_var

    if (present(with_var) .and. with_var) then
      if (associated(this%problem%scaling_func)) then
        call indiv%print_wv(unit, this%problem%scaling_func)
      else
        call indiv%print_wv(unit)
      end if
    else
      if (associated(this%problem%scaling_func)) then
        call indiv%print(unit, this%problem%scaling_func)
      else
        call indiv%print(unit)
      end if
    end if
  end subroutine print_indiv

  function dec(this, variables) result(res)
    implicit none
    class(TOptimizer), intent(in) :: this
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    if (associated(this%problem%scaling_func)) then
      res = this%problem%scaling_func(variables)
    else
      res = variables
    end if
  end function dec


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_instance(this)
    implicit none
    type(TOptimizer), intent(inout) :: this

    if (allocated(this%problem)) deallocate(this%problem)
    if (allocated(this%prototype)) deallocate(this%prototype)
  end subroutine destroy_instance
end module basic_optimizer
