module rga
  use interface
  use util
  use individual
  use moead_unit
  use ga_unit
  use problem
  use soga
  use nsga2c
  implicit none

  private
  public :: TRGA

  type, extends(TNSGA2C) :: TRGA
    integer :: subpop_size
    class(TRGA), pointer :: subpop(:)

    contains

    generic :: save_result => save_result_con
    generic :: save_history => save_history_con

    procedure :: run
    procedure :: save_result_con => save_result
    procedure :: save_history_con => save_history

    procedure :: evaluate_single
    procedure :: calc_fitness
    procedure :: count_feasible
  end type TRGA

  contains


  ! ============================================================================
  ! calculation body
  ! ============================================================================

  subroutine run(this, num_generation)
    implicit none
    class(TRGA), intent(inout) :: this
    integer, intent(in) :: num_generation
    integer :: offset, i

    call this%alloc_history(num_generation, offset)

    print "(2a/)", "Start: ", this%problem_type

    do i = 1, num_generation
      call this%evolution
      this%history(:, i + offset) = this%population
      if (mod(i, 50) == 0 .and. i > 0) then
        print "(a, i0, a)", "Progress: ", i, " steps finished"
        call this%count_feasible
      end if
    end do

    print "(2a/)", "End: ", this%problem_type
  end subroutine run

  subroutine evaluate_single(this, indiv)
    implicit none
    class(TRGA), intent(inout) :: this
    type(TIndiv), intent(inout) :: indiv

    call this%problem%call(indiv%variables, indiv%objectives, indiv%constraints, indiv%feasible)
  end subroutine evaluate_single

  subroutine calc_fitness(this, population)
    implicit none
    class(TRGA), intent(in) :: this
    type(TPopulation), intent(inout) :: population(:)
    real(8), allocatable :: objectives(:, :), constraints(:)
    logical, allocatable :: feasible(:)
    real(8) :: a = 0.1d0
    integer :: i

    objectives = reshape([(population(i)%indiv%objectives, i = 1, this%pop_size)], &
                         [this%num_obj, this%pop_size])
    constraints = [(scalar_constraints(population(i)%indiv%constraints), i = 1, this%pop_size)]
    feasible = population%indiv%feasible
    population%rank = rank_pareto(objectives, constraints, feasible, dominated_p)

    if (this%sharing) then
      ! population%crowding = crowding_distance(objectives, population%rank)
      ! population%crowding = min(crowding_distance(objectives, population%rank), 1d0)
      population%crowding = tanh(crowding_distance(objectives, population%rank))
      population%fitness = a * (1.0d0 - a) ** (population%rank - population%crowding)
      ! population%fitness = 1.0d0 / (population%rank + 1 - population%crowding)
    else
      population%fitness = a * (1.0d0 - a) ** (population%rank - 1)
      ! population%fitness = 1.0d0 / population%rank
    end if
  end subroutine calc_fitness

  logical function dominated_p(val, con, feasible) result(res)
    implicit none
    real(8), intent(in) :: val(:, :), con(:)
    logical, intent(in) :: feasible(:)

    res = (.not. feasible(1) .and. feasible(2)) .or.                        &
          (.not. (feasible(1) .or. feasible(2)) .and. con(1) > con(2)) .or. &
          (feasible(1) .and. feasible(2) .and. all(val(:, 1) > val(:, 2)))
  end function dominated_p


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine save_result(this, filename, mode1, mode2)
    implicit none
    class(TRGA), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: mode1, mode2
    character(:), allocatable :: format
    type(string), allocatable :: s(:)
    integer :: unit, i

    format = "(" // str(this%num_obj) // "(e15.5',')"    &
                 // str(this%num_var) // "(','e15.5)','" &
                 // str(this%num_con) // "(','e15.5)',,'i0,2(','f0.5))"

    open(newunit=unit, file=filename)
      s = [(string("obj" // str(i)), i = 1, this%num_obj), string(""), &
           (string("var" // str(i)), i = 1, this%num_var), string(""), &
           (string("con" // str(i)), i = 1, this%num_con), string(""), &
           string("rank"), string("fitness"), string("crowding-dist")]
      write(unit, "(a)") join(s, ",")

      do i = 1, this%pop_size
        if ((mode1 == 0 .or. xor(mode1 == 2, this%population(i)%rank == 1)) .and. &
            (mode2 == 0 .or. xor(mode2 == 2, this%population(i)%indiv%feasible))) &
              write(unit, format) this%population(i)%indiv%objectives,            &
                                  this%dec(this%population(i)%indiv%variables),   &
                                  this%population(i)%indiv%constraints,           &
                                  this%population(i)%rank,                        &
                                  this%population(i)%fitness,                     &
                                  this%population(i)%crowding
      end do
    close(unit)
  end subroutine save_result

  subroutine save_history(this, filename, mode1, mode2)
    implicit none
    class(TRGA), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: mode1, mode2
    character(:), allocatable :: format
    type(string), allocatable :: s(:)
    integer :: unit, i, j

    format = "(" // str(this%num_obj) // "(e15.5',')"    &
                 // str(this%num_var) // "(','e15.5)','" &
                 // str(this%num_con) // "(','e15.5)',,'i0,2(','f0.5))"

    open(newunit=unit, file=filename)
      s = [(string("obj" // str(i)), i = 1, this%num_obj), string(""), &
           (string("var" // str(i)), i = 1, this%num_var), string(""), &
           (string("con" // str(i)), i = 1, this%num_con), string(""), &
           string("rank"), string("fitness"), string("crowding-dist")]
      write(unit, "(a)") join(s, ",")

      do j = 1, size(this%history, dim=2)
        do i = 1, this%pop_size
          if ((mode1 == 0 .or. xor(mode1 == 2, this%history(i, j)%rank == 1)) .and. &
              (mode2 == 0 .or. xor(mode2 == 2, this%history(i, j)%indiv%feasible))) &
                write(unit, format) this%history(i, j)%indiv%objectives,            &
                                    this%dec(this%history(i, j)%indiv%variables),   &
                                    this%history(i, j)%indiv%constraints,           &
                                    this%history(i, j)%rank,                        &
                                    this%history(i, j)%fitness,                     &
                                    this%history(i, j)%crowding
        end do
      end do
    close(unit)
  end subroutine save_history


  ! ============================================================================
  ! scalar function (temp)
  ! ============================================================================

  real(8) function scalar_constraints(constraints) result(res)
    implicit none
    real(8), intent(in) :: constraints(:)

    res = -sum(min(constraints, 0d0))
  end function scalar_constraints


  ! ============================================================================
  ! other
  ! ============================================================================

  subroutine count_feasible(this)
    implicit none
    class(TRGA), intent(in) :: this

    print "('  ->feasible: 'i0'/'i0)", count(this%population%indiv%feasible), this%pop_size
  end subroutine count_feasible
end module rga
