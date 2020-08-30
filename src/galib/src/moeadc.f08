module moeadc
  use interface
  use util
  use individual
  use problem
  use moead_unit
  use ga_unit
  use moead
  implicit none

  private
  public :: TMOEADC

  type, extends(TMOEAD) :: TMOEADC
    contains

    generic :: save_result => save_result_feasible
    generic :: save_history => save_history_feasible

    procedure :: logger
    procedure :: save_result_feasible
    procedure :: save_history_feasible
    procedure :: update_population
  end type TMOEADC

  contains


  ! ============================================================================
  ! calculation body
  ! ============================================================================

  subroutine update_population(this, indiv, table, idx)
    implicit none
    class(TMOEADC), intent(inout) :: this
    class(TIndiv), intent(in) :: indiv
    integer, intent(in) :: table(:), idx
    real(8) :: f_old, f_new
    logical :: flag
    integer :: i, j

    associate (y => indiv)
      do i = 1, size(table)
        flag = .false.
        j = table(i)
        if (j /= idx) cycle
        associate (x => this%population(j)%indiv)
          if (y%feasible .and. .not. x%feasible) then
            flag = .true.
          else if (.not. (y%feasible .or. x%feasible)) then
            f_old = scalar_constraints(y%constraints)
            f_new = scalar_constraints(x%constraints)
            flag = f_new < f_old
          else if (y%feasible .and. x%feasible) then
            f_old = this%scalar_value(x, j)
            f_new = this%scalar_value(y, j)
            flag = f_new < f_old
          end if
          if (flag) x = y
        end associate
      end do
    end associate
  end subroutine update_population


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine logger(this, n, total)
    implicit none
    class(TMOEADC), intent(in) :: this
    integer, intent(in) :: n, total

    if (n == 0) then
      print "(2a)", "Start: ", this%problem_type
    else if (n == -1) then
      print "(2a/)", "End: ", this%problem_type
    else if (mod(n, 100) == 0) then
      print "(a i0 a)", "Progress: ", n, " steps finished"
      print "('  ->feasible: 'i0'/'i0)", count_feasible(this%population), this%pop_size
    end if
  end subroutine logger

  subroutine save_result_feasible(this, filename, feasible)
    implicit none
    class(TMOEADC), intent(in) :: this
    character(*), intent(in) :: filename, feasible
    integer :: unit, i

    open(newunit=unit, file=filename)
      call this%population(1)%indiv%print_header_wv(unit)
      write(unit, *)

      do i = 1, this%pop_size
        if (feasible == "all" .or. xor(feasible == "only", .not. this%population(i)%indiv%feasible)) then
          call this%print_indiv(this%population(i)%indiv, unit, .true.)
          write(unit, *)
        end if
      end do
    close(unit)
  end subroutine save_result_feasible

  subroutine save_history_feasible(this, filename, feasible)
    implicit none
    class(TMOEADC), intent(in) :: this
    character(*), intent(in) :: filename, feasible
    integer :: unit, i, j

    open(newunit=unit, file=filename)
      write(unit, "(a)", advance='no') "step,"
      call this%population(1)%indiv%print_header(unit)
      write(unit, *)

      do j = 1, size(this%history, dim=2)
        do i = 1, this%pop_size
          if (feasible == "all" .or. xor(feasible == "only", .not. this%history(i, j)%indiv%feasible)) then
            write(unit, "(i0',')", advance='no') j - 1
            call this%print_indiv(this%history(i, j)%indiv, unit)
            write(unit, *)
          end if
        end do
      end do
    close(unit)
  end subroutine save_history_feasible


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

  integer function count_feasible(population) result(res)
    implicit none
    type(TPopulation), intent(in) :: population(:)
    integer :: i

    res = 0
    do i = 1, size(population)
      if (population(i)%indiv%feasible) res = res + 1
    end do
  end function count_feasible
end module moeadc
