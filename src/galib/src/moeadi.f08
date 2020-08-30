module moeadi
  use interface
  use util
  use individual
  use problem
  use moead_unit
  use ga_unit
  use moead
  implicit none

  private
  public :: TMOEADI

  type, extends(TMOEAD) :: TMOEADI
    contains

    procedure :: set_prototype1
    procedure :: reproduce
  end type TMOEADI

  contains

  subroutine set_prototype1(this)
    implicit none
    class(TMOEADI), intent(inout) :: this

    if (.not. allocated(this%prototype)) allocate(TIndivI::this%prototype)
  end subroutine set_prototype1

  subroutine reproduce(this, index, child)
    implicit none
    class(TMOEADI), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: child
    integer, allocatable :: parents_value(:, :), children_value(:, :)
    integer :: i

    allocate(child, source=this%prototype)

    parents_value = reshape([(this%population(index(i))%indiv%ivariables, i = 1, 2)], &
                            [this%num_var, 2])

    call this%crossover%call(parents_value, children_value)
    call child%set_variables(children_value(:, random(2)))
    call this%mutation%call(child%ivariables)
  end subroutine reproduce
end module moeadi
