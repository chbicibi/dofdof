module nsga2i
  use individual
  use nsga2
  implicit none

  private
  public :: TNSGA2I

  type, extends(TNSGA2) :: TNSGA2I
    contains

    procedure :: set_prototype1
    procedure :: reproduce
  end type TNSGA2I

  contains


  subroutine set_prototype1(this)
    implicit none
    class(TNSGA2I), intent(inout) :: this

    if (.not. allocated(this%prototype)) allocate(TIndivI::this%prototype)
  end subroutine set_prototype1

  subroutine reproduce(this, index, children)
    implicit none
    class(TNSGA2I), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: children(:)
    integer, allocatable :: parents_value(:, :), children_value(:, :)
    integer :: i

    parents_value = reshape([(this%population(index(i))%indiv%ivariables, i = 1, 2)], &
                            [this%num_var, 2])

    call this%crossover%call(parents_value, children_value)
    allocate(children(size(children_value, dim=2)), source=this%prototype)

    do i = 1, size(children)
      call children(i)%set_variables(children_value(:, i))
      call this%mutation%call(children(i)%ivariables)
    end do
  end subroutine reproduce
end module nsga2i
