module tnsdm
  use interface
  use util
  use individual
  use problem
  use ga_unit
  use soga
  use nsga2c
  implicit none

  private
  public :: TTNSDM

  type, extends(TNSGA2C) :: TTNSDM
    integer, allocatable :: rank_con(:), rank_obj(:)
    logical, allocatable :: feasible(:), dominated(:, :), mask(:)

    contains

    procedure :: calc_fitness
    procedure :: select_parents
  end type TTNSDM

  contains


  ! ============================================================================
  ! calculation body
  ! ============================================================================

  subroutine calc_fitness(this, population)
    implicit none
    class(TTNSDM), intent(inout) :: this
    type(TPopulation), intent(inout) :: population(:)
    real(8), allocatable :: objectives(:, :), constraints(:, :)
    logical, allocatable :: feasible(:)!, dominated(:, :), mask(:)
    ! integer, allocatable :: rank_con(:), rank_obj(:), order(:), index(:)
    integer, allocatable :: rank(:), order(:)
    ! integer, allocatable :: rank2(:)
    real(8) :: a = 0.1d0
    integer :: pop_size, i, j

    pop_size = size(population)

    ! ------------------------------------------------------
    ! 作業用変数割付
    ! ------------------------------------------------------
    objectives = reshape([(population(i)%indiv%objectives, i = 1, this%pop_size)], &
                         [this%num_obj, this%pop_size])
    constraints = reshape([(max(-population(i)%indiv%constraints, 0d0), i = 1, this%pop_size)], &
                          [this%num_obj, this%pop_size])
    feasible = [(this%population(i)%indiv%feasible, i = 1, this%pop_size)]

    this%rank_con = rank_pareto(constraints)
    this%rank_obj = rank_pareto(objectives, this%rank_con)

    associate (n => this%pop_size, val => objectives)
      this%dominated = reshape([((dominatedp1(val(:, i), val(:, j)), j = 1, n), i = 1, n)], [n, n])
    end associate


    allocate(rank(this%pop_size), source=this%rank_con*this%pop_size+this%rank_obj)
    allocate(order(this%pop_size), source=sort(rank))
    order(order) = integers(this%pop_size)

    ! print "(100i5/)", this%rank_con
    ! print "(100i5/)", this%rank_obj
    ! print "(100i5/)", rank_pareto(objectives)
    ! ! print "(100i5/)", rank
    ! print "(100i5/)", order

    ! rank2 = rank_pareto(objectives, constraints, feasible, dominatedp)
    ! print *, feasible
    ! print "(100i5/)", rank2
    ! stop

    if (this%sharing) then
      ! population%crowding = crowding_distance(objectives, rank)
      ! population%crowding = min(crowding_distance(objectives, rank), 1d0)
      population%crowding = tanh(crowding_distance(objectives, rank))
      population%fitness = a * (1.0d0 - a) ** (order - population%crowding)
      ! population%fitness = 1.0d0 / (order + 1 - population%crowding)
    else
      population%fitness = a * (1.0d0 - a) ** (order - 1)
      ! population%fitness = 1.0d0 / order
    end if

    population%rank = rank_pareto(objectives, constraints, feasible, dominatedp)

    if (allocated(this%mask)) then
      this%mask = .true.
    else
      allocate(this%mask(this%pop_size), source=.true.)
    end if

    where (order > this%pop_size / 2) this%mask = .false.
  end subroutine calc_fitness

  subroutine select_parents(this, parents_index)
    implicit none
    class(TTNSDM), intent(in) :: this
    integer, intent(out), allocatable :: parents_index(:)
    ! integer, allocatable :: order(:), index(:)
    logical, allocatable :: mask1(:), mask2(:)
    integer :: parent1, parent2

    ! allocate(order(this%pop_size), source=sort(this%rank_con * this%pop_size + this%rank_obj))

    allocate(mask1(this%pop_size), source=this%mask)

    ! ------------------------------------------------------
    ! 親1を選択
    ! ------------------------------------------------------
    call this%selection%call(this%population%fitness, mask1, parent1)

    ! parent1 = order(random(this%pop_size / 5)) ! => ランキング or トーナメント
    ! print *, "parent0", order(1)
    ! print *, objectives(:, order(1))
    ! print *, constraints(:, order(1))
    ! print *, "parent1", parent1
    ! print *, objectives(:, parent1)
    ! print *, constraints(:, parent1)


    ! ------------------------------------------------------
    ! 親2を選択
    ! ------------------------------------------------------

    mask2 = this%dominated(:, parent1)
    if (count(mask2) < 1) then
      call this%selection%call(this%population%fitness, mask1, parent2)
    else
      call this%selection%call(this%population%fitness, mask2, parent2)
    end if

    ! index = mask_index(mask)
    ! parent2 = index(random(size(index))) ! ランキング or トーナメント
    ! print *, "parent2", parent2
    ! stop
    if (parent1 == parent2) then
      print *, "parent1 == parent2"
      stop
    end if
    parents_index = [parent1, parent2]
  end subroutine select_parents

  subroutine select_neighbor(this, index1, index2)
    implicit none
    class(TTNSDM), intent(in) :: this
    integer, intent(in) :: index1
    integer, intent(out) :: index2
    logical, allocatable :: mask(:)
    integer :: i

    mask = [(this%population(i)%indiv%feasible, i = 1, this%pop_size)] .and. &
           this%population%rank == this%population(index1)%rank
    if (count(mask) < 1) then
      index2 = 0
    else
      call this%selection%call(this%population%fitness, mask, index2)
    end if
  end subroutine select_neighbor

  logical function dominatedp(val, con, feasible) result(res)
    implicit none
    real(8), intent(in) :: val(:, :), con(:, :)
    logical, intent(in) :: feasible(:)

    res = (.not. feasible(1) .and. feasible(2)) .or.                                          &
          (.not. (feasible(1) .or. feasible(2)) .and. dominatedp1(con(:, 1), con(:, 2))) .or. &
          (feasible(1) .and. feasible(2) .and. dominatedp1(val(:, 1), val(:, 2)))
  end function dominatedp

  logical function dominatedp1(val1, val2, mode) result(res)
    implicit none
    real(8), intent(in) :: val1(:), val2(:)
    integer, intent(in), optional :: mode

    if (present(mode) .and. mode > 0) then
      res = all(val1 > val2)
    else
      res = all(val1 >= val2) .and. any(val1 /= val2)
    end if
  end function dominatedp1


  ! ============================================================================
  ! scalar function (temp)
  ! ============================================================================

  real(8) function scalar_constraints(constraints) result(res)
    implicit none
    real(8), intent(in) :: constraints(:)

    res = -sum(min(constraints, 0d0))
  end function scalar_constraints

end module tnsdm

