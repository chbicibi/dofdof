module nsga2c_dp
  use interface
  use util
  use individual
  use problem
  use ga_unit
  use soga
  use nsga2
  use nsga2c
  implicit none

  private
  public :: TNSGA2C_DP

  type, extends(TNSGA2C) :: TNSGA2C_DP
    contains

    procedure :: preserve_elite
    procedure :: evolution
  end type TNSGA2C_DP

  contains

  subroutine preserve_elite(this, new_population, num_elite)
    implicit none
    class(TNSGA2C_DP), intent(in) :: this
    type(TPopulation), intent(inout) :: new_population(:)
    integer, intent(out) :: num_elite
    integer :: i

    if (this%elite_preservation) then
      ! num_elite = count(this%population%rank == 1)
      num_elite = 0
      do i = 1, this%pop_size
        if (this%population(i)%rank /= 1) cycle
        if (.not. this%population(i)%indiv%feasible) cycle
        if (num_elite > 0 .and. .not. this%check_uniq(this%population(i)%indiv, new_population(1:num_elite))) cycle
        num_elite = num_elite + 1
        new_population(num_elite) = this%population(i)
      end do
      ! new_population(1:num_elite) = pack(this%population, this%population%rank == 1)
    else
      num_elite = 0
    end if
  end subroutine preserve_elite

  subroutine evolution(this)
    implicit none
    class(TNSGA2C_DP), intent(inout) :: this
    type(TPopulation), allocatable :: next_population(:)
    class(TIndiv), allocatable :: children(:)
    logical, allocatable :: mask(:)
    integer, allocatable :: parents_index(:)
    integer :: num_elite, quota, pop_size_t, p, c, i, j

    pop_size_t = this%pop_size

    allocate(next_population(pop_size_t))
    allocate(mask(pop_size_t))
    call this%preserve_elite(next_population, num_elite)
    i = num_elite
    c = 1

    do while (i < pop_size_t)
      mask = [(this%population(i)%indiv%feasible, i = 1, pop_size_t)]
      call this%selection%call(this%population%fitness, mask, parents_index)
      ! print *, parents_index
      ! stop
      ! parents_index = shuffle(this%pop_size, 2)
      call this%reproduce(parents_index, children)

      if (this%dup_rejection) then
        quota = 0
        do j = 1, size(children)
          p = i + quota
          if (p >= pop_size_t) exit
          if (p > 0 .and. .not. this%check_uniq(children(j), next_population(1:p))) cycle
          allocate(next_population(p+1)%indiv, source=children(j))
          quota = quota + 1
        end do
      else
        quota = min(size(children), pop_size_t - i)
        do j = 1, quota
          allocate(next_population(i+j)%indiv, source=children(j))
        end do
      end if

      i = i + quota
      if (quota == 0) c = c + 1
    end do

    call this%evaluate(next_population)
    call this%calc_fitness(next_population)

    ! call this%selection%call(next_population%fitness, parents_index)
    ! parents_index = reverse(sort(next_population%fitness))
    ! this%population = next_population(parents_index(1:this%pop_size))
    ! call this%calc_fitness(this%population)

    ! print *, size(parents_index)
    ! print "(10i6)", parents_index(sort(parents_index))
    ! print *, "uniq:", count([(any(i == parents_index), i = 1, pop_size_t)])
    ! stop
    ! return

    call move_alloc(from=next_population, to=this%population)
  end subroutine evolution
end module nsga2c_dp
