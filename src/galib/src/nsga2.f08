module nsga2
  use interface
  use util
  use individual
  use problem
  use ga_unit
  use soga
  implicit none

  private
  public :: TNSGA2

  type, extends(TSOGA) :: TNSGA2
    contains

    generic :: initialize => initialize_instance3, initialize_instance4
    generic :: best => best_indiv_n, best_dvar_n, best_ivar_n

    procedure :: initialize_instance3, initialize_instance4
    procedure :: set_prototype1
    procedure :: calc_fitness
    procedure :: preserve_elite
    procedure :: save_result_elite
    procedure :: save_history_elite
    procedure :: best_indiv_n, best_dvar_n, best_ivar_n
  end type TNSGA2

  interface
    module subroutine initialize_instance3(this, nx, m, N, f, selection, crossover, mutation)
      class(TNSGA2), intent(inout) :: this
      integer, intent(in) :: nx, m, N
      character(*), intent(in) :: f, selection, crossover, mutation
    end subroutine initialize_instance3

    module subroutine initialize_instance4(this, nx, m, N, selection, crossover, mutation)
      class(TNSGA2), intent(inout) :: this
      integer, intent(in) :: nx, m, N
      character(*), intent(in) :: selection, crossover, mutation
    end subroutine initialize_instance4
  end interface

  contains


  subroutine set_prototype1(this)
    implicit none
    class(TNSGA2), intent(inout) :: this

    if (.not. allocated(this%prototype)) allocate(TIndiv::this%prototype)
  end subroutine set_prototype1


  ! ============================================================================
  ! calculation body
  ! ============================================================================

  subroutine calc_fitness(this, population)
    implicit none
    class(TNSGA2), intent(inout) :: this
    type(TPopulation), intent(inout) :: population(:)
    real(8), allocatable :: objectives(:, :)
    real(8) :: a
    integer :: pop_size, i

    a = this%fitness_weight
    pop_size = size(population)

    objectives = reshape([(population(i)%indiv%objectives(1:this%num_obj), i = 1, pop_size)], &
                         [this%num_obj, pop_size])

    population%rank = rank_pareto(objectives)

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

  subroutine preserve_elite(this, new_population, num_elite)
    implicit none
    class(TNSGA2), intent(in) :: this
    type(TPopulation), intent(inout) :: new_population(:)
    integer, intent(out) :: num_elite
    integer :: i

    if (this%elite_preservation) then
      ! num_elite = count(this%population%rank == 1)
      num_elite = 0
      do i = 1, this%pop_size
        if (this%population(i)%rank /= 1) cycle
        if (num_elite > 0 .and. .not. this%check_uniq(this%population(i)%indiv, new_population(1:num_elite))) cycle
        num_elite = num_elite + 1
        new_population(num_elite) = this%population(i)
      end do
      ! new_population(1:num_elite) = pack(this%population, this%population%rank == 1)
    else
      num_elite = 0
    end if
  end subroutine preserve_elite


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine save_result_elite(this, filename, elite)
    implicit none
    class(TNSGA2), intent(in) :: this
    character(*), intent(in) :: filename, elite
    integer :: unit, i

    open(newunit=unit, file=filename)
      call this%population(1)%indiv%print_header_wv(unit)
      write(unit, "(a)") ",rank,fitness,crowding-dist"

      do i = 1, this%pop_size
        if (elite == "all" .or. xor(elite == "only", this%population(i)%rank > 1)) then
          call this%print_indiv(this%population(i)%indiv, unit, .true.)
          write(unit, "(','i0,2(','es15.8))") this%population(i)%rank,    &
                                              this%population(i)%fitness, &
                                              this%population(i)%crowding
        end if
      end do
    close(unit)
  end subroutine save_result_elite

  subroutine save_history_elite(this, filename, elite)
    implicit none
    class(TNSGA2), intent(in) :: this
    character(*), intent(in) :: filename, elite
    integer :: unit, i, j

    open(newunit=unit, file=filename)
      write(unit, "(a)", advance='no') "step,"
      call this%population(1)%indiv%print_header(unit)
      write(unit, "(a)") ",rank,fitness,crowding-dist,pid1,pid2"

      outer: do j = 1, size(this%history, dim=2)
        inner: do i = 1, this%pop_size
          if (.not. this%history(i, j)%init) exit outer

          if (elite == "all" .or. xor(elite == "only", this%history(i, j)%rank > 1)) then
            write(unit, "(i0',')", advance='no') j - 1
            call this%print_indiv(this%history(i, j)%indiv, unit)

            write(unit, "(','i0,2(','es15.8),2(','i0))") &
              this%history(i, j)%rank,      &
              this%history(i, j)%fitness,   &
              this%history(i, j)%crowding,  &
              this%history(i, j)%indiv%parents_id
          end if
        end do inner
      end do outer
    close(unit)
  end subroutine save_history_elite


  ! ============================================================================
  ! other
  ! ============================================================================

  subroutine best_indiv_n(this, n, best)
    implicit none
    class(TNSGA2), intent(in) :: this
    integer, intent(in) :: n
    class(TIndiv), intent(out), allocatable :: best
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(n), i = 1, this%pop_size)], dim=1)
    best = this%population(index)%indiv
  end subroutine best_indiv_n


  subroutine best_dvar_n(this, n, best)
    implicit none
    class(TNSGA2), intent(in) :: this
    integer, intent(in) :: n
    real(8), intent(out), allocatable :: best(:)
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(n), i = 1, this%pop_size)], dim=1)
    best = this%dec(this%population(index)%indiv%dvariables)
  end subroutine best_dvar_n

  subroutine best_ivar_n(this, n, best)
    implicit none
    class(TNSGA2), intent(in) :: this
    integer, intent(in) :: n
    integer, intent(out), allocatable :: best(:)
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(n), i = 1, this%pop_size)], dim=1)
    best = this%population(index)%indiv%ivariables
  end subroutine best_ivar_n
end module nsga2
