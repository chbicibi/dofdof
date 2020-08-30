module soga
  use interface
  use util
  use individual
  use problem
  use ga_unit
  use basic_optimizer
  implicit none

  private
  public :: TSOGA, TPopulation

  type :: TPopulation
    class(TIndiv), allocatable :: indiv
    logical :: init     = .false.
    integer :: rank     = -1
    real(8) :: fitness  = 0d0
    real(8) :: crowding = 0d0

    contains

    final :: destroy_population
  end type TPopulation

  type, extends(TOptimizer) :: TSOGA
    character(:), allocatable :: selection_type
    character(:), allocatable :: crossover_type
    character(:), allocatable :: mutation_type
    character(:), allocatable :: fitness_type

    logical :: elite_preservation = .false.
    logical :: dup_rejection = .false.
    logical :: sharing = .false.
    logical :: history_preservation = .true.

    real(8) :: fitness_weight = 0.1d0

    type(TPopulation), allocatable :: population(:)
    type(TPopulation), allocatable :: history(:, :)

    class(TSelection), allocatable :: selection
    class(TCrossover), allocatable :: crossover
    class(TMutation), allocatable :: mutation

    contains

    generic :: initialize => initialize_instance1, initialize_instance2

    generic :: run => run_default, run_with_hook
    generic :: evaluate => evaluate_pop
    generic :: save_result => save_result_default, save_result_elite
    generic :: save_history => save_history_default, save_history_elite
    generic :: best => best_indiv, best_dvar, best_ivar

    procedure :: initialize_instance1, initialize_instance2
    procedure :: set_parameter
    procedure :: set_operator
    procedure :: set_selection
    procedure :: set_crossover
    procedure :: set_mutation
    procedure :: set_fitness_type

    procedure :: prepare_calculation
    procedure :: init_population

    procedure :: run_default, run_with_hook
    procedure :: evolution
    procedure :: reproduce
    procedure :: evaluate_pop
    procedure :: calc_fitness
    procedure :: preserve_elite

    procedure :: save_result_default, save_result_elite
    procedure :: save_history_default, save_history_elite
    procedure :: alloc_history, keep_population
    procedure :: best_indiv, best_dvar, best_ivar, best_obj
    procedure :: check_uniq

    final :: destroy_instance
  end type TSOGA

  integer, parameter :: LIM_APPEND_CYCLE = 1000

  interface
    module subroutine initialize_instance1(this, nx, N, f, selection, crossover, mutation)
      class(TSOGA), intent(inout) :: this
      integer, intent(in) :: nx, N
      character(*), intent(in) :: f, selection, crossover, mutation
    end subroutine initialize_instance1

    module subroutine initialize_instance2(this, nx, N, selection, crossover, mutation)
      class(TSOGA), intent(inout) :: this
      integer, intent(in) :: nx, N
      character(*), intent(in) :: selection, crossover, mutation
    end subroutine initialize_instance2

    module subroutine set_parameter(this, nx, m, N, f, selection, crossover, mutation)
      class(TSOGA), intent(inout) :: this
      integer, intent(in) :: nx, m, N
      character(*), intent(in) :: f, selection, crossover, mutation
    end subroutine set_parameter

    module subroutine set_selection(this, selection_type, num_selection, num_tournament)
      class(TSOGA), intent(inout) :: this
      character(*), intent(in) :: selection_type
      integer, intent(in) :: num_selection, num_tournament
    end subroutine set_selection

    module subroutine set_crossover(this, crossover_type, rate, param)
      class(TSOGA), intent(inout) :: this
      character(*), intent(in) :: crossover_type
      real(8), intent(in) :: rate, param
    end subroutine set_crossover

    module subroutine set_mutation(this, mutation_type, rate, param)
      class(TSOGA), intent(inout) :: this
      character(*), intent(in) :: mutation_type
      real(8), intent(in) :: rate, param
    end subroutine set_mutation

    module subroutine set_fitness_type(this, fitness_type)
      class(TSOGA), intent(inout) :: this
      character(*), intent(in) :: fitness_type
    end subroutine set_fitness_type
  end interface

  contains


  ! ============================================================================
  ! setter
  ! ============================================================================

  subroutine set_operator(this) ! it doesn't work in submodule
    implicit none
    class(TSOGA), intent(inout) :: this

    allocate(TSelection::this%selection)
    allocate(TCrossover::this%crossover)
    allocate(TMutation::this%mutation)

    call this%selection%initialize(this%selection_type, 2, 2)
    call this%crossover%initialize(this%crossover_type, 1d0, 0.5d0)
    call this%mutation%initialize(this%mutation_type, 0.1d0, 20d0)
  end subroutine set_operator


  ! ============================================================================
  ! preprocessing
  ! ============================================================================

  subroutine prepare_calculation(this)
    implicit none
    class(TSOGA), intent(inout) :: this

    this%current_step = 0
    call system_clock(this%start_time)

    call this%init_population
    call this%evaluate(this%population)
    call this%calc_fitness(this%population)
  end subroutine prepare_calculation

  subroutine init_population(this)
    implicit none
    class(TSOGA), intent(inout) :: this
    integer :: i

    call this%set_prototype
    if (allocated(this%population)) deallocate(this%population)
    allocate(this%population(this%pop_size))

    do i = 1, this%pop_size
      allocate(this%population(i)%indiv, source=this%prototype)
      allocate(this%population(i)%indiv%parents_id, source=[-1, -1])
      call this%init_indiv(this%population(i)%indiv)
    end do
    this%population%init = .true.
  end subroutine init_population


  ! ============================================================================
  ! calculation body
  ! ============================================================================

  subroutine run_default(this, num_generation)
    implicit none
    class(TSOGA), intent(inout) :: this
    integer, intent(in) :: num_generation
    integer :: offset, i

    call this%alloc_history(num_generation, offset)
    call this%logger(0, num_generation)

    do i = 1, num_generation
      this%current_step = i
      call this%evolution
      call this%keep_population(i + offset)
      call this%logger(i, num_generation)
    end do

    call this%logger(-1, num_generation)
  end subroutine run_default

  subroutine run_with_hook(this, num_generation, proc)
    implicit none

    interface
      module subroutine proc(args)
        class(TSOGA), intent(inout) :: args
      end subroutine proc
    end interface

    class(TSOGA), intent(inout) :: this
    integer, intent(in) :: num_generation
    integer :: offset, i

    call this%alloc_history(num_generation, offset)
    call this%logger(0, num_generation)
    call proc(this)

    do i = 1, num_generation
      this%current_step = i
      call this%evolution
      call this%keep_population(i + offset)
      call this%logger(i, num_generation)
      call proc(this)
    end do

    call this%logger(-1, num_generation)
  end subroutine run_with_hook

  subroutine evolution(this)
    implicit none
    class(TSOGA), intent(inout) :: this
    type(TPopulation), allocatable :: next_population(:)
    class(TIndiv), allocatable :: children(:)
    integer, allocatable :: parents_index(:)
    integer :: num_elite, quota, pop_size_t, p, c, i, j

    pop_size_t = this%pop_size

    allocate(next_population(pop_size_t))
    call this%preserve_elite(next_population, num_elite)
    i = num_elite
    c = 1

    do while (i < pop_size_t)
      call this%selection%call(this%population%fitness, parents_index)
      ! print *, parents_index
      ! stop
      ! parents_index = shuffle(this%pop_size, 2)
      call this%reproduce(parents_index, children)

      if (this%dup_rejection .and. c <= LIM_APPEND_CYCLE) then
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

    next_population%init = .true.
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

  subroutine reproduce(this, index, children)
    implicit none
    class(TSOGA), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: children(:)
    real(8), allocatable :: parents_value(:, :), children_value(:, :)
    integer, allocatable :: parents_id(:)
    integer :: i

    parents_value = reshape([(this%population(index(i))%indiv%dvariables, i = 1, 2)], &
                            [this%num_var, 2])
    parents_id = [(this%population(index(i))%indiv%id, i = 1, 2)]

    call this%crossover%call(parents_value, children_value)
    allocate(children(size(children_value, dim=2)), source=this%prototype)

    do i = 1, size(children)
      allocate(children(i)%parents_id, source=parents_id)
      call children(i)%set_variables(children_value(:, i))
      call children(i)%clamp_variables(lower=0d0, upper=1d0)
      call this%mutation%call(children(i)%dvariables)
      call children(i)%clamp_variables(lower=0d0, upper=1d0)
    end do
  end subroutine reproduce

  subroutine evaluate_pop(this, population)
    implicit none
    class(TSOGA), intent(inout) :: this
    type(TPopulation), intent(inout) :: population(:)
    integer :: i

    do i = 1, size(population)
      print *, "evaluate:", i, "/", size(population)
      call this%evaluate(population(i)%indiv)
    end do
  end subroutine evaluate_pop

  subroutine calc_fitness(this, population) ! TODO
    implicit none
    class(TSOGA), intent(inout) :: this
    type(TPopulation), intent(inout) :: population(:)
    integer, allocatable :: index(:), order(:)
    real(8), allocatable :: objective(:)
    logical, allocatable :: feasible(:)
    real(8) :: a
    integer :: pop_size, i

    a = this%fitness_weight
    pop_size = size(population)

    allocate(objective(pop_size), source=[(population(i)%indiv%obj(), i = 1, pop_size)])
    allocate(feasible(pop_size), source=[(population(i)%indiv%feasible, i = 1, pop_size)])

    if (count(feasible) <= 1) then
      index = integers(pop_size)
    else
      index = pack(integers(pop_size), feasible)
    end if
    order = sort(objective, index)

    population(order)%rank = index

    select case (this%fitness_type)
    case ("VALUE")
      where (feasible)
        population%fitness = 1d0 / (1d0 + objective)
        ! population%fitness = (sum(objective) - objective) / sum(objective)
      else where
        population%fitness = 0d0
      end where
      ! do i = 1, 10
      !   print *, i, population(order(i))%rank, population(order(i))%indiv%obj(), population(order(i))%fitness
      ! end do
      ! do i = 5990, 6000
      !   print *, i, population(order(i))%rank, population(order(i))%indiv%obj(), population(order(i))%fitness
      ! end do
      ! stop
    case ("RANK")
      where (feasible)
        population%fitness = a * (1.0d0 - a) ** (population%rank - 1)
      else where
        population%fitness = 0d0
      end where
    case default
      write(0, "(3a/ a)") "Error: unknown fitness_type '", this%fitness_type, "'", &
                          "in subroutine 'calc_fitness'"
      call exit(1)
    end select
  end subroutine calc_fitness

  subroutine preserve_elite(this, new_population, num_elite)
    implicit none
    class(TSOGA), intent(in) :: this
    type(TPopulation), intent(inout) :: new_population(:)
    integer, intent(out) :: num_elite

    if (this%elite_preservation) then
      new_population(1) = this%population(maxloc(this%population%fitness, dim=1))
      num_elite = 1
    else
      num_elite = 0
    end if
  end subroutine preserve_elite


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine save_result_default(this, filename)
    implicit none
    class(TSOGA), intent(in) :: this
    character(*), intent(in) :: filename

    call this%save_result(filename, elite="all")
  end subroutine save_result_default

  subroutine save_result_elite(this, filename, elite)
    implicit none
    class(TSOGA), intent(in) :: this
    character(*), intent(in) :: filename, elite
    integer :: unit, i

    open(newunit=unit, file=filename)
      call this%population(1)%indiv%print_header_wv(unit)
      write(unit, "(a)") ",fitness"

      do i = 1, this%pop_size
        if (elite == "all" .or. xor(elite == "only", this%population(i)%rank > 1)) then
          call this%print_indiv(this%population(i)%indiv, unit, .true.)
          write(unit, "(','es15.8)") this%population(i)%fitness
        end if
      end do
    close(unit)
  end subroutine save_result_elite

  subroutine save_history_default(this, filename)
    implicit none
    class(TSOGA), intent(in) :: this
    character(*), intent(in) :: filename

    call this%save_history(filename, elite="all")
  end subroutine save_history_default

  subroutine save_history_elite(this, filename, elite)
    implicit none
    class(TSOGA), intent(in) :: this
    character(*), intent(in) :: filename, elite
    integer :: unit, index, i, j

    open(newunit=unit, file=filename)
      write(unit, "(a)", advance='no') "step,"
      call this%population(1)%indiv%print_header(unit)
      write(unit, "(a)") ",fitness"

      do j = 1, size(this%history, dim=2)
        if (.not. this%history(1, j)%init) exit

        if (elite == "best") then
          index = minloc(this%history(:, j)%rank, dim=1)
          write(unit, "(i0',')", advance='no') j - 1
          call this%print_indiv(this%history(index, j)%indiv, unit)
          write(unit, "(','es15.8)") this%history(index, j)%fitness
        else
          do i = 1, this%pop_size
            if (elite == "all" .or. xor(elite == "only", this%history(i, j)%rank > 1)) then
              write(unit, "(i0',')", advance='no') j - 1
              call this%print_indiv(this%history(i, j)%indiv, unit)
              write(unit, "(','es15.8)") this%history(i, j)%fitness
            end if
          end do
        end if
      end do
    close(unit)
  end subroutine save_history_elite


  ! ============================================================================
  ! history operation
  ! ============================================================================

  subroutine alloc_history(this, len, cols)
    implicit none
    class(TSOGA), intent(inout) :: this
    integer, intent(in) :: len
    integer, intent(out) :: cols
    type(TPopulation), allocatable :: temp(:, :)
    integer :: rows

    if(.not. this%history_preservation) then
      if (.not. allocated(this%history)) allocate(this%history(0, 0))
      return
    end if

    if (allocated(this%history)) then
      rows = size(this%history, dim=1)
      cols = size(this%history, dim=2)
      allocate(temp(rows, cols + len))
      temp(:, 1:cols) = this%history
      call move_alloc(from=temp, to=this%history)
    else
      allocate(this%history(this%pop_size, len + 1))
      this%history(:, 1) = this%population
      cols = 1
    end if
  end subroutine alloc_history

  subroutine keep_population(this, n)
    implicit none
    class(TSOGA), intent(inout) :: this
    integer, intent(in) :: n

    if (this%history_preservation) this%history(:, n) = this%population
  end subroutine keep_population


  ! ============================================================================
  ! other
  ! ============================================================================

  subroutine best_indiv(this, best)
    implicit none
    class(TSOGA), intent(in) :: this
    class(TIndiv), intent(out), allocatable :: best
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(), i = 1, this%pop_size)], dim=1)
    best = this%population(index)%indiv
  end subroutine best_indiv

  subroutine best_dvar(this, best)
    implicit none
    class(TSOGA), intent(in) :: this
    real(8), intent(out), allocatable :: best(:)
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(), i = 1, this%pop_size)], dim=1)
    best = this%dec(this%population(index)%indiv%dvariables)
  end subroutine best_dvar

  subroutine best_ivar(this, best)
    implicit none
    class(TSOGA), intent(in) :: this
    integer, intent(out), allocatable :: best(:)
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(), i = 1, this%pop_size)], dim=1)
    best = this%population(index)%indiv%ivariables
  end subroutine best_ivar

  subroutine best_obj(this, best)
    implicit none
    class(TSOGA), intent(in) :: this
    real(8), intent(out) :: best
    integer :: i

    best = minval([(this%population(i)%indiv%obj(), i = 1, this%pop_size)], dim=1)
  end subroutine best_obj

  logical function check_uniq(this, indiv, population) result(flag)
    implicit none
    class(TSOGA), intent(in) :: this
    class(TIndiv), intent(in) :: indiv
    type(TPopulation), intent(in) :: population(:)
    integer :: i

    do i = 1, size(population)
      if (indiv%is_same(population(i)%indiv)) then
        flag = .false.
        return
      end if
    end do
    flag = .true.
  end function check_uniq


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_instance(this)
    implicit none
    type(TSOGA), intent(inout) :: this

    if (allocated(this%population)) deallocate(this%population)
    if (allocated(this%history)) deallocate(this%history)
    if (allocated(this%selection)) deallocate(this%selection)
    if (allocated(this%crossover)) deallocate(this%crossover)
    if (allocated(this%mutation)) deallocate(this%mutation)
  end subroutine destroy_instance

  elemental subroutine destroy_population(this)
    implicit none
    type(TPopulation), intent(inout) :: this

    if (allocated(this%indiv)) deallocate(this%indiv)
  end subroutine destroy_population
end module soga
