module moead
  use interface
  use util
  use individual
  use problem
  use moead_unit
  use ga_unit
  use basic_optimizer
  implicit none

  private
  public :: TMOEAD, TPopulation

  type :: TPopulation
    class(TIndiv), allocatable :: indiv
    real(8), allocatable :: weight(:)
    integer, allocatable :: table(:)
    real(8), allocatable :: fitness(:)

    contains

    final :: destroy_population
  end type TPopulation

  type, extends(TOptimizer) :: TMOEAD
    integer :: table_size

    character(:), allocatable :: selection_type
    character(:), allocatable :: crossover_type
    character(:), allocatable :: mutation_type
    character(:), allocatable :: fitness_type
    character(:), allocatable :: scalar_func_type

    real(8), allocatable :: reference_point(:)
    type(TPopulation), allocatable :: population(:)
    type(TPopulation), allocatable :: history(:, :)

    logical :: history_preservation = .true.

    class(TSelection), allocatable :: selection
    class(TCrossover), allocatable :: crossover
    class(TMutation), allocatable :: mutation

    procedure(func_1d_1d_1d_d), nopass, pointer :: scalar_func

    contains

    generic :: initialize => initialize_instance1, initialize_instance2
    generic :: save_result => save_result_base
    generic :: save_history => save_history_base
    generic :: best => best_dvar, best_ivar, best_var_n

    procedure :: save_result_base => save_result
    procedure :: save_history_base => save_history
    procedure :: dump => dump_instance
    procedure :: load => load_instance

    procedure :: initialize_instance1, initialize_instance2
    procedure :: set_parameter
    procedure :: set_operator

    procedure :: prepare_calculation
    procedure :: init_weight
    procedure :: init_neighbourhood
    procedure :: init_population

    procedure :: update_reference
    procedure :: update_population
    procedure :: run
    procedure :: evolution
    procedure :: reproduce
    procedure :: calc_fitness
    procedure :: scalar_value

    procedure :: alloc_history, keep_population
    procedure :: best_dvar, best_ivar, best_var_n, best_obj

    final :: destroy_instance
  end type TMOEAD

  interface
    module subroutine initialize_instance1(this, nx, m, N, T, f, g, selection, crossover, mutation)
      class(TMOEAD), intent(inout) :: this
      integer, intent(in) :: nx, m, N, T
      character(*), intent(in) :: f, g, selection, crossover, mutation
    end subroutine initialize_instance1

    module subroutine initialize_instance2(this, nx, m, N, T, g, selection, crossover, mutation)
      class(TMOEAD), intent(inout) :: this
      integer, intent(in) :: nx, m, N, T
      character(*), intent(in) :: g, selection, crossover, mutation
    end subroutine initialize_instance2

    module subroutine set_parameter(this, nx, m, N, T, f, g, selection, crossover, mutation)
      class(TMOEAD), intent(inout) :: this
      integer, intent(in) :: nx, m, N, T
      character(*), intent(in) :: f, g, selection, crossover, mutation
    end subroutine set_parameter

    module subroutine set_selection(this, selection_type, num_selection, num_tournament)
      class(TMOEAD), intent(inout) :: this
      character(*), intent(in) :: selection_type
      integer, intent(in) :: num_selection, num_tournament
    end subroutine set_selection

    module subroutine set_crossover(this, crossover_type, rate, param)
      class(TMOEAD), intent(inout) :: this
      character(*), intent(in) :: crossover_type
      real(8), intent(in) :: rate, param
    end subroutine set_crossover

    module subroutine set_mutation(this, mutation_type, rate, param)
      class(TMOEAD), intent(inout) :: this
      character(*), intent(in) :: mutation_type
      real(8), intent(in) :: rate, param
    end subroutine set_mutation

    module subroutine set_fitness_type(this, fitness_type)
      class(TMOEAD), intent(inout) :: this
      character(*), intent(in) :: fitness_type
    end subroutine set_fitness_type
  end interface

  contains


  ! ============================================================================
  ! setter
  ! ============================================================================

  subroutine set_operator(this)
    implicit none
    class(TMOEAD), intent(inout) :: this

    allocate(TSelection::this%selection)
    allocate(TCrossover::this%crossover)
    allocate(TMutation::this%mutation)

    this%scalar_func => get_scalar_func(this%scalar_func_type)

    call this%selection%initialize(this%selection_type, 2, 2)
    call this%crossover%initialize(this%crossover_type, 1d0, 0.5d0)
    call this%mutation%initialize(this%mutation_type, 0.1d0, 20d0)
  end subroutine set_operator


  ! ============================================================================
  ! preprocessing
  ! ============================================================================

  subroutine prepare_calculation(this)
    implicit none
    class(TMOEAD), intent(inout) :: this

    this%reference_point = filled(this%num_obj, huge(0d0))

    if (allocated(this%population)) deallocate(this%population)
    allocate(this%population(this%pop_size))

    call this%init_weight        ! setting weight vectors (== lambda)
    call this%init_neighbourhood ! detecting closest weight vectors to each weight vector
    call this%init_population    ! setting variables and evaluating
  end subroutine prepare_calculation


  ! ******************************************************************
  ! initializing weight vectors
  ! ******************************************************************

  subroutine init_weight(this)
    implicit none
    class(TMOEAD), intent(inout) :: this

    select case (this%num_obj)
    case (2)
      call init_weight_2d(this)
    case (3)
      call init_weight_3d(this)
    case default
      write(0, "(ai0a/a)") "Error: wrong number of objectives (given ", &
                           this%num_obj, ", expected 2 or 3)",          &
                           "in subroutine 'init_weight'"
      call exit(1)
    end select
  end subroutine init_weight

  subroutine init_weight_2d(this)
    implicit none
    class(TMOEAD), intent(inout) :: this
    integer :: i

    do i = 1, this%pop_size
      this%population(i)%weight = [i, this%pop_size - i + 1]
    end do
  end subroutine init_weight_2d

  subroutine init_weight_3d(this)
    implicit none
    class(TMOEAD), intent(inout) :: this
    integer :: n, i, j, k

    n = ceiling(sqrt(2d0 * this%pop_size + 0.25d0) - 0.5d0)
    k = 0

    do i = 1, n
      do j = 1, min(i, this%pop_size - k)
        this%population(k + j)%weight = [n - i + 1, j, i - j + 1]
      end do
      k = k + i
    end do
  end subroutine init_weight_3d


  ! ******************************************************************
  ! detecting T closest vectors (T == this%table_size)
  ! ******************************************************************

  subroutine init_neighbourhood(this)
    implicit none
    class(TMOEAD), intent(inout) :: this
    real(8), allocatable :: dist(:)
    integer :: i, j

    do i = 1, this%pop_size
      dist = [(distance(this%population(i)%weight,  &
                        this%population(j)%weight), &
               j = 1, this%pop_size)]
      this%population(i)%table = sort(dist, this%table_size)
    end do
  end subroutine init_neighbourhood


  ! ******************************************************************
  ! initializing population randomly
  ! ******************************************************************

  subroutine init_population(this)
    implicit none
    class(TMOEAD), intent(inout) :: this
    integer :: i

    call this%set_prototype
    ! if (allocated(this%population)) deallocate(this%population)
    ! allocate(this%population(this%pop_size))

    do i = 1, this%pop_size
      allocate(this%population(i)%indiv, source=this%prototype)
      call this%init_indiv(this%population(i)%indiv)
      call this%evaluate(this%population(i)%indiv)
      call this%update_reference(this%population(i)%indiv)
    end do
    call this%calc_fitness(this%population)
  end subroutine init_population


  ! ============================================================================
  ! calculation body
  ! ============================================================================

  subroutine update_reference(this, indiv)
    implicit none
    class(TMOEAD), intent(inout) :: this
    class(TIndiv), intent(in) :: indiv

    this%reference_point = min(this%reference_point, indiv%objectives(1:this%num_obj))
  end subroutine update_reference

  subroutine update_population(this, indiv, table, idx)
    implicit none
    class(TMOEAD), intent(inout) :: this
    class(TIndiv), intent(in) :: indiv
    integer, intent(in) :: table(:), idx
    real(8) :: f_old, f_new
    integer :: i, j

    do i = 1, size(table)
      j = table(i)
      if (j /= idx) cycle
      f_old = this%scalar_value(this%population(j)%indiv, j)
      f_new = this%scalar_value(indiv, j)
      if (f_new < f_old) this%population(j)%indiv = indiv
    end do
  end subroutine update_population


  ! ============================================================================
  ! ============================================================================

  subroutine run(this, num_generation)
    implicit none
    class(TMOEAD), intent(inout) :: this
    integer, intent(in) :: num_generation
    integer :: offset, i

    call this%alloc_history(num_generation, offset)
    call this%logger(0, num_generation)

    do i = 1, num_generation
      call this%evolution
      call this%keep_population(i + offset)
      call this%logger(i, num_generation)
    end do

    call this%logger(-1, num_generation)
  end subroutine run

  subroutine evolution(this)
    implicit none
    class(TMOEAD), intent(inout) :: this
    class(TIndiv), allocatable :: child
    integer, allocatable :: parents_index(:)
    integer :: i

    do i = 1, this%pop_size
      associate (table=>this%population(i)%table)
        ! parents_index = table(shuffle(this%table_size, 2))
        call this%selection%call(this%population(i)%fitness, parents_index)
        call this%reproduce(table(parents_index), child)
        call this%evaluate(child)
        call this%update_reference(child)
        call this%update_population(child, table, i)
      end associate
    end do
    call this%calc_fitness(this%population)
  end subroutine evolution

  subroutine reproduce(this, index, child)
    implicit none
    class(TMOEAD), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: child
    real(8), allocatable :: parents_value(:, :)
    integer :: i

    allocate(child, source=this%prototype)

    parents_value = reshape([(this%population(index(i))%indiv%dvariables, i = 1, 2)], &
                            [this%num_var, 2])

    call this%crossover%call(parents_value, child%dvariables)
    call child%clamp_variables(lower=0d0, upper=1d0)
    call this%mutation%call(child%dvariables)
    call child%clamp_variables(lower=0d0, upper=1d0)
  end subroutine reproduce

  subroutine calc_fitness(this, population)
    implicit none
    class(TMOEAD), intent(inout) :: this
    type(TPopulation), intent(inout) :: population(:)
    real(8), allocatable :: objective(:)
    integer :: pop_size, i, j

    pop_size = size(population)
    allocate(objective(this%table_size))
    do j = 1, pop_size
      associate (table=>population(j)%table)
        objective = [(this%scalar_value(population(table(i))%indiv, j), i = 1, this%table_size)]
        population(j)%fitness = 1d0 / (1d0 + objective)
      end associate
    end do
  end subroutine calc_fitness

  real(8) function scalar_value(this, indiv, idx) result(res)
    implicit none
    class(TMOEAD), intent(in) :: this
    class(TIndiv), intent(in) :: indiv
    integer, intent(in) :: idx

    res = this%scalar_func(indiv%objectives(1:this%num_obj), &
                           this%population(idx)%weight,      &
                           this%reference_point)
  end function scalar_value


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine save_result(this, filename)
    implicit none
    class(TMOEAD), intent(in) :: this
    character(*), intent(in) :: filename
    integer :: unit, i

    open(newunit=unit, file=filename)
      call this%population(1)%indiv%print_header_wv(unit)
      write(unit, *)

      do i = 1, this%pop_size
        call this%print_indiv(this%population(i)%indiv, unit, .true.)
        write(unit, *)
      end do
    close(unit)
  end subroutine save_result

  subroutine save_history(this, filename)
    implicit none
    class(TMOEAD), intent(in) :: this
    character(*), intent(in) :: filename
    integer :: unit, i, j

    open(newunit=unit, file=filename)
      write(unit, "(a)", advance='no') "step,"
      call this%population(1)%indiv%print_header(unit)
      write(unit, *)

      do j = 1, size(this%history, dim=2)
        do i = 1, this%pop_size
          write(unit, "(i0',')", advance='no') j - 1
          call this%print_indiv(this%history(i, j)%indiv, unit)
          write(unit, *)
        end do
      end do
    close(unit)
  end subroutine save_history

  subroutine dump_instance(this, filename)
    implicit none
    class(TMOEAD), intent(in) :: this
    character(*), intent(in) :: filename
    integer :: unit, i

    open(newunit=unit, file=filename)
      write(unit, "(i0)") this%num_var
      write(unit, "(i0)") this%num_obj
      write(unit, "(i0)") this%pop_size
      write(unit, "(i0)") this%table_size
      write(unit, "(a)") this%problem_type
      write(unit, "(a)") this%scalar_func_type
      write(unit, "(a)") this%crossover_type
      write(unit, "(a)") this%mutation_type
      do i = 1, this%pop_size
        write(unit, *) this%population(i)%indiv%dvariables, &
                       this%population(i)%indiv%objectives
      end do
    close(unit)
  end subroutine dump_instance

  subroutine load_instance(this, filename)
    implicit none
    class(TMOEAD), intent(inout) :: this
    character(*), intent(in) :: filename
    integer :: nx, m, N, T
    character(64) :: f, g, crossover, mutation
    integer :: unit, i

    open(newunit=unit, file=filename)
      read(unit, *) nx
      read(unit, *) m
      read(unit, *) N
      read(unit, *) T
      read(unit, *) f
      read(unit, *) g
      read(unit, *) crossover
      read(unit, *) mutation
      call this%initialize(nx, m, N, T, f, g, crossover, mutation)
      do i = 1, this%pop_size
        read(unit, *) this%population(i)%indiv%dvariables, &
                      this%population(i)%indiv%objectives
      end do
    close(unit)
  end subroutine load_instance


  ! ============================================================================
  ! ============================================================================

  subroutine alloc_history(this, len, cols)
    implicit none
    class(TMOEAD), intent(inout) :: this
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
    class(TMOEAD), intent(inout) :: this
    integer, intent(in) :: n

    if (this%history_preservation) this%history(:, n) = this%population
  end subroutine keep_population


  ! ============================================================================
  ! other
  ! ============================================================================

  subroutine best_dvar(this, best)
    implicit none
    class(TMOEAD), intent(in) :: this
    real(8), intent(out), allocatable :: best(:)
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(), i = 1, this%pop_size)], dim=1)
    best = this%dec(this%population(index)%indiv%dvariables)
  end subroutine best_dvar

  subroutine best_ivar(this, best)
    implicit none
    class(TMOEAD), intent(in) :: this
    integer, intent(out), allocatable :: best(:)
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(), i = 1, this%pop_size)], dim=1)
    best = this%population(index)%indiv%ivariables
  end subroutine best_ivar

  subroutine best_var_n(this, n, best)
    implicit none
    class(TMOEAD), intent(in) :: this
    integer, intent(in) :: n
    real(8), intent(out), allocatable :: best(:)
    integer :: index, i

    index = minloc([(this%population(i)%indiv%obj(n), i = 1, this%pop_size)], dim=1)
    best = this%dec(this%population(index)%indiv%dvariables)
  end subroutine best_var_n

  subroutine best_obj(this, best)
    implicit none
    class(TMOEAD), intent(in) :: this
    real(8), intent(out) :: best
    integer :: i

    best = minval([(this%population(i)%indiv%obj(), i = 1, this%pop_size)], dim=1)
  end subroutine best_obj


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_instance(this)
    implicit none
    type(TMOEAD), intent(inout) :: this

    if (allocated(this%reference_point)) deallocate(this%reference_point)
    if (allocated(this%population)) deallocate(this%population)
    if (allocated(this%crossover)) deallocate(this%crossover)
    if (allocated(this%mutation)) deallocate(this%mutation)
    if (associated(this%scalar_func)) nullify(this%scalar_func)
  end subroutine destroy_instance

  elemental subroutine destroy_population(this)
    implicit none
    type(TPopulation), intent(inout) :: this

    if (allocated(this%indiv)) deallocate(this%indiv)
    if (allocated(this%weight)) deallocate(this%weight)
    if (allocated(this%table)) deallocate(this%table)
  end subroutine destroy_population
end module moead
