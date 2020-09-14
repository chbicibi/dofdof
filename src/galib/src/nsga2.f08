module nsga2
    use util, only: reverse, sort
    use interface
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
        procedure :: advance
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

    subroutine advance(this, next_population)
        implicit none
        class(TNSGA2), intent(inout) :: this
        type(TPopulation), intent(inout), allocatable :: next_population(:)
        type(TPopulation), allocatable :: temp_population(:)
        integer, allocatable :: rank_index(:)
        integer :: popsize, popsize_total

        ! call move_alloc(from=next_population, to=this%population)
        popsize = size(this%population)

        popsize_total = popsize + size(next_population)

        allocate(temp_population(popsize_total), source=[this%population, next_population])
        call this%calc_fitness(temp_population)
        rank_index = reverse(sort(temp_population%fitness))
        ! print *, rank_index(1:10)

        this%population = temp_population(rank_index(1:popsize))
        call this%calc_fitness(this%population)
        deallocate(next_population)
        deallocate(temp_population)
    end subroutine advance

    subroutine calc_fitness(this, population)
        implicit none
        class(TNSGA2), intent(inout) :: this
        type(TPopulation), intent(inout) :: population(:)
        real(8), allocatable :: objectives(:, :)
        real(8) :: a
        integer :: pop_size, i

        a = this%fitness_weight
        pop_size = size(population)

        objectives = reshape([(population(i)%indiv%objectives(1:this%num_obj), &
                               i = 1, pop_size)],                              &
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

        num_elite = 0
        return

        if (this%elite_preservation) then
            ! num_elite = count(this%population%rank == 1)
            num_elite = 0
            do i = 1, this%pop_size
                if (this%population(i)%rank /= 1) cycle
                if (num_elite > 0 .and. &
                    .not. this%check_uniq(this%population(i)%indiv, &
                                          new_population(1:num_elite))) cycle

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

        if (.not. this%population(1)%init) return

        open(newunit=unit, file=filename)
            ! print csv header
            call this%print_header(this%population(1)%indiv, unit, .true.)
            write(unit, "(a)") ",rank,fitness,crowding-dist"

            ! print csv body
            do i = 1, this%pop_size
                associate (p => this%population(i))
                    if (elite == "only" .and. p%rank > 1) cycle
                    if (elite == "not" .and. p%rank == 1) cycle

                    call this%print_indiv(p%indiv, unit, .true.)
                    write(unit, "(',',i0,2(',',es15.8))") p%rank, p%fitness, p%crowding
                end associate
            end do
        close(unit)
    end subroutine save_result_elite

    subroutine save_history_elite(this, filename, elite)
        implicit none
        class(TNSGA2), intent(in) :: this
        character(*), intent(in) :: filename, elite
        integer :: unit, i, j

        if (.not. this%history(1, 1)%init) return

        open(newunit=unit, file=filename)
            ! print csv header
            write(unit, "(a)", advance='no') "step,"
            call this%print_header(this%history(1, 1)%indiv, unit, .false.)
            write(unit, "(a)") ",rank,fitness,crowding-dist,pid1,pid2"

            ! print csv body
            outer: do j = 1, size(this%history, dim=2)
                inner: do i = 1, this%pop_size
                    associate (h => this%history(i, j))
                        if (.not. h%init) exit outer

                        if (elite == "only" .and. h%rank > 1) cycle
                        if (elite == "not" .and. h%rank == 1) cycle

                        write(unit, "(i0,',')", advance='no') j - 1
                        call this%print_indiv(h%indiv, unit, .false.)

                        write(unit, "(',',i0,2(',',es15.8),2(',',i0))") &
                            h%rank, h%fitness, h%crowding, h%indiv%parents_id
                    end associate
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

        index = minloc([(this%population(i)%indiv%obj(n), &
                         i = 1, this%pop_size)], dim=1)
        best = this%population(index)%indiv
    end subroutine best_indiv_n


    subroutine best_dvar_n(this, n, best)
        implicit none
        class(TNSGA2), intent(in) :: this
        integer, intent(in) :: n
        real(8), intent(out), allocatable :: best(:)
        integer :: index, i

        index = minloc([(this%population(i)%indiv%obj(n), &
                         i = 1, this%pop_size)], dim=1)
        best = this%dec(this%population(index)%indiv%dvariables)
    end subroutine best_dvar_n

    subroutine best_ivar_n(this, n, best)
        implicit none
        class(TNSGA2), intent(in) :: this
        integer, intent(in) :: n
        integer, intent(out), allocatable :: best(:)
        integer :: index, i

        index = minloc([(this%population(i)%indiv%obj(n), &
                         i = 1, this%pop_size)], dim=1)
        best = this%population(index)%indiv%ivariables
    end subroutine best_ivar_n
end module nsga2
