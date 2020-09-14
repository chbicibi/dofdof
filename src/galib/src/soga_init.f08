submodule (soga) soga_init
    ! use ga_unit
    implicit none

    contains


    ! ============================================================================
    ! constructor
    ! ============================================================================

    subroutine initialize_instance1(this, nx, N, f, selection, crossover, mutation)
        implicit none
        class(TSOGA), intent(inout) :: this
        integer, intent(in) :: nx, N
        character(*), intent(in) :: f, selection, crossover, mutation

        call this%set_parameter(nx, 1, N, f, selection, crossover, mutation)
        call this%set_operator
        call this%set_problem
        ! call this%prepare_calculation
    end subroutine initialize_instance1

    subroutine initialize_instance2(this, nx, N, selection, crossover, mutation)
        implicit none
        class(TSOGA), intent(inout) :: this
        integer, intent(in) :: nx, N
        character(*), intent(in) :: selection, crossover, mutation

        call this%set_parameter(nx, 1, N, "USER", selection, crossover, mutation)
        call this%set_operator
        call this%set_problem
    end subroutine initialize_instance2


    ! ============================================================================
    ! setter
    ! ============================================================================

    subroutine set_parameter(this, nx, m, N, f, selection, crossover, mutation)
        implicit none
        class(TSOGA), intent(inout) :: this
        integer, intent(in) :: nx, m, N
        character(*), intent(in) :: f, selection, crossover, mutation

        this%num_var = nx
        this%num_obj = m
        this%pop_size = N
        this%problem_type = trim(f)
        this%selection_type = trim(selection)
        this%crossover_type = trim(crossover)
        this%mutation_type = trim(mutation)
        this%fitness_type = "VALUE"
    end subroutine set_parameter

    subroutine set_selection(this, selection_type, num_selection, num_tournament)
        implicit none
        class(TSOGA), intent(inout) :: this
        character(*), intent(in) :: selection_type
        integer, intent(in) :: num_selection, num_tournament

        this%selection_type = selection_type
        call this%selection%initialize(selection_type, num_selection, num_tournament)
    end subroutine set_selection

    subroutine set_crossover(this, crossover_type, rate, param)
        implicit none
        class(TSOGA), intent(inout) :: this
        character(*), intent(in) :: crossover_type
        real(8), intent(in) :: rate, param

        this%crossover_type = crossover_type
        call this%crossover%initialize(crossover_type, rate, param)
    end subroutine set_crossover

    subroutine set_mutation(this, mutation_type, rate, param)
        implicit none
        class(TSOGA), intent(inout) :: this
        character(*), intent(in) :: mutation_type
        real(8), intent(in) :: rate, param

        this%mutation_type = mutation_type
        call this%mutation%initialize(mutation_type, rate, param)
    end subroutine set_mutation

    subroutine set_fitness_type(this, fitness_type)
        implicit none
        class(TSOGA), intent(inout) :: this
        character(*), intent(in) :: fitness_type

        this%fitness_type = fitness_type
    end subroutine set_fitness_type
end submodule soga_init
