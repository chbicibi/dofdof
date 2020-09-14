submodule (nsga2) nsga2_init
    implicit none

    contains


    ! ============================================================================
    ! constructor
    ! ============================================================================

    subroutine initialize_instance3(this, nx, m, N, f, selection, crossover, mutation)
        implicit none
        class(TNSGA2), intent(inout) :: this
        integer, intent(in) :: nx, m, N
        character(*), intent(in) :: f, selection, crossover, mutation

        call this%set_parameter(nx, m, N, f, selection, crossover, mutation)
        call this%set_operator
        call this%set_problem
        ! call this%prepare_calculation
    end subroutine initialize_instance3

    subroutine initialize_instance4(this, nx, m, N, selection, crossover, mutation)
        implicit none
        class(TNSGA2), intent(inout) :: this
        integer, intent(in) :: nx, m, N
        character(*), intent(in) :: selection, crossover, mutation

        call this%set_parameter(nx, m, N, "USER", selection, crossover, mutation)
        call this%set_operator
        call this%set_problem
        ! call this%prepare_calculation
    end subroutine initialize_instance4
end submodule nsga2_init
