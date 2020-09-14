module nsga2c
    use util, only: elapsed_seconds
    use interface
    use individual
    use problem
    use ga_unit
    use soga
    use nsga2
    implicit none

    private
    public :: TNSGA2C

    type, extends(TNSGA2) :: TNSGA2C
        contains

        generic :: save_result => save_result_elite_feasible
        generic :: save_history => save_history_elite_feasible

        procedure :: set_prototype1
        procedure :: logger
        procedure :: save_result_elite_feasible, save_history_elite_feasible
        procedure :: calc_fitness
    end type TNSGA2C

    contains


    subroutine set_prototype1(this)
        implicit none
        class(TNSGA2C), intent(inout) :: this

        if (.not. allocated(this%prototype)) allocate(TIndivC::this%prototype)
    end subroutine set_prototype1


    ! ============================================================================
    ! calculation body
    ! ============================================================================

    subroutine calc_fitness(this, population)
        implicit none
        class(TNSGA2C), intent(inout) :: this
        type(TPopulation), intent(inout) :: population(:)
        real(8), allocatable :: objectives(:, :), constraints(:)
        logical, allocatable :: feasible(:)
        real(8) :: a = 0.1d0
        integer :: pop_size, i

        pop_size = size(population)
        objectives = reshape([(population(i)%indiv%objectives, i = 1, pop_size)], &
                             [this%num_obj, pop_size])
        constraints = [(scalar_constraints(population(i)%indiv%constraints), i = 1, pop_size)]
        feasible = [(population(i)%indiv%feasible, i = 1, pop_size)]
        population%rank = rank_pareto(objectives, constraints, feasible, dominatedp)

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

    logical function dominatedp(val, con, feasible) result(res)
        implicit none
        real(8), intent(in) :: val(:, :), con(:)
        logical, intent(in) :: feasible(:)

        res = (.not. feasible(1) .and. feasible(2)) .or.                        &
              (.not. (feasible(1) .or. feasible(2)) .and. con(1) > con(2)) .or. &
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
    ! IO
    ! ============================================================================

    subroutine logger(this, n, total)
        implicit none
        class(TNSGA2C), intent(in) :: this
        integer, intent(in) :: n, total
        integer, save :: last_time = 0
        real(8) :: d

        if (n == 0) then
            print "(2a)", "Start: ", this%problem_type

        else if (n == -1) then
            print "(2a/)", "End: ", this%problem_type

        ! else if (mod(n, 100) == 0) then
        else if (last_time == 0 .or. elapsed_seconds(last_time) > 10) then
            call system_clock(last_time)
            d = elapsed_seconds(this%start_time)

            print "(a,5(i0,a))",                                       &
                  "Progress: ", n, "steps, Feasible: ",                &
                  count_feasible(this%population), "/", this%pop_size, &
                  ", Elapsed: ", int(d), "s, Remain: ",                &
                  int(d * (total-n) / (n+1)), "s"
        end if
    end subroutine logger

    subroutine save_result_elite_feasible(this, filename, elite, feasible)
        implicit none
        class(TNSGA2C), intent(in) :: this
        character(*), intent(in) :: filename, elite, feasible
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
                    if (feasible == "only" .and. .not. p%indiv%feasible) cycle
                    if (feasible == "not" .and. p%indiv%feasible) cycle

                    call this%print_indiv(p%indiv, unit, .true.)
                    write(unit, "(',',i0,2(',',es15.8))") p%rank, p%fitness, p%crowding
                end associate
            end do
        close(unit)
    end subroutine save_result_elite_feasible

    subroutine save_history_elite_feasible(this, filename, elite, feasible)
        implicit none
        class(TNSGA2C), intent(in) :: this
        character(*), intent(in) :: filename, elite, feasible
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
                        if (feasible == "only" .and. .not. h%indiv%feasible) cycle
                        if (feasible == "not" .and. h%indiv%feasible) cycle

                        write(unit, "(i0,',')", advance='no') j - 1
                        call this%print_indiv(h%indiv, unit, .false.)

                        write(unit, "(',',i0,2(',',es15.8),2(',',i0))") &
                            h%rank, h%fitness, h%crowding, h%indiv%parents_id
                    end associate
                end do inner
            end do outer
        close(unit)
    end subroutine save_history_elite_feasible


    ! ============================================================================
    ! scalar function (temp)
    ! ============================================================================

    real(8) function scalar_constraints(constraints) result(res)
        implicit none
        real(8), intent(in) :: constraints(:)

        res = -sum(min(constraints, 0d0))
    end function scalar_constraints


    ! ============================================================================
    ! other
    ! ============================================================================

    integer function count_feasible(population) result(res)
        implicit none
        type(TPopulation), intent(in) :: population(:)
        integer :: i

        res = 0
        do i = 1, size(population)
            if (population(i)%indiv%feasible) res = res + 1
        end do
    end function count_feasible
end module nsga2c
