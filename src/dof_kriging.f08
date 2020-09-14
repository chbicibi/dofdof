module mod_kriging
    use util
    use individual
    use kriging
    use soga
    implicit none

    private
    public :: init_krig, estimate

    type(TKriging), allocatable :: model(:)


    contains

    ! 2019.6.17追加
    subroutine init_krig(kriging_file)
        character(*), intent(in) :: kriging_file
        character(:), allocatable :: input_files(:), theta_files(:)
        integer, allocatable :: data_cols(:, :)
        real(8), allocatable :: data_range(:, :)
        integer :: n_model_total, n_model, i_model, n_var, is_init
        integer :: unit, i

        open(newunit=unit, file=kriging_file, status="old")
            read(unit, *) n_model_total
            read(unit, *) n_model
            read(unit, *) i_model
            read(unit, *) n_var
            read(unit, *) is_init

            if (i_model == 1) then
                if (allocated(model)) deallocate(model)
                allocate(model(n_model_total))
            end if

            allocate(character(100) :: input_files(n_model))
            allocate(character(100) :: theta_files(n_model))
            allocate(data_cols(2, n_model))
            allocate(data_range(2, n_var))

            do i = 1, n_model
                read(unit, *) input_files(i)
                read(unit, *) theta_files(i)
                read(unit, *) data_cols(:, i)
            end do

            do i = 1, n_var
                read(unit, *) data_range(:, i)
            end do

            if (is_init == 1) then
                call make_krig(n_var, input_files, theta_files, data_cols, data_range, unit, model(i_model:i_model+n_model-1))
            else
                call load_krig(n_var, input_files, theta_files, data_cols, data_range, model(i_model:i_model+n_model-1))
            end if
        close(unit)
    end subroutine init_krig


    subroutine make_krig(n_var, input_files, theta_files, data_cols, data_range, unit, model_local)
        integer, intent(in) :: n_var, unit
        character(:), intent(in), allocatable :: input_files(:), theta_files(:)
        integer, intent(in) :: data_cols(:, :)
        real(8), intent(in) :: data_range(:, :)
        type(TKriging), intent(inout) :: model_local(:)

        type(TSOGA) :: optimizer
        real(8), allocatable :: variables(:, :), objectives(:), theta(:), scale(:)
        integer, allocatable :: idum1(:), idum2(:)
        integer :: n_sample, pop_size, steps
        integer :: unit1, i, j

        allocate(scale(n_var))
        read(unit, *) pop_size, steps
        read(unit, *) scale

        call optimizer%initialize(nx=n_var, N=pop_size, selection="T", crossover="BLX", mutation="PM")
        optimizer%elite_preservation = .true.
        optimizer%sharing = .true.

        do i = 1, size(model_local)
            open(newunit=unit1, file=input_files(i), status="old")
                read(unit1, *) n_sample

                allocate(variables(n_var, n_sample))
                allocate(objectives(n_sample))
                allocate(theta(n_var))
                allocate(idum1(data_cols(1, i) - 1))
                allocate(idum2(data_cols(2, i) - 1))

                do j = 1, n_sample
                    read(unit1, *) idum1, variables(:, j), idum2, objectives(j)
                end do

                deallocate(idum1, idum2)
            close(unit1)

            call model_local(i)%initialize(variables, objectives)
            call model_local(i)%set_bounds(data_range(1, :), data_range(2, :))
            call model_local(i)%set_scale(scale)

            call optimizer%set_problem(model_local(i))
            call optimizer%prepare_calculation
            call optimizer%run(steps)

            ! index = minloc([(optimizer%population(i)%indiv%o1(), i = 1, pop_size)], dim=1)
            ! theta = optimizer%population(index)%indiv%variables
            call optimizer%best(theta)

            open(newunit=unit1, file=theta_files(i), status="replace")
                write(unit1, *) model_local(i)%scale
                write(unit1, *) theta
            close(unit1)

            call model_local(i)%calc_likelihood(theta)
            deallocate(variables, objectives, theta)
        end do
    end subroutine make_krig


    subroutine load_krig(n_var, input_files, theta_files, data_cols, data_range, model_local)
        integer, intent(in) :: n_var
        character(:), intent(in), allocatable :: input_files(:), theta_files(:)
        integer, intent(in) :: data_cols(:, :)
        real(8), intent(in) :: data_range(:, :)
        type(TKriging), intent(inout) :: model_local(:)

        real(8), allocatable :: variables(:, :), objectives(:), theta(:), scale(:)
        integer, allocatable :: idum1(:), idum2(:)
        integer :: n_sample
        integer :: unit, i, j

        do i = 1, size(model_local)
            ! print *, i, input_files(i)

            open(newunit=unit, file=input_files(i), status="old")
                read(unit, *) n_sample

                allocate(variables(n_var, n_sample))
                allocate(objectives(n_sample))
                allocate(theta(n_var))
                allocate(scale(n_var))

                allocate(idum1(data_cols(1, i) - 1))
                allocate(idum2(data_cols(2, i) - 1))

                do j = 1, n_sample
                    read(unit, *) idum1, variables(:, j), idum2, objectives(j)
                end do
            close(unit)

            call model_local(i)%initialize(variables, objectives)
            call model_local(i)%set_bounds(data_range(1, :), data_range(2, :))

            open(newunit=unit, file=theta_files(i), status="old")
                read(unit, *) scale
                read(unit, *) theta
            close(unit)

            call model_local(i)%set_scale(scale)
            call model_local(i)%calc_likelihood(theta)

            deallocate(variables, objectives, theta, scale, idum1, idum2)
        end do
    end subroutine load_krig


    real(8) function estimate(model_no, params) result(res)
        integer, intent(in) :: model_no
        real(8), intent(in) :: params(:)

        res = model(model_no)%estimate(params(:))
    end function estimate

end module mod_kriging
