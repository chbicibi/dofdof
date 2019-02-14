module dof_kriging
  use util
  use individual
  use kriging
  use soga
  implicit none

  private
  public :: model, init_krig, load_krig, estimate, estimate2

  type(TKriging) :: model(3)

  contains

  subroutine init_krig(num_var, pop_size, steps, scale)
    implicit none
    integer, intent(in) :: num_var, pop_size, steps
    real(8), intent(in) :: scale(:)
    type(TSOGA) :: optimizer
    character(:), allocatable :: input_files(:), theta_files(:)
    real(8), allocatable :: variables(:, :), objectives(:), theta(:)
    integer :: num_sample
    integer :: unit, i, j, idum

    input_files = ["table_xl.txt", "table_ym.txt", "table_zn.txt"]
    theta_files = ["btheta_cx", "btheta_cm", "btheta_cz"]

    call optimizer%initialize(nx=num_var, N=pop_size, selection="R", crossover="BLX", mutation="PM")
    optimizer%elite_preservation = .true.
    optimizer%sharing = .true.

    do i = 1, 3
      open(newunit=unit, file=input_files(i), status="old")
        read(unit, *) num_sample
        allocate(variables(num_var, num_sample))
        allocate(objectives(num_sample))
        allocate(theta(num_var))
        do j = 1, num_sample
          if (i == 2) then
            read(unit, *) idum, variables(:, j), idum, objectives(j)
          else
            read(unit, *) idum, variables(:, j), objectives(j)
          end if
        end do
      close(unit)

      call model(i)%initialize(variables, objectives)
      call model(i)%set_bounds([-5d0, -40d0, 0.1d0], [10d0, 40d0, 0.5d0])
      call model(i)%set_scale(scale)

      call optimizer%set_problem(model(i))
      call optimizer%prepare_calculation
      call optimizer%run(steps)

      ! index = minloc([(optimizer%population(i)%indiv%o1(), i = 1, pop_size)], dim=1)
      ! theta = optimizer%population(index)%indiv%variables
      call optimizer%best(theta)

      open(newunit=unit, file=theta_files(i), status="replace")
        write(unit, *) model(i)%scale
        write(unit, *) theta
      close(unit)

      call model(i)%calc_likelihood(theta)
      deallocate(variables, objectives, theta)
    end do
  end subroutine init_krig

  subroutine load_krig(num_var)
    implicit none
    integer, intent(in) :: num_var

    character(:), allocatable :: input_files(:), theta_files(:)

    real(8), allocatable :: variables(:, :), objectives(:), theta(:), scale(:)
    integer :: num_sample
    integer :: unit, i, j, idum

    input_files = ["table_xl.txt", "table_ym.txt", "table_zn.txt"]
    theta_files = ["btheta_cx", "btheta_cm", "btheta_cz"]

    do i = 1, 3
      open(newunit=unit, file=input_files(i), status="old")
        read(unit, *) num_sample
        allocate(variables(num_var, num_sample))
        allocate(objectives(num_sample))
        allocate(theta(num_var))
        allocate(scale(num_var))
        do j = 1, num_sample
          if (i == 2) then
            read(unit, *) idum, variables(:, j), idum, objectives(j)
          else
            read(unit, *) idum, variables(:, j), objectives(j)
          end if
        end do
      close(unit)

      call model(i)%initialize(variables, objectives)
      call model(i)%set_bounds([-5d0, -40d0, 0.1d0], [10d0, 40d0, 0.5d0])

      open(newunit=unit, file=theta_files(i), status="old")
        read(unit, *) scale
        read(unit, *) theta
      close(unit)

      call model(i)%set_scale(scale)
      call model(i)%calc_likelihood(theta)
      deallocate(variables, objectives, theta, scale)
    end do
  end subroutine load_krig

  real(8) function estimate(model_no, aa, elv, ma) result(krig_out)
    implicit none
    integer, intent(in) :: model_no
    real(8), intent(in) :: aa, elv, ma

    krig_out = model(model_no)%estimate(dble([aa, elv, ma]))
  end function estimate

  real(8) function estimate2(model_no, aa, elv, ma) result(krig_out)
    implicit none
    integer, intent(in) :: model_no
    real(8), intent(in) :: aa, elv, ma
    character(:), allocatable :: input_files(:), theta_files(:)
    real(8) :: dum
    integer :: unit

    input_files = ["table_xl.txt", "table_ym.txt", "table_zn.txt"]
    theta_files = ["btheta_cx", "btheta_cm", "btheta_cz"]

    open(newunit=unit, file='inp', status='replace')
    write(unit, *) aa
    write(unit, *) elv
    write(unit, *) ma
    close(unit)

    call system('cp -r ' // input_files(model_no) // ' table.txt')

    open(newunit=unit, file='ei_inp', status='replace')
    write(unit, *) 1
    close(unit)

    call system('cp -r ' // theta_files(model_no) // ' ./btheta')
    call system('eis.exe<ei_inp')
    open(newunit=unit, file='out')
    read(unit, *) dum, krig_out, dum
    close(unit)
    call system('rm -rf btheta')
  end function estimate2
end module dof_kriging
