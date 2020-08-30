!=======================================================================
!         6DoF EOM solver
!=======================================================================
module dof_kriging
  use util
  use individual
  use kriging
  use soga
  implicit none

  private
  public :: model, init_krig, load_krig, estimate

  type(TKriging) :: model(3)

  contains

  subroutine init_krig(num_var, pop_size, steps, scale)  !蛻昴ａ縺ｦ繧ｯ繝ｪ繧ｮ繝ｳ繧ｰ縺吶ｋ縺ｨ縺阪↓菴ｿ縺�繝ｫ繝ｼ繝√Φ
    implicit none
    integer, intent(in) :: num_var, pop_size, steps
    real(8), intent(in) :: scale(:)
    type(TSOGA) :: optimizer
    character(:), allocatable :: input_files(:), theta_files(:)
    real(8), allocatable :: variables(:, :), objectives(:), theta(:)
    integer :: num_sample
    integer :: unit, i, j, idum

    input_files = ["table_xl.txt", "table_ym.txt", "table_zn.txt"]
    theta_files = ["btheta_cx_2", "btheta_cm_2", "btheta_cz_2"]

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
      call model(i)%set_bounds([-5d0, -50d0, 0.1d0], [21d0, 10d0, 0.5d0])
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

  subroutine load_krig(num_var) !菴懊▲縺溘け繝ｪ繧ｮ繝ｳ繧ｰ繧偵Ο繝ｼ繝峨☆繧九Ν繝ｼ繝√Φ
    implicit none
    integer, intent(in) :: num_var

    character(:), allocatable :: input_files(:), theta_files(:)

    real(8), allocatable :: variables(:, :), objectives(:), theta(:), scale(:)
    integer :: num_sample
    integer :: unit, i, j, idum

    input_files = ["table_xl.txt", "table_ym.txt", "table_zn.txt"]
    theta_files = ["btheta_cx_2", "btheta_cm_2", "btheta_cz_2"]

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
      call model(i)%set_bounds([-5d0, -50d0, 0.1d0], [21d0, 10d0, 0.5d0])

      open(newunit=unit, file=theta_files(i), status="old")
        read(unit, *) scale
        read(unit, *) theta
      close(unit)

      call model(i)%set_scale(scale)
      call model(i)%calc_likelihood(theta)
      deallocate(variables, objectives, theta, scale)
    end do
  end subroutine load_krig

  real function estimate(model_no, aa, elv, ma) result(krig_out)
    implicit none
    integer, intent(in) :: model_no
    real, intent(in) :: aa, elv, ma
    krig_out = model(model_no)%estimate(dble([aa, elv, ma]))
  end function estimate
end module dof_kriging

module para
  implicit none

  real rad, deg, pi

  integer n, div, n_last
  real delt, tend

  real grav, onsoku, rho
  real m_body, s_ref, b_span, mac
  real roll, yaw, pitch

  real u_start, v_start, w_start
  real p_start, q_start, r_start
  real phi_start, tht_start, psi_start
  real x_start, y_start, z_start

  real dum

  real ix, iy, iz
  real ixy, ixz
  real iyz

  real ei, cx, cy, cz, cl, cm, cn, ini_elv, ini_thrust

  real ang_thrust
  !real gamt
  real kw, x1, xl !蠑ｷ縺輔�髟ｷ縺輔∽ｸｭ蠢�

  integer end_flag
  integer counter
  real elv_lower, elv_upper, elv_bound
  real thrust_lower, thrust_upper, thrust_bound
  character(64) :: file_out

  real, allocatable :: tt(:)
  real, allocatable :: uu(:), vv(:), ww(:), v_ini(:), ma(:), elv(:), thrust(:)
  real, allocatable :: ue(:), ve(:), we(:)
  real, allocatable :: jj(:), ll(:), li(:)

  real, allocatable :: pp(:), qq(:), rr(:), gamt(:)

  real, allocatable :: phip(:), thtt(:), psip(:), aa(:)

  real, allocatable :: xx(:), yy(:), zz(:)
  real, allocatable :: cxh(:), cmh(:), czh(:)
end module para

module dof_problem
  use para
  use individual
  use problem
  use dof_kriging
  implicit none

  type, extends(TProblem) :: TDoF
    integer :: start, end
    real :: lead_elm, lead_elm2

    contains

    generic :: initialize => initialize_dof
    procedure :: initialize_dof
    procedure :: call_indiv, calc_route4
  end type TDoF

  type(TDoF) :: problem

  contains

  subroutine initialize_dof(this, start, end, lead_elm, lead_elm2)
    implicit none
    class(TDoF), intent(inout) :: this
    integer, intent(in) :: start, end
    real, intent(in) :: lead_elm, lead_elm2

    this%start = start
    this%end = end
    this%lead_elm =  min(max((lead_elm + 20) / 40, 0.0), 1.0)
    this%lead_elm2 = min(max((lead_elm2),0.0), 1.0)

  end subroutine initialize_dof

  subroutine call_indiv(this, indiv)
    implicit none
    class(TDoF), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv

    if (indiv%evaluated) return

    call this%calc_route4(indiv%dvariables, indiv%objectives, indiv%constraints, indiv%feasible)  !calc_route 縺ｮ蜻ｼ縺ｳ蜃ｺ縺�

    indiv%evaluated = .true.
  end subroutine call_indiv

  subroutine calc_route4(this, variables, objectives, constraints, feasible)
    implicit none
    class(TDoF), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out), allocatable :: objectives(:), constraints(:)
    logical, intent(out) :: feasible
    real, allocatable :: variables_2(:,:)

    counter = counter + 1

    variables_2 = reshape(variables,[2,size(variables)/2])

    call set_elv(this%start, this%end, real(variables_2(1,:)))
    call set_thrust(this%start, this%end, real(variables_2(2,:)))
    call flight(this%start, this%end)

    allocate(objectives(2))
    !逶ｮ逧�髢｢謨ｰ
    allocate(constraints(2))
    !蛻ｶ邏�驕募渚

    objectives(1) = ga_obj_func(this%start, this%end)
    objectives(2) = sum(thrust)

    constraints(1) = max(aa(n_last + 1) - 20, 0.0, -(aa(n_last + 1) - (0)))
    constraints(2) = abs(xx(n_last + 1) - 4770)

    feasible = end_flag == 1 .and. xx(n_last) > 4769.5 .and. xx(n_last) < 4770.5
    !蛻ｶ邏�譚｡莉ｶ

  end subroutine calc_route4

  subroutine set_elv(start, end, elv_input)
    implicit none
    integer, intent(in) :: start, end
    real, intent(in) :: elv_input(:)
    real :: r, d_elv
    integer :: in_size, out_size, i

    in_size = size(elv_input)
    out_size = end - start + 1
    r = real(in_size) / out_size
    print *, start, end, in_size

    do i = 0, out_size - 1
      if (start + i == 0) cycle
       d_elv = (2 * elv_input(int(r * i) + 1) - 1) * elv_bound
       elv(start + i) = min(max(elv(start + i - 1) + (d_elv * delt), elv_lower), elv_upper)
      end do
  end subroutine set_elv

  subroutine set_thrust(start, end, thrust_input)
    implicit none
    integer, intent(in) :: start, end
    real, intent(in) :: thrust_input(:)
    real :: r, d_thrust
    integer :: in_size, out_size, i

    in_size = size(thrust_input)
    out_size = end - start + 1
    r = real(in_size) / out_size
    print *, start, end, in_size

    do i = 0, out_size - 1
      if (start + i == 0) cycle
       d_thrust = (2 * thrust_input(int(r * i) + 1) - 1) * thrust_bound
       thrust(start + i) = min(max(thrust(start + i - 1) + (d_thrust * delt), thrust_lower), thrust_upper)
      end do
  end subroutine set_thrust

  subroutine flight(start, end)
    implicit none
    integer, intent(in) :: start, end

    end_flag = 0

    do n = start, end
      if (mod(n, 5000) == 0 .and. n > start) print*, counter, n, "step calculating..."

      cx = estimate(1, aa(n), elv(n), ma(n))
      cm = estimate(2, aa(n), elv(n), ma(n))
      cz = estimate(3, aa(n), elv(n), ma(n))
      cxh(n) = cx
      cmh(n) = cm
      czh(n) = cz
      call rk4

      call rot_iner
      call check_conv(n, end_flag)
      if (end_flag > 0) exit
      n_last = n
    end do
  end subroutine flight

  real(8) function ga_obj_func(start, end) result(obj) !繧ｳ繧ｹ繝磯未謨ｰ縺ｮ螳夂ｾｩ
    implicit none
    integer, intent(in) :: start, end

    do n = start, n_last
      li(n) = 2 * ((((zz(n) - zz(0)) * cos(aa(0)) - ((xx(n) - xx(0)) * sin(aa(0)))))**2) &
            + 0.04 * ((elv(n) - elv(0)) **2)
      ll(n + 1) =  ll(n) + li(n) * (tt(n+1) - tt(n))
    end do

     obj = abs(gamt(n_last)) + ll(n_last)

  end function ga_obj_func
end module dof_problem

module dof_optimizer
  use util
  use soga
  use nsga2
  use nsga2c
  use moead
  use dof_problem
  implicit none

  ! type(TSOGA) :: optimizer
  ! type(TNSGA2) :: optimizer
  type(TNSGA2C) :: optimizer
  ! type(TMOEAD) :: optimizer
  ! integer :: num_var_, pop_size_
  integer :: num_step_

  contains

  subroutine init_optimizer(num_var, pop_size, num_step)
    implicit none
    integer, intent(in) :: num_var, pop_size, num_step


    ! num_var_ = num_var
    ! pop_size_ = pop_size
    num_step_ = num_step

    ! call optimizer%initialize(nx=num_var, N=pop_size, selection="T", crossover="BLX", mutation="PM") ! soga
    ! call optimizer%initialize(nx=num_var, m=2, N=pop_size, selection="T", crossover="BLX", mutation="PM") ! nsga2
    ! call optimizer%initialize(nx=num_var, m=2, N=pop_size, T=5, g="CH", crossover="BLX", mutation="PM") ! moead
     call optimizer%initialize(nx=num_var, m=2,N=pop_size, selection="T", crossover="BLX", mutation="PM") !nsga2c
    ! call optimizer%set_fitness_type("VALUE")
    optimizer%elite_preservation = .false.
    optimizer%sharing = .true.
    optimizer%history_preservation = .true.

  end subroutine init_optimizer

  subroutine optimize(problem, result)
    implicit none
    type(TDoF), intent(in) :: problem
    real, intent(out), allocatable :: result(:)
    real(8), allocatable :: best(:)

    call optimizer%set_problem(problem)
    call optimizer%prepare_calculation
    call optimizer%run(num_step_)
    call optimizer%best(best)

    result = real(best, kind=4)
  end subroutine optimize

  subroutine save_result
    implicit none

    call optimizer%save_result("out_so_X_1.csv", elite="only", feasible="all")
    call optimizer%save_result("out_so_X_2.csv", elite="all", feasible="all")
    call optimizer%save_history("out_so_X_3.csv", elite="only", feasible="all")
    call optimizer%save_history("out_so_X_4.csv", elite="all", feasible="all")
    call optimizer%save_history("out_so_X_5.csv", elite="all", feasible="only")
  end subroutine save_result
end module dof_optimizer

!=======================================================================
program main     ! 6DoF main !
!=======================================================================
  use para
  use dof_kriging
  ! use dof_problem
  ! use dof_optimizer
  use util
  implicit none

  integer :: start_time

  ! ============================================================================
  ! initialize
  ! ============================================================================

  if (.false.) then
    ! ﾎｸ蛻晄悄蛹也畑
    call init_krig(num_var=3, pop_size=100, steps=100, scale=[100d0, 150d0, 10d0])
    stop
  else
    call load_krig(num_var=3)
    !call chech_kriging
  end if

  pi=4*atan(1.0) !definition of pi
  rad = pi/180        !convert [deg]ﾂ�ﾂｨ[rad]
  deg = 180/pi        !convert [rad]ﾂ�ﾂｨ[deg]

  call ini_cond

  div = tend / delt
  n_last = 0

  allocate (tt(0:div))
  allocate (uu(0:div), vv(0:div), ww(0:div))
  allocate (v_ini(0:div), ma(0:div))
  allocate (ue(0:div), ve(0:div), we(0:div))
  allocate (pp(0:div), qq(0:div), rr(0:div))
  allocate (phip(0:div), thtt(0:div), psip(0:div), aa(0:div))
  allocate (xx(0:div), yy(0:div), zz(0:div))
  allocate (gamt(0:div))
  allocate (cxh(0:div), cmh(0:div), czh(0:div))
  allocate (elv(0:div), thrust(0:div))
  allocate (jj(0:div), li(0:div), ll(0:div))


  tt(0) = 0

  uu(0) = u_start
  vv(0) = v_start
  ww(0) = w_start
  v_ini(0) = sqrt(uu(0)**2+vv(0)**2+ww(0)**2)

  pp(0) = p_start
  qq(0) = q_start
  rr(0) = r_start

  phip(0) = phi_start
  thtt(0) = tht_start
  psip(0) = psi_start

  xx(0) = x_start
  yy(0) = y_start
  zz(0) = z_start

  call cal_init

  ! ============================================================================
  ! main loop
  ! ============================================================================

  call system_clock(start_time)

  call main_2(input="input_elv.csv")

  call elapsed_time(start_time)

  contains

  subroutine main_2(input)
    implicit none
    character(*), intent(in) :: input
    !character(*) :: input_thrust.csv
    integer :: start, end

    call ini_cond
    tt(0) = 0
    xx(0) = x_start
    zz(0) = z_start

    call load_elv(input)
    call load_thrust(input)
    start = 0
    end = 50
    call flight(start, end)
    call out_total(start, end - 1)
  end subroutine main_2

  ! subroutine load_ga_result(file, lines, num_offset, z_bound)
  !   implicit none
  !   character(*), intent(in) :: file
  !   integer, intent(in) :: lines, num_offset
  !   real, intent(in) :: z_bound
  !   real, allocatable :: offset(:), var(:), var_2(:,:)
  !   integer :: unit, start, end, num_var, i

  !   start = 0
  !   end = div
  !   num_var = 800

  !   call ini_cond
  !   tt(0)=0

  !   allocate(offset(num_offset))
  !   allocate(var(num_var))

  !   open(newunit=unit, file=file)
  !     do i = 1, lines - 1
  !       read(unit, *)
  !     end do
  !     read(unit, *) offset, var
  !   close(unit)

  !   var_2 = reshape(var,[2, size(var)/2])
  !   call set_elv(start, end, var_2(1,:))
  !   call set_thrust(start, end, var_2(2,:))
  !   call flight(start, end)
  !   call out_total(start, end - 1)
  ! end subroutine load_ga_result

  ! subroutine ga_landing(start, end)
  !   implicit none
  !   integer, intent(in) :: start, end
  !   real, allocatable :: best(:), best_2(:,:)

  !   if (start >= div) return
  !   counter = 0 ! 陦ｨ遉ｺ逕ｨ

  !   call problem%initialize(start, end, elv(start), thrust(start))
  !   call optimize(problem, best)
  !   print *,best
  !   best_2 = reshape(best,[2, size(best)/2])
  !   call set_elv(start, end, real(best_2(1,:)))
  !   call set_thrust(start, end, real(best_2(2,:)))
  !   call flight(start, end)
  !   call elapsed_time(start_time)
  !   call out_total(start, end - 1)
  !   call save_result
  ! end subroutine ga_landing

  subroutine flight(start, end)
    implicit none
    integer, intent(in) :: start, end

    end_flag = 0

    do n = start, end
      if (mod(n, 5000) == 0 .and. n > start) print*, counter, n, "step calculating..."

      cx = estimate(1, aa(n), elv(n), ma(n))
      cm = estimate(2, aa(n), elv(n), ma(n))
      cz = estimate(3, aa(n), elv(n), ma(n))
      cxh(n) = cx
      cmh(n) = cm
      czh(n) = cz
      call rk4

      call rot_iner
      call check_conv(n, end_flag)
      if (end_flag > 0) exit
      n_last = n
    end do
  end subroutine flight
end program main
!=======================================================================

! 蛻晄悄蛟､險ｭ螳�
subroutine ini_cond
  use para
  implicit none

  open(55, file='initial_cond.inp', status="old")
    read(55, *) dum
    read(55, *) rho, grav
    read(55, *) onsoku
    read(55, *) m_body
    read(55, *) s_ref
    read(55, *) mac, b_span
    read(55, *) ix, iy, iz
    read(55, *) ixy, ixz
    read(55, *) iyz
    read(55, *) tend, delt
    read(55, *) u_start, v_start, w_start
    read(55, *) p_start, q_start, r_start
    read(55, *) phi_start, tht_start, psi_start
    read(55, *) x_start, y_start, z_start
    read(55, *) file_out
    read(55, *) ini_elv
    read(55, *) ini_thrust
  close(55)

  elv_lower = -50.0
  elv_upper = 30.0
  elv_bound = 2.0

  thrust_lower = -140000
  thrust_upper =  140000
  thrust_bound =  1000

end subroutine ini_cond

!=======================================================================
subroutine cal_init
  use para
  implicit none

  aa(0)=deg*atan(ww(0)/uu(0))
  gamt(0)=deg*atan(we(0)/ue(0))
  ma(0)=v_ini(0)/onsoku
  elv = 0
  elv(0) = ini_elv
  thrust = 0
  thrust(0) = ini_thrust

  cx = 0.0
  cy = 0.0
  cz = 0.0
  cl = 0.0
  cm = 0.0
  cn = 0.0
  cxh = 0.0
  cmh = 0.0
  czh = 0.0

  counter = 0 ! 陦ｨ遉ｺ逕ｨ

end subroutine cal_init

subroutine load_elv(file)
  use para
  implicit none
  character(*), intent(in) :: file
  real, allocatable :: elv_input(:)
  integer, allocatable :: steps(:)
  integer :: lines, unit, i, j, k

  open(newunit=unit, file=file, status="old")
    read(unit, *) lines
    allocate(elv_input(lines), steps(lines))
    do i = 1, lines
      read(unit, *) steps(i), elv_input(i)
    end do
  close(unit)

  do i = 1, lines
    if (i < lines) then
      k = min(steps(i+1)-1, div)
    else
      k = div
    end if
    do j = steps(i), k
      if (j == 0) cycle
      elv(j) = min(max(elv(j - 1) + elv_input(i) * delt, elv_lower), elv_upper)
    end do
  end do
end subroutine load_elv

subroutine load_thrust(file)
  use para
  implicit none
  character(*), intent(in) :: file
  real, allocatable :: thrust_input(:)
  integer, allocatable :: steps(:)
  integer :: lines, unit, i, j, k

  open(newunit=unit, file=file, status="old")
    read(unit, *) lines
    allocate(thrust_input(lines), steps(lines))
    do i = 1, lines
      read(unit, *) steps(i), thrust_input(i)
    end do
  close(unit)

  do i = 1, lines
    if (i < lines) then
      k = min(steps(i+1)-1, div)
    else
      k = div
    end if
    do j = steps(i), k
      if (j == 0) cycle
      thrust(j) = min(max(thrust(j - 1) + thrust_input(i) * delt, thrust_lower), thrust_upper)
    end do
  end do
end subroutine load_thrust

!=======================================================================
!     function of EOM ( 6 degrees of freedom )
!=======================================================================

!-------translational motion
real function fu(t, q, w, r, v, tht)
  use para
  implicit none
  real t, q, w, r, v, tht
  fu = (((-q*w/deg)+(r*v/deg))-(grav*sin(tht*rad)) &
  & +(((thrust(n))/m_body)*cos(ang_thrust*rad))&
  & +((rho*(v_ini(n)**2)*s_ref*cx)/(2*m_body)))
end function fu

real function fv(t, r, u, p, w, tht, phi)
  use para
  implicit none
  real t, r, u, p, w, tht, phi
  fv = ((-(r*u/deg)+(p*w/deg))+(grav*cos(tht*rad)*sin(phi*rad)) &
  & +((rho*(v_ini(n)**2)*s_ref*cy)/(2*m_body)))
end function fv

real function fw(t, p, v, q, u, tht, phi)
  use para
  implicit none
  real t, p, v, q, u, tht, phi
  fw = ((-(p*v/deg)+(q*u/deg))+(grav*cos(tht*rad)*cos(phi*rad)) &
  & +(((thrust(n))/m_body)*sin(ang_thrust*rad))&
  & +((rho*(v_ini(n)**2)*s_ref*cz)/(2*m_body)))
end function fw
!-------angular acceleration
real function fp(t, p, q, r)
  use para
  implicit none
  real t, p, q, r
  roll = ((iy-iz)*q*r/deg)+(ixz*p*q/deg) &
  & +((rho*(v_ini(n)**2)*s_ref*b_span*deg*cl)/2)
  yaw = ((ix-iy)*p*q/deg)-(ixz*q*r/deg) &
  & +((rho*(v_ini(n)**2)*s_ref*b_span*deg*cn)/2)
  fp = (((roll/ix)+(ixz/ix)*(yaw/iz))/(1.d0-((ixz**2)/(ix*iz))))
end function fp

real function fq(t, p, r)
  use para
  implicit none
  real t, r, p
  pitch = (((iz-ix)*r*p/deg)+((ixz*(r**2-p**2))/deg)) &
  & +((rho*(v_ini(n)**2)*s_ref*mac*deg*cm)/2)
  fq = (pitch/iy)
end function fq

real function fr(t, p, q, r)
  use para
  implicit none
  real t, p, q, r
  roll = ((iy-iz)*q*r/deg)+(ixz*p*q/deg) &
  & +((rho*(v_ini(n)**2)*s_ref*b_span*deg*cl)/2)
  yaw = ((ix-iy)*p*q/deg)-(ixz*q*r/deg) &
  & +((rho*(v_ini(n)**2)*s_ref*b_span*deg*cn)/2)
  fr = (((yaw/iz)+(ixz/iz)*(roll/ix))/(1.d0-((ixz**2)/(ix*iz))))
end function fr
!-------attitude angle
real function fpsi(t, q, r, phi, tht)
  use para
  implicit none
  real t, r, q, phi, tht
  fpsi = ((r*cos(phi*rad)+q*sin(phi*rad))/(cos(tht*rad)))
end function fpsi

real function ftht(t, q, r, phi, tht)
  use para
  implicit none
  real t, r, q, phi, tht
  ftht = (q*cos(phi*rad)-r*sin(phi*rad))
end function ftht

real function fphi(t, p, q, r, phi, tht)
  use para
  implicit none
  real t, p, q, r, phi, tht
  fphi =(p+(r*cos(phi*rad)+q*sin(phi*rad))*tan(tht*rad))
end function fphi

!=======================================================================

subroutine rk4
  use para
  implicit none

  real t, fu, fv, fw, fp, fq, fr, fphi, ftht, fpsi
  real u, v, w
  real p, q, r
  real phi, tht, psi
  real xxx, zzz

  real k1_u, k1_v, k1_w
  real k2_u, k2_v, k2_w
  real k3_u, k3_v, k3_w
  real k4_u, k4_v, k4_w

  real k1_p, k1_q, k1_r
  real k2_p, k2_q, k2_r
  real k3_p, k3_q, k3_r
  real k4_p, k4_q, k4_r

  real k1_phi, k1_tht, k1_psi
  real k2_phi, k2_tht, k2_psi
  real k3_phi, k3_tht, k3_psi
  real k4_phi, k4_tht, k4_psi

  t = tt(n)

  xxx = xx(n)
  zzz = zz(n)

  p = pp(n)
  q = qq(n)
  r = rr(n)

  phi = phip(n)
  tht = thtt(n)
  psi = psip(n)

  u = uu(n)
  v = vv(n)
  w = ww(n)

  !-------------------k1
  k1_u = delt * fu(t, q, w, r, v, tht)
  k1_v = delt * fv(t, r, u, p, w, tht, phi)
  k1_w = delt * fw(t, p, v, q, u, tht, phi)

  k1_p = delt * fp(t, p, q, r)
  k1_q = delt * fq(t, p, r)
  k1_r = delt * fr(t, p, q, r)

  k1_psi = delt * fpsi(t, q, r, phi, tht)
  k1_tht = delt * ftht(t, q, r, phi, tht)
  k1_phi = delt * fphi(t, p, q, r, phi, tht)
  !-------------------k2
  k2_u = delt * fu(t+delt/2, q+k1_q/2, w+k1_w/2, r+k1_r/2, v+k1_v/2, tht+k1_tht/2)
  k2_v = delt * fv(t+delt/2, r+k1_r/2, u+k1_u/2, p+k1_p/2, w+k1_w/2, tht+k1_tht/2, phi+k1_phi/2)
  k2_w = delt * fw(t+delt/2, p+k1_p/2, v+k1_v/2, q+k1_q/2, u+k1_u/2, tht+k1_tht/2, phi+k1_phi/2)

  k2_p = delt * fp(t+delt/2, p+k1_p/2, q+k1_q/2, r+k1_r/2)
  k2_q = delt * fq(t+delt/2, p+k1_p/2, r+k1_r/2)
  k2_r = delt * fr(t+delt/2, p+k1_p/2, q+k1_q/2, r+k1_r/2)

  k2_psi = delt * fpsi(t+delt/2, q+k1_q/2, r+k1_r/2, phi+k1_phi/2, tht+k1_tht/2)
  k2_tht = delt * ftht(t+delt/2, q+k1_q/2, r+k1_r/2, phi+k1_phi/2, tht+k1_tht/2)
  k2_phi = delt * fphi(t+delt/2, p+k1_p/2, q+k1_q/2, r+k1_r/2, phi+k1_phi/2, tht+k1_tht/2)
  !-------------------k3
  k3_u = delt * fu(t+delt/2, q+k2_q/2, w+k2_w/2, r+k2_r/2, v+k2_v/2, tht+k2_tht/2)
  k3_v = delt * fv(t+delt/2, r+k2_r/2, u+k2_u/2, p+k2_p/2, w+k2_w/2, tht+k2_tht/2, phi+k2_phi/2)
  k3_w = delt * fw(t+delt/2, p+k2_p/2, v+k2_v/2, q+k2_q/2, u+k2_u/2, tht+k2_tht/2, phi+k2_phi/2)

  k3_p = delt * fp(t+delt/2, p+k2_p/2, q+k2_q/2, r+k2_r/2)
  k3_q = delt * fq(t+delt/2, p+k2_p/2, r+k2_r/2)
  k3_r = delt * fr(t+delt/2, p+k2_p/2, q+k2_q/2, r+k2_r/2)

  k3_psi = delt * fpsi(t+delt/2, q+k2_q/2, r+k2_r/2, phi+k2_phi/2, tht+k2_tht/2)
  k3_tht = delt * ftht(t+delt/2, q+k2_q/2, r+k2_r/2, phi+k2_phi/2, tht+k2_tht/2)
  k3_phi = delt * fphi(t+delt/2, p+k2_p/2, q+k2_q/2, r+k2_r/2, phi+k2_phi/2, tht+k2_tht/2)
  !-------------------k4
  k4_u = delt * fu(t+delt, q+k3_q, w+k3_w, r+k3_r, v+k3_v, tht+k3_tht)
  k4_v = delt * fv(t+delt, r+k3_r, u+k3_u, p+k3_p, w+k3_w, tht+k3_tht, phi+k3_phi)
  k4_w = delt * fw(t+delt, p+k3_p, v+k3_v, q+k3_q, u+k3_u, tht+k3_tht, phi+k3_phi)

  k4_p = delt * fp(t+delt, p+k3_p, q+k3_q, r+k3_r)
  k4_q = delt * fq(t+delt, p+k3_p, r+k3_r)
  k4_r = delt * fr(t+delt, p+k3_p, q+k3_q, r+k3_r)

  k4_psi = delt * fpsi(t+delt, q+k3_q, r+k3_r, phi+k3_phi, tht+k3_tht)
  k4_tht = delt * ftht(t+delt, q+k3_q, r+k3_r, phi+k3_phi, tht+k3_tht)
  k4_phi = delt * fphi(t+delt, p+k3_p, q+k3_q, r+k3_r, phi+k3_phi, tht+k3_tht)
  !-------------------output for eis-------------------------------------

  u = u + (k1_u+(2*k2_u)+(2*k3_u)+k4_u)/6.d0
  v = v + (k1_v+(2*k2_v)+(2*k3_v)+k4_v)/6.d0
  w = w + (k1_w+(2*k2_w)+(2*k3_w)+k4_w)/6.d0

  p = p + (k1_p+(2*k2_p)+(2*k3_p)+k4_p)/6.d0
  q = q + (k1_q+(2*k2_q)+(2*k3_q)+k4_q)/6.d0
  r = r + (k1_r+(2*k2_r)+(2*k3_r)+k4_r)/6.d0

  phi = phi + (k1_phi+(2*k2_phi)+(2*k3_phi)+k4_phi)/6.d0
  tht = tht + (k1_tht+(2*k2_tht)+(2*k3_tht)+k4_tht)/6.d0
  psi = psi + (k1_psi+(2*k2_psi)+(2*k3_psi)+k4_psi)/6.d0

  tt(n+1) = t+delt

  uu(n+1) = u
  vv(n+1) = v
  ww(n+1) = w
  v_ini(n+1) = sqrt((u**2)+(v**2)+(w**2))

  pp(n+1) = p
  qq(n+1) = q
  rr(n+1) = r

  phip(n+1) = phi
  thtt(n+1) = tht
  psip(n+1) = psi

  aa(n+1) = deg*atan(ww(n+1)/uu(n+1))    !for next inp(eis)
  ma(n+1) = v_ini(n+1)/onsoku !mach num for next step

end subroutine rk4

subroutine check_conv(i, flag)
  use para
  implicit none
  integer, intent(in) :: i       ! => the number of calculated steps
  integer, intent(out) :: flag   ! => TRUE then calculation finish

  flag = 0
  if (zz(i) <= 0) then
    flag = 1
    print *, "zz= ", zz(i)
    print *, "we= ", we(i)
    print *, "ue= ", ue(i)
    print *, "xx= ", xx(i)
    print *, "ma= ", ma(i)
  else if (aa(i) < 0 .or. aa(i) > 20) then
    flag = 2
    print *, "aa= ", aa(i)
  else if (ma(i) < 0.1 .or. ma(i) > 0.5) then
    flag = 3
    print *, "ma= ", ma(i)
  else if (elv(i) < -50 .or. elv(i) > 30) then
    flag = 4
    print *, "elv= ", elv(i)
  else
    flag = 0
  end if
end subroutine check_conv

!======================================================================

subroutine rot_iner
  use para
  implicit none

    ue(n) = (cos(thtt(n)*rad)*cos(psip(n)*rad))*uu(n) &
    & +(sin(phip(n)*rad)*sin(thtt(n)*rad)*cos(psip(n)*rad) &
    & -cos(phip(n)*rad)*sin(psip(n)*rad))*vv(n) &
    & +(cos(phip(n)*rad)*sin(thtt(n)*rad)*cos(psip(n)*rad) &
    & +sin(phip(n)*rad)*sin(psip(n)*rad))*ww(n)

    ve(n) = (cos(thtt(n)*rad)*sin(psip(n)*rad))*uu(n) &
    & +(sin(phip(n)*rad)*sin(thtt(n)*rad)*sin(psip(n)*rad) &
    & +cos(phip(n)*rad)*cos(psip(n)*rad))*vv(n) &
    & +(cos(phip(n)*rad)*sin(thtt(n)*rad)*sin(psip(n)*rad) &
    & -sin(phip(n)*rad)*cos(psip(n))*rad)*ww(n)

    we(n) = (sin(thtt(n)*rad))*uu(n) &
    & +(-sin(phip(n)*rad)*cos(thtt(n)*rad))*vv(n) &
    & +(-cos(phip(n)*rad)*cos(thtt(n)*rad))*ww(n)

    gamt(n) = deg*atan(we(n) / ue(n))
    !======================================================================
    !for trajectory

    xx(n+1) = xx(n) + ue(n)*delt
    yy(n+1) = yy(n) + ve(n)*delt
    zz(n+1) = zz(n) + we(n)*delt

end subroutine rot_iner

!======================================================================

subroutine out_total(start, end)
  use para
  implicit none
  integer, intent(in) :: start, end
  character(:), allocatable :: format1,format2
  integer :: unit1,unit2
  integer :: las

  format1 = "(17(es15.8',')es15.8)"

  if (start == 0) then
    open(newunit=unit1, file=file_out, status='replace')
    write(unit1, "(a)") "time,theta,alpha,gamt,elv,ma,x,z,ue,we,uu,ww,cx,cm,cz,q,thrust"
  else
    open(newunit=unit1, file=file_out, position='append')
  end if

  do n = start, n_last
    li(n) = 2 * ((((zz(n) - zz(0)) * cos(aa(0)) - ((xx(n) - xx(0)) * sin(aa(0)))))**2) &
          + 0.04 * ((elv(n) - elv(0)) **2)
    ll(n+1) =  ll(n) + li(n) * (tt(n+1) - tt(n))
    jj(n) = abs(gamt(n) *rad) + ll(n)
  end do

  las = min(end, n_last)
  do n = start, n_last
    if (mod(n, 1) == 0) then
    write(unit1, format1) tt(n), thtt(n), aa(n), gamt(n), elv(n), ma(n), xx(n), zz(n), &
                        ue(n), we(n), uu(n), ww(n), cxh(n), cmh(n), czh(n), qq(n), thrust(n)
  end if
  end do
  close(unit1)

  format2 = "(3(es15.8','))"
  open(newunit=unit2, file='J.csv')
  write(unit2,  '(a)') "time,J,gamt"
  do n = start, las
    if (mod(n, 50) == 0 .or. n >= int(las / 50) * 50) then
    write(unit2, format2) tt(n), ll(n), gamt(n)
  end if
  end do

  close(unit2)
end subroutine out_total

subroutine chech_kriging
  use dof_kriging
  implicit none

  real :: a, e, m, y1, y2, y3
  integer :: unit1, unit2, unit3, i, j

  open(newunit=unit1, file="test1.csv", status="replace")
  open(newunit=unit2, file="test2.csv", status="replace")
  open(newunit=unit3, file="test3.csv", status="replace")
    do i = 0, 100
      if (mod(i, 10) == 0) print *, i
      do j = 0, 10
        a = -5 + i * 0.26
        e = 5
        m = j * 0.04 + 0.1
        y1 = estimate(1, a, e, m)
        y2 = estimate(2, a, e, m)
        y3 = estimate(3, a, e, m)
        write(unit1, "(3(f0.10','))") m, a, y1
        write(unit2, "(3(f0.10','))") m, a, y2
        write(unit3, "(3(f0.10','))") m, a, y3
      end do
    end do
  close(unit3)
  close(unit2)
  close(unit1)

  stop
end subroutine chech_kriging

