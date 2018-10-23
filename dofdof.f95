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
  public :: model, init_krig, load_krig, estimate, estimate2

  type(TKriging), allocatable :: model(:)
  character(*), parameter :: krig_input_file = "kriging.inp"
  integer, parameter :: ga_pop_size = 100
  integer, parameter :: ga_steps = 10

  contains

  subroutine load_inputfile(num_model, num_var, scales, lbounds, ubounds, table_files, theta_files)
    implicit none
    integer, intent(out) :: num_model, num_var
    real(8), intent(out), allocatable :: scales(:), lbounds(:), ubounds(:)
    character(:), intent(out), allocatable :: table_files(:), theta_files(:)
    integer :: unit, i

    open(newunit=unit, file=krig_input_file, status="old")
    read(unit, *) num_model
    allocate(model(num_model))
    allocate(character(60) :: table_files(num_model))
    allocate(character(60) :: theta_files(num_model))
    do i = 1, num_model
      read(unit, *) table_files(i)
    end do
    do i = 1, num_model
      read(unit, *) theta_files(i)
    end do
    read(unit, *) num_var
    allocate(scales(num_var))
    allocate(lbounds(num_var))
    allocate(ubounds(num_var))
    do i = 1, num_var
      read(unit, *) scales(i)
    end do
    do i = 1, num_var
      read(unit, *) lbounds(i), ubounds(i)
    end do
  end subroutine load_inputfile

  subroutine init_krig
    implicit none
    type(TSOGA) :: optimizer
    character(:), allocatable :: table_files(:), theta_files(:)
    real(8), allocatable :: variables(:, :), objectives(:), theta(:)
    real(8), allocatable :: scales(:), lbounds(:), ubounds(:)
    integer :: num_model, num_var, num_sample
    integer :: unit, i, j, k, idummy, nskip

    call load_inputfile(num_model, num_var, scales, lbounds, ubounds, table_files, theta_files)
    call optimizer%initialize(nx=num_var, N=ga_pop_size, selection="R", crossover="BLX", mutation="PM")
    optimizer%elite_preservation = .true.
    optimizer%sharing = .true.
    allocate(theta(num_var))

    do i = 1, num_model
      open(newunit=unit, file=table_files(i), status="old")
        read(unit, *) num_sample
        allocate(variables(num_var, num_sample))
        allocate(objectives(num_sample))
        do j = 1, num_sample
          if (i == 2) then
            nskip = 1
          else
            nskip = 0
          end if
          read(unit, *) idummy, variables(:, j), (idummy, k=1, nskip), objectives(j)
        end do
      close(unit)

      call model(i)%initialize(variables, objectives)
      call model(i)%set_bounds(lbounds, ubounds)
      call model(i)%set_scale(scales)

      call optimizer%set_problem(model(i))
      call optimizer%prepare_calculation
      call optimizer%run(ga_steps)
      call optimizer%best(theta)

      open(newunit=unit, file=theta_files(i), status="replace")
        write(unit, *) model(i)%scale
        write(unit, *) theta
      close(unit)

      call model(i)%calc_likelihood(theta)
      deallocate(variables, objectives)
    end do
  end subroutine init_krig

  subroutine load_krig
    implicit none
    character(:), allocatable :: table_files(:), theta_files(:)
    real(8), allocatable :: variables(:, :), objectives(:), theta(:)
    real(8), allocatable :: scales(:), lbounds(:), ubounds(:), dummy(:)
    integer :: num_model, num_var, num_sample
    integer :: unit, i, j, k, idummy, nskip

    call load_inputfile(num_model, num_var, dummy, lbounds, ubounds, table_files, theta_files)
    allocate(theta(num_var))
    allocate(scales(num_var))

    do i = 1, num_model
      open(newunit=unit, file=table_files(i), status="old")
        read(unit, *) num_sample
        print *, num_sample
        allocate(variables(num_var, num_sample))
        allocate(objectives(num_sample))
        do j = 1, num_sample
          if (i == 2) then
            nskip = 1
          else
            nskip = 0
          end if
          read(unit, *) idummy, variables(:, j), (idummy, k=1, nskip), objectives(j)
        end do
      close(unit)

      call model(i)%initialize(variables, objectives)
      call model(i)%set_bounds(lbounds, ubounds)

      open(newunit=unit, file=theta_files(i), status="old")
        read(unit, *) scales
        read(unit, *) theta
      close(unit)

      call model(i)%set_scale(scales)
      call model(i)%calc_likelihood(theta)
      deallocate(variables, objectives)
    end do
  end subroutine load_krig

  real function estimate(model_no, variables) result(krig_out)
    implicit none
    integer, intent(in) :: model_no
    real, intent(in) :: variables(:)

    krig_out = model(model_no)%estimate(dble(variables))
  end function estimate

  real function estimate2(model_no, aa, elv, ma) result(krig_out)
    implicit none
    integer, intent(in) :: model_no
    real, intent(in) :: aa, elv, ma
    character(:), allocatable :: table_files(:), theta_files(:)
    real :: dum
    integer :: unit

    table_files = ["table_xl.txt", "table_ym.txt", "table_zn.txt"]
    theta_files = ["btheta_cx", "btheta_cm", "btheta_cz"]

    open(newunit=unit, file='inp', status='replace')
    write(unit, *) aa
    write(unit, *) elv
    write(unit, *) ma
    close(unit)

    call system('cp -r ' // table_files(model_no) // ' table.txt')

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

  real ei, cx, cy, cz, cl, cm, cn, ini_elv
  real cx_pred, cy_pred, cz_pred
  real cl_pred, cm_pred, cn_pred
  logical end_flag
  integer counter
  real tar_thtt, elv_lower, elv_upper, elv_bound
  character(64) :: file_out

  real, allocatable :: tt(:)
  real, allocatable :: uu(:), vv(:), ww(:), v_ini(:), ma(:), elv(:)
  real, allocatable :: ue(:), ve(:), we(:)

  real, allocatable :: pp(:), qq(:), rr(:)

  real, allocatable :: phip(:), thtt(:), psip(:), aa(:)

  real, allocatable :: xx(:), yy(:), zz(:)
  real, allocatable :: cxh(:), cmh(:), czh(:)
end module para


!=======================================================================
program main     ! 6DoF main !
!=======================================================================
  use para
  use dof_kriging
  use util
  implicit none

  integer :: start_time
  integer :: is_init_krig

  ! ============================================================================
  ! initialize
  ! ============================================================================

  print *, 'init? 0=read btheta 1=init btheta 2=check kriging'
  read *, is_init_krig

  if (is_init_krig == 0) then
    call load_krig
  else if (is_init_krig == 1) then
    call init_krig
  else if (is_init_krig == 2) then
    call load_krig
    call chech_kriging
    stop
  else
    print *, "Error: invalid number"
    call exit(1)
  end if

  pi=4*atan(1.0) !definition of pi
  rad = pi/180        !convert [deg]¨[rad]
  deg = 180/pi        !convert [rad]¨[deg]

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
  allocate (cxh(0:div), cmh(0:div), czh(0:div))
  allocate (elv(0:div))


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
  call system_clock(start_time)
  call main_1
  call elapsed_time(start_time)

  contains

  subroutine main_1
    implicit none

    call flight(0, div)
    call out_total(0, div - 1)
  end subroutine main_1

  subroutine flight(start, end)
    implicit none
    integer, intent(in) :: start, end

    end_flag = .false.

    do n = start, end
      if (mod(n, 10000) == 0 .and. n > start) print*, counter, n, "step calculating..."
      ! call coef_pred
      cx = estimate(1, [aa(n), elv(n), ma(n)])
      cm = estimate(2, [aa(n), elv(n), ma(n)])
      cz = estimate(3, [aa(n), elv(n), ma(n)])
      cxh(n) = cx
      cmh(n) = cm
      czh(n) = cz
      call rk4

      call rot_iner
      call check_conv(n, end_flag)
      if (end_flag) exit
      n_last = n
    end do
  end subroutine flight

end program main
!=======================================================================

! 初期値設定
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
    read(55, *) tar_thtt
    read(55, *) file_out
    read(55, *) ini_elv
  close(55)

  elv_lower = -40.0
  elv_upper = 40.0
  elv_bound = 7.0

end subroutine ini_cond

!=======================================================================
subroutine cal_init
  use para
  implicit none

  aa(0)=deg*atan(ww(0)/uu(0))
  ma(0)=v_ini(0)/onsoku
  elv = 0
  elv(0) = ini_elv

  cx = 0.0
  cy = 0.0
  cz = 0.0
  cl = 0.0
  cm = 0.0
  cn = 0.0
  cxh = 0.0
  cmh = 0.0
  czh = 0.0

  counter = 0 ! 表示用

end subroutine cal_init

subroutine load_elv(file)
  use para
  implicit none
  character(*), intent(in) :: file
  real, allocatable :: elv_input(:)
  integer, allocatable :: steps(:)
  integer :: lines, unit, i

  ! elv = 0
  open(newunit=unit, file=file, status="old")
    read(unit, *) lines
    allocate(elv_input(lines), steps(lines))
    do i = 1, lines
      read(unit, *) steps(i), elv_input(i)
    end do
  close(unit)

  do i = 1, lines
    if (i < lines) then
      elv(steps(i):min(steps(i+1)-1,div)) = elv_input(i)
      if (steps(i+1) > div) exit
    else
      elv(steps(lines):) = elv_input(lines)
    end if
  end do

  elv(0) = ini_elv
end subroutine load_elv

subroutine load_elv2(file)
  use para
  implicit none
  character(*), intent(in) :: file
  real, allocatable :: elv_input(:)
  integer, allocatable :: steps(:)
  integer :: lines, unit, i, j

  open(newunit=unit, file=file, status="old")
    read(unit, *) lines
    allocate(elv_input(lines), steps(lines))
    do i = 1, lines
      read(unit, *) steps(i), elv_input(i)
    end do
  close(unit)

  do i = 1, lines
    if (i < lines) then
      do j = steps(i), min(steps(i+1)-1, div)
        if (j == 0) cycle
        elv(j) = min(max(elv(j - 1) + elv_input(i) * delt, elv_lower), elv_upper)
      end do
      if (steps(i+1) > div) exit
    else
      do j = steps(lines), div
        elv(j) = min(max(elv(j - 1) + elv_input(i) * delt, elv_lower), elv_upper)
      end do
    end if
  end do
end subroutine load_elv2

!=======================================================================
!     function of EOM ( 6 degrees of freedom )
!=======================================================================

!-------translational motion
real function fu(t, q, w, r, v, tht)
  use para
  implicit none
  real t, q, w, r, v, tht
  fu = (((-q*w/deg)+(r*v/deg))-(grav*sin(tht*rad)) &
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

  u = uu(n)
  v = vv(n)
  w = ww(n)

  p = pp(n)
  q = qq(n)
  r = rr(n)

  phi = phip(n)
  tht = thtt(n)
  psi = psip(n)

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

  ! elv(n+1) = 0

  !      aa(n+1)=0
  !      re(n+1)=0

  ! call calc_elv(n)

  ! open(201, file='inp', status='replace')
  ! write(201, *) aa(n+1), 'tes2'
  ! write(201, *) elv(n+1)
  ! write(201, *) ma(n+1), 'tes2'
  ! close (201)
end subroutine rk4

subroutine calc_elv(i)
  use para
  implicit none
  integer, intent(in) :: i

  elv(i) = uu(i)
end subroutine calc_elv

subroutine check_conv(i, flag)
  use para
  implicit none
  integer, intent(in) :: i       ! => the number of calculated steps
  logical, intent(out) :: flag   ! => TRUE then calculation finish

  if (zz(i) <= 0) then
    flag = .true.
    print *, "zz= ", zz(i)
  ! else if (aa(i) < -5 .or. aa(i) > 10) then
  !   flag = .true.
  !   print *, "aa= ", aa(i)
  ! else if (ma(i) < 0.1 .or. ma(i) > 0.5) then
  !   flag = .true.
  !   print *, "ma= ", ma(i)
  ! else if (thtt(i) < -90 .or. thtt(i) > 90) then
  !   flag = .true.
  !   print *, "thtt= ", thtt(i)
  else
    flag = .false.
  end if
end subroutine check_conv

!======================================================================

subroutine rot_iner
  use para
  implicit none

  ! do n = 0, n_last-1
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



    !======================================================================
    !for trajectory

    xx(n+1) = xx(n) + ue(n)*delt
    yy(n+1) = yy(n) + ve(n)*delt
    zz(n+1) = zz(n) + we(n)*delt
  ! end do
end subroutine rot_iner

!======================================================================

subroutine out_total(start, end)
  use para
  implicit none
  integer, intent(in) :: start, end
  character(:), allocatable :: format
  integer :: unit

  format = "(17(es15.5',')es15.5)"

  if (start == 0) then
    open(newunit=unit, file=file_out, status='replace')
    write(unit, "(a)") "time,theta,alpha,elv,ma,x,y,z,ue,ve,we,uu,vv,ww,cx,cm,cz,q"
  else
    open(newunit=unit, file=file_out, position='append')
  end if

  do n = start, min(end, n_last), 50
    write(unit, format) tt(n), thtt(n), aa(n), elv(n), ma(n), xx(n), yy(n), zz(n), &
                        ue(n), ve(n), we(n), ue(n), vv(n), ww(n), cxh(n), cmh(n), czh(n),qq(n)
  end do
  close(unit)
end subroutine out_total

subroutine chech_kriging
  use dof_kriging
  implicit none

  real :: a, e, m, y1, y2, y3
  integer :: unit1, unit2, unit3, i, j

  ! call init_krig(num_var=3, ga_pop_size=50)
  ! call load_krig(num_var=3)

  open(newunit=unit1, file="test1.csv", status="replace")
  open(newunit=unit2, file="test2.csv", status="replace")
  open(newunit=unit3, file="test3.csv", status="replace")
    do i = 0, 100
      if (mod(i, 10) == 0) print *, i
      do j = 0, 100
        a = i * 0.15 - 5
         ! a = 1.5
        e = j * 0.8 - 40
        ! e = 5
        ! m = j * 0.004 + 0.1
        m = 0.4
        y1 = estimate(1, [a, e, m])
        y2 = estimate(2, [a, e, m])
        y3 = estimate(3, [a, e, m])
        write(unit1, "(3(f0.10','))") a, e, y1
        write(unit2, "(3(f0.10','))") a, e, y2
        write(unit3, "(3(f0.10','))") a, e, y3
      end do
    end do
  close(unit3)
  close(unit2)
  close(unit1)

  ! open(newunit=unit, file="test2.csv", status="replace")
  !   do i = 1, 90
  !     if (mod(i, 10) == 0) print *, i
  !     a = -1
  !     e = i - 40
  !     y = estimate2(1, a, e, 0.3)
  !     write(unit, "(i0, 2(','f0.10))") i, e, y
  !   end do
  ! close(unit)

  stop
end subroutine chech_kriging
