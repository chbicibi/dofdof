!=======================================================================
program main     ! 6DoF main !
!=======================================================================
  use util
  use dof_kriging
  use para
  use dof_flight
  implicit none

  integer :: start_time

  call initialize

  ! ============================================================================
  ! main loop
  ! ============================================================================

  call system_clock(start_time) ! 時間計測開始

  call main_1

  call elapsed_time(start_time) ! 計測時間表示

  contains

  ! ============================================================================
  ! プログラム初期化
  ! ============================================================================

  subroutine initialize
    implicit none

    ! ******************************************************
    ! krigingモデル初期化
    ! ******************************************************
    ! if (.false.) then
    !   ! θ初期化用 普段は使わない
    !   call init_krig(num_var=3, pop_size=100, steps=100, scale=[100d0, 100d0, 100d0])
    !   stop
    ! else
    !   call load_krig(num_var=3)
    !   ! call chech_kriging
    ! end if

    ! ******************************************************
    ! 定数設定
    ! ******************************************************
    pi=4*atan(1.0) !definition of pi
    rad = pi/180        !convert [deg]¨[rad]
    deg = 180/pi        !convert [rad]¨[deg]

    ! ******************************************************
    ! ファイル入力
    ! ******************************************************
    call read_input

    div = tend / delt
    n_last = 0

    ! ******************************************************
    ! 配列割付
    ! ******************************************************
    allocate (tt(0:div))
    allocate (uu(0:div), vv(0:div), ww(0:div))
    allocate (v_ini(0:div), ma(0:div))
    allocate (ue(0:div), ve(0:div), we(0:div))
    allocate (pp(0:div), qq(0:div), rr(0:div))
    allocate (phip(0:div), thtt(0:div), psip(0:div), aa(0:div))
    allocate (xx(0:div), yy(0:div), zz(0:div))
    allocate (cxh(0:div), cmh(0:div), czh(0:div))
    allocate (elv(0:div))
    allocate (x_thrust_ar(0:div))


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
  end subroutine initialize

  subroutine main_1
    implicit none
    integer :: start, end

    start = 0
    end = div

    call flight_const(start, end)
    call out_total(start, end)
  end subroutine main_1
end program main
