module dof_flight
  use dof_kriging
  use para
  use dof_rk4
  implicit none

  contains

  ! ============================================================================
  ! ファイル入力
  ! ============================================================================

  subroutine read_input
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
      read(55, *) cx_g, cm_g, cz_g
      read(55, *) x_thrust
      read(55, *) file_out
    close(55)
    print *, 'rho =', rho
  end subroutine read_input

  ! ============================================================================
  ! ファイル出力
  ! ============================================================================

  subroutine out_total(start, end)
    implicit none
    integer, intent(in) :: start, end
    character(:), allocatable :: format
    integer :: unit

    format = "(17(es15.8',')es15.8)"

    if (start == 0) then
      open(newunit=unit, file=file_out, status='replace')
      write(unit, "(a)") "time,theta,alpha,elv,ma,x,y,z,ue,ve,we,uu,vv,ww,cx,cm,cz,q"
    else
      open(newunit=unit, file=file_out, position='append')
    end if

    do n = start, min(end, n_last)
      write(unit, format) tt(n), thtt(n), aa(n), elv(n), ma(n), xx(n), yy(n), zz(n), &
                          ue(n), ve(n), we(n), uu(n), vv(n), ww(n), cxh(n), cmh(n), czh(n), qq(n)
    end do
    close(unit)
  end subroutine out_total

  ! ============================================================================
  ! 変数初期化・入力
  ! ============================================================================

  subroutine cal_init
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

  ! ============================================================================
  ! 運動計算ループ
  ! ============================================================================

  subroutine flight(start, end)
    implicit none
    integer, intent(in) :: start, end

    end_flag = .false.

    do n = start, end
      if (mod(n, 10000) == 0 .and. n > start) print*, counter, n, "step calculating..."
      ! call coef_pred
      cx = estimate(1, aa(n), elv(n), ma(n))
      cm = estimate(2, aa(n), elv(n), ma(n))
      cz = estimate(3, aa(n), elv(n), ma(n))
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

  subroutine flight_const(start, end)
    implicit none
    integer, intent(in) :: start, end

    end_flag = .false.

    do n = start, end
      if (mod(n, 10000) == 0 .and. n > start) print*, counter, n, "step calculating..."
      ! call coef_pred
      cx = cx_g
      cm = cm_g
      cz = cz_g
      call rk4

      call rot_iner
      call check_conv(n, end_flag)
      if (end_flag) exit
      n_last = n
    end do
  end subroutine flight_const

  subroutine check_conv(i, flag)
    implicit none
    integer, intent(in) :: i       ! => the number of calculated steps
    logical, intent(out) :: flag   ! => TRUE then calculation finish

    if (zz(i) <= 0) then
      flag = .true.
      print *, "zz= ", zz(i)
    else
      flag = .false.
    end if
  end subroutine check_conv
end module dof_flight
