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
  real cx_g, cm_g, cz_g
  real x_thrust
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
  real, allocatable :: x_thrust_ar(:)
end module para

module dof_rk4
  use para
  implicit none

  contains

  !-------translational motion
  real function fu(t, q, w, r, v, tht)
    implicit none
    real t, q, w, r, v, tht
    fu = (((-q*w/deg)+(r*v/deg))-(grav*sin(tht*rad)) &
    & +((rho*(v_ini(n)**2)*s_ref*cx)/(2*m_body)))+x_thrust
  end function fu

  real function fv(t, r, u, p, w, tht, phi)
    implicit none
    real t, r, u, p, w, tht, phi
    fv = ((-(r*u/deg)+(p*w/deg))+(grav*cos(tht*rad)*sin(phi*rad)) &
    & +((rho*(v_ini(n)**2)*s_ref*cy)/(2*m_body)))
  end function fv

  real function fw(t, p, v, q, u, tht, phi)
    implicit none
    real t, p, v, q, u, tht, phi
    fw = ((-(p*v/deg)+(q*u/deg))+(grav*cos(tht*rad)*cos(phi*rad)) &
    & +((rho*(v_ini(n)**2)*s_ref*cz)/(2*m_body)))
  end function fw
  !-------angular acceleration
  real function fp(t, p, q, r)
    implicit none
    real t, p, q, r
    roll = ((iy-iz)*q*r/deg)+(ixz*p*q/deg) &
    & +((rho*(v_ini(n)**2)*s_ref*b_span*deg*cl)/2)
    yaw = ((ix-iy)*p*q/deg)-(ixz*q*r/deg) &
    & +((rho*(v_ini(n)**2)*s_ref*b_span*deg*cn)/2)
    fp = (((roll/ix)+(ixz/ix)*(yaw/iz))/(1.d0-((ixz**2)/(ix*iz))))
  end function fp

  real function fq(t, p, r)
    implicit none
    real t, r, p
    pitch = (((iz-ix)*r*p/deg)+((ixz*(r**2-p**2))/deg)) &
    & +((rho*(v_ini(n)**2)*s_ref*mac*deg*cm)/2)
    fq = (pitch/iy)
  end function fq

  real function fr(t, p, q, r)
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
    implicit none
    real t, r, q, phi, tht
    fpsi = ((r*cos(phi*rad)+q*sin(phi*rad))/(cos(tht*rad)))
  end function fpsi

  real function ftht(t, q, r, phi, tht)
    implicit none
    real t, r, q, phi, tht
    ftht = (q*cos(phi*rad)-r*sin(phi*rad))
  end function ftht

  real function fphi(t, p, q, r, phi, tht)
    implicit none
    real t, p, q, r, phi, tht
    fphi =(p+(r*cos(phi*rad)+q*sin(phi*rad))*tan(tht*rad))
  end function fphi

  subroutine rk4
    implicit none

    real t!, fu, fv, fw, fp, fq, fr, fphi, ftht, fpsi
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
  end subroutine rk4

  subroutine rot_iner
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

  subroutine calc_x_thrust
    x_thrust = 
  end subroutine calc_x_thrust
end module dof_rk4
