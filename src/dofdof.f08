!=======================================================================
!         6DoF EOM solver
!=======================================================================

module global
    implicit none

    ! 定数
    real(8), parameter :: pi = 4 * atan(1d0)
    real(8), parameter :: radians = pi / 180d0
    real(8), parameter :: degrees = 180d0 / pi

    integer :: n_step
    real(8) :: time_end, delta_t

    real(8) :: grav, s_ref, m_body, rho, mac, b_span
    real(8) :: cx, cy, cz, cl, cm, cn
    real(8) :: ix, iy, iz, ixy, ixz, iyz

    real(8) :: sonic_speed

    real(8) x_start, y_start, z_start
    real(8) u_start, v_start, w_start
    real(8) p_start, q_start, r_start
    real(8) phi_start, tht_start, psi_start
    real(8) cx_g, cm_g, cz_g
    real(8) x_thrust

    real(8) ei
    real(8) cx_pred, cy_pred, cz_pred
    real(8) cl_pred, cm_pred, cn_pred
end module global


module equation
    use global
    implicit none

    contains

    !=======================================================================
    !     function of EOM ( 6 degrees of freedom )
    !=======================================================================

    !-------translational motion
    real(8) function fu(velocity, q, w, r, v, tht)
        real(8), intent(in) :: velocity, q, w, r, v, tht

        fu = -q * w + r * v - grav * sin(tht) &
           + rho * velocity ** 2 * s_ref * cx / (2 * m_body)
    end function fu

    real(8) function fv(velocity, r, u, p, w, tht, phi)
        real(8), intent(in) :: velocity, r, u, p, w, tht, phi

        fv = -r * u + p * w + grav * cos(tht) * sin(phi) &
           + rho * velocity ** 2 * s_ref * cy / (2 * m_body)
    end function fv


    real(8) function fw(velocity, p, v, q, u, tht, phi)
        real(8), intent(in) :: velocity, p, v, q, u, tht, phi

        fw = -p * v + q * u + grav * cos(tht) * cos(phi) &
           + rho * velocity ** 2 * s_ref * cz / (2 * m_body)
    end function fw


    !-------angular acceleration
    real(8) function fp(velocity, p, q, r)
        real(8), intent(in) :: velocity, p, q, r
        real(8) :: roll, yaw

        roll = ((iy - iz) * q * r + ixz * p * q) &
             + rho * velocity ** 2 * s_ref * b_span * cl / 2
        yaw = ((ix - iy) * p * q - ixz * q * r) &
            + rho * velocity ** 2 * s_ref * b_span * cn / 2
        fp = (roll / ix + ixz / ix * yaw / iz) / (1 - ixz ** 2 / (ix * iz))
    end function fp


    real(8) function fq(velocity, p, r)
        real(8), intent(in) :: velocity, p, r
        real(8) :: pitch

        pitch = (iz - ix) * r * p + ixz * (r ** 2 - p ** 2) &
              + rho * velocity ** 2 * s_ref * mac * cm / 2
        fq = pitch / iy
    end function fq


    real(8) function fr(velocity, p, q, r)
        real(8), intent(in) :: velocity, p, q, r
        real(8) :: roll, yaw

        roll = (iy - iz) * q * r + ixz * p * q &
             + rho * velocity ** 2 * s_ref * b_span * cl / 2
        yaw = (ix - iy) * p * q - ixz * q * r &
            + rho * velocity ** 2 * s_ref * b_span * cn / 2
        fr = yaw / iz + ixz / iz * roll / ix / (1 - ixz ** 2 / (ix * iz))
    end function fr


    !-------attitude angle
    real(8) function fpsi(q, r, phi, tht)
        real(8), intent(in) :: q, r, phi, tht

        fpsi = (r * cos(phi) + q * sin(phi)) / cos(tht)
    end function fpsi


    real(8) function ftht(q, r, phi, tht)
        real(8), intent(in) :: q, r, phi, tht

        ftht = q * cos(phi) - r * sin(phi)
    end function ftht


    real(8) function fphi(p, q, r, phi, tht)
        real(8), intent(in) :: p, q, r, phi, tht

        fphi = (p + r * cos(phi) + q * sin(phi)) * tan(tht)
    end function fphi


    subroutine rk4(u, v, w, p, q, r, phi, tht, psi,     &
        u_next, v_next, w_next, p_next, q_next, r_next, &
        phi_next, tht_next, psi_next)

        real(8), intent(in) :: u, v, w, p, q, r, phi, tht, psi
        real(8), intent(out) :: u_next, v_next, w_next
        real(8), intent(out) :: p_next, q_next, r_next
        real(8), intent(out) :: phi_next, tht_next, psi_next

        real(8) :: k1_u, k1_v, k1_w
        real(8) :: k2_u, k2_v, k2_w
        real(8) :: k3_u, k3_v, k3_w
        real(8) :: k4_u, k4_v, k4_w

        real(8) :: k1_p, k1_q, k1_r
        real(8) :: k2_p, k2_q, k2_r
        real(8) :: k3_p, k3_q, k3_r
        real(8) :: k4_p, k4_q, k4_r

        real(8) :: k1_phi, k1_tht, k1_psi
        real(8) :: k2_phi, k2_tht, k2_psi
        real(8) :: k3_phi, k3_tht, k3_psi
        real(8) :: k4_phi, k4_tht, k4_psi

        real(8) :: velocity

        velocity = sqrt(u ** 2 + v ** 2 + w ** 2)

        !-------------------k1
        k1_u = delta_t * fu(velocity, q, w, r, v, tht)
        k1_v = delta_t * fv(velocity, r, u, p, w, tht, phi)
        k1_w = delta_t * fw(velocity, p, v, q, u, tht, phi)

        k1_p = delta_t * fp(velocity, p, q, r)
        k1_q = delta_t * fq(velocity, p, r)
        k1_r = delta_t * fr(velocity, p, q, r)

        k1_psi = delta_t * fpsi(q, r, phi, tht)
        k1_tht = delta_t * ftht(q, r, phi, tht)
        k1_phi = delta_t * fphi(p, q, r, phi, tht)
        !-------------------k2
        k2_u = delta_t * fu(velocity, q+k1_q/2, w+k1_w/2, r+k1_r/2, v+k1_v/2, tht+k1_tht/2)
        k2_v = delta_t * fv(velocity, r+k1_r/2, u+k1_u/2, p+k1_p/2, w+k1_w/2, tht+k1_tht/2, phi+k1_phi/2)
        k2_w = delta_t * fw(velocity, p+k1_p/2, v+k1_v/2, q+k1_q/2, u+k1_u/2, tht+k1_tht/2, phi+k1_phi/2)

        k2_p = delta_t * fp(velocity, p+k1_p/2, q+k1_q/2, r+k1_r/2)
        k2_q = delta_t * fq(velocity, p+k1_p/2, r+k1_r/2)
        k2_r = delta_t * fr(velocity, p+k1_p/2, q+k1_q/2, r+k1_r/2)

        k2_psi = delta_t * fpsi(q+k1_q/2, r+k1_r/2, phi+k1_phi/2, tht+k1_tht/2)
        k2_tht = delta_t * ftht(q+k1_q/2, r+k1_r/2, phi+k1_phi/2, tht+k1_tht/2)
        k2_phi = delta_t * fphi(p+k1_p/2, q+k1_q/2, r+k1_r/2, phi+k1_phi/2, tht+k1_tht/2)
        !-------------------k3
        k3_u = delta_t * fu(velocity, q+k2_q/2, w+k2_w/2, r+k2_r/2, v+k2_v/2, tht+k2_tht/2)
        k3_v = delta_t * fv(velocity, r+k2_r/2, u+k2_u/2, p+k2_p/2, w+k2_w/2, tht+k2_tht/2, phi+k2_phi/2)
        k3_w = delta_t * fw(velocity, p+k2_p/2, v+k2_v/2, q+k2_q/2, u+k2_u/2, tht+k2_tht/2, phi+k2_phi/2)

        k3_p = delta_t * fp(velocity, p+k2_p/2, q+k2_q/2, r+k2_r/2)
        k3_q = delta_t * fq(velocity, p+k2_p/2, r+k2_r/2)
        k3_r = delta_t * fr(velocity, p+k2_p/2, q+k2_q/2, r+k2_r/2)

        k3_psi = delta_t * fpsi(q+k2_q/2, r+k2_r/2, phi+k2_phi/2, tht+k2_tht/2)
        k3_tht = delta_t * ftht(q+k2_q/2, r+k2_r/2, phi+k2_phi/2, tht+k2_tht/2)
        k3_phi = delta_t * fphi(p+k2_p/2, q+k2_q/2, r+k2_r/2, phi+k2_phi/2, tht+k2_tht/2)
        !-------------------k4
        k4_u = delta_t * fu(velocity, q+k3_q, w+k3_w, r+k3_r, v+k3_v, tht+k3_tht)
        k4_v = delta_t * fv(velocity, r+k3_r, u+k3_u, p+k3_p, w+k3_w, tht+k3_tht, phi+k3_phi)
        k4_w = delta_t * fw(velocity, p+k3_p, v+k3_v, q+k3_q, u+k3_u, tht+k3_tht, phi+k3_phi)

        k4_p = delta_t * fp(velocity, p+k3_p, q+k3_q, r+k3_r)
        k4_q = delta_t * fq(velocity, p+k3_p, r+k3_r)
        k4_r = delta_t * fr(velocity, p+k3_p, q+k3_q, r+k3_r)

        k4_psi = delta_t * fpsi(q+k3_q, r+k3_r, phi+k3_phi, tht+k3_tht)
        k4_tht = delta_t * ftht(q+k3_q, r+k3_r, phi+k3_phi, tht+k3_tht)
        k4_phi = delta_t * fphi(p+k3_p, q+k3_q, r+k3_r, phi+k3_phi, tht+k3_tht)
        !-------------------output for eis-------------------------------------

        u_next = u + (k1_u+(2*k2_u)+(2*k3_u)+k4_u)/6d0
        v_next = v + (k1_v+(2*k2_v)+(2*k3_v)+k4_v)/6d0
        w_next = w + (k1_w+(2*k2_w)+(2*k3_w)+k4_w)/6d0

        p_next = p + (k1_p+(2*k2_p)+(2*k3_p)+k4_p)/6d0
        q_next = q + (k1_q+(2*k2_q)+(2*k3_q)+k4_q)/6d0
        r_next = r + (k1_r+(2*k2_r)+(2*k3_r)+k4_r)/6d0

        phi_next = phi + (k1_phi+(2*k2_phi)+(2*k3_phi)+k4_phi)/6d0
        tht_next = tht + (k1_tht+(2*k2_tht)+(2*k3_tht)+k4_tht)/6d0
        psi_next = psi + (k1_psi+(2*k2_psi)+(2*k3_psi)+k4_psi)/6d0
    end subroutine rk4

end module equation


module flight
    use global
    implicit none

    type :: aircraft
        real(8) :: x, y, z ! 位置
        real(8) :: u, v, w ! 速度
        real(8) :: p, q, r ! 角速度
        real(8) :: phi, theta, psi ! 角度

        real(8) :: u_earth, v_earth, w_earth
        real(8) :: velocity
        real(8) :: mach_number
        real(8) :: angle_of_attack

        real(8) :: elevator
        real(8) :: cx, cm, cz

        contains

        procedure :: set_state
    end type aircraft

    contains

    subroutine set_state(self)
        class(aircraft), intent(inout) :: self
        real(8) :: velocity, angle_of_attack, mach_number
        real(8) :: u_earth, v_earth, w_earth
        real(8) :: tht, phi, psi

        phi = self%phi
        tht = self%theta
        psi = self%psi

        velocity        = sqrt(self%u ** 2 + self%v ** 2 + self%w ** 2)
        angle_of_attack = degrees * atan(self%w / self%u)
        mach_number     = velocity / sonic_speed

        u_earth = cos(tht) * cos(psi) * self%u   &
                + (sin(phi) * sin(tht) * cos(psi)  &
                   - cos(phi) * sin(psi)) * self%v &
                + (cos(phi) * sin(tht) * cos(psi)  &
                   + sin(phi) * sin(psi)) * self%w

        v_earth = cos(tht) * sin(psi) * self%u   &
                + (sin(phi) * sin(tht) * sin(psi)  &
                   + cos(phi) * cos(psi)) * self%v &
                + (cos(phi) * sin(tht) * sin(psi)  &
                   - sin(phi) * cos(psi)) * self%w

        w_earth = sin(tht) * self%u              &
                - sin(phi) * cos(tht) * self%v  &
                - cos(phi) * cos(tht) * self%w

        self%u_earth = u_earth
        self%v_earth = v_earth
        self%w_earth = w_earth
        self%velocity = velocity
        self%mach_number = mach_number
        self%angle_of_attack = angle_of_attack
    end subroutine set_state

end module flight


! program test
!     use global
!     use flight
!     use equation
!     use dof_kriging
!     implicit none

!     type(aircraft), allocatable :: ac(:)
!     integer :: i

!     call init

!     print *, 1, ac(1)%x, ac(1)%z
!     do i = 2, n_step

!         cx = estimate(1, ac(i)%angle_of_attack, ac(i)%elevator, ac(i)%mach_number)
!         cm = estimate(2, ac(i)%angle_of_attack, ac(i)%elevator, ac(i)%mach_number)
!         cz = estimate(3, ac(i)%angle_of_attack, ac(i)%elevator, ac(i)%mach_number)

!         call calc_motion(ac(i-1), ac(i))
!         call ac(i)%set_state
!         print *, i, ac(i)%x, ac(i)%z
!     end do

!     contains

!     ! ======================================================
!     ! 初期化
!     ! ======================================================

!     subroutine init
!         call read_input
!         n_step = time_end / delta_t
!         allocate(ac(n_step))

!         call load_krig(num_var=3)

!         call init_aircraft(ac(1))
!         call ac(1)%set_state
!     end subroutine init


!     subroutine read_input
!         integer :: unit

!         open(newunit=unit, file='initial_cond.inp', status="old")
!         read(unit, *)
!         read(unit, *) rho, grav
!         read(unit, *) sonic_speed
!         read(unit, *) m_body
!         read(unit, *) s_ref
!         read(unit, *) mac, b_span
!         read(unit, *) ix, iy, iz
!         read(unit, *) ixy, ixz
!         read(unit, *) iyz
!         read(unit, *) time_end, delta_t
!         read(unit, *) x_start, y_start, z_start
!         read(unit, *) u_start, v_start, w_start
!         read(unit, *) p_start, q_start, r_start
!         read(unit, *) phi_start, tht_start, psi_start
!         close(unit)

!         u_start = radians * u_start
!         v_start = radians * v_start
!         w_start = radians * w_start
!         phi_start = radians * phi_start
!         tht_start = radians * tht_start
!         psi_start = radians * psi_start
!     end subroutine read_input


!     subroutine init_aircraft(ac)
!         type(aircraft), intent(inout) :: ac
!         ac%x = x_start
!         ac%y = y_start
!         ac%z = z_start
!         ac%u = u_start
!         ac%v = v_start
!         ac%w = w_start
!         ac%p = p_start
!         ac%q = q_start
!         ac%r = r_start
!         ac%phi = phi_start
!         ac%theta = tht_start
!         ac%psi = psi_start
!         ac%elevator = 0
!     end subroutine init_aircraft


!     subroutine calc_motion(ac_in, ac_out)
!         type(aircraft), intent(in) :: ac_in
!         type(aircraft), intent(out) :: ac_out

!         call rk4(ac_in%u, ac_in%v, ac_in%w, ac_in%p, ac_in%q, ac_in%r, &
!                  ac_in%phi, ac_in%theta, ac_in%psi, &
!                  ac_out%u, ac_out%v, ac_out%w, ac_out%p, ac_out%q, ac_out%r, &
!                  ac_out%phi, ac_out%theta, ac_out%psi)

!         call ac_out%set_state

!         ac_out%x = ac_in%x + ac_out%u_earth * delta_t
!         ac_out%y = ac_in%y + ac_out%v_earth * delta_t
!         ac_out%z = ac_in%z + ac_out%w_earth * delta_t
!     end subroutine calc_motion

! end program test


module dof
    use iso_c_binding
    use global
    use flight
    use equation
    use dof_kriging
    implicit none

    type(aircraft), allocatable :: ac(:)

    contains

    ! ======================================================
    ! 初期化
    ! ======================================================

    subroutine read_input
        integer :: unit

        open(newunit=unit, file='initial_cond.inp', status="old")
        read(unit, *)
        read(unit, *) rho, grav
        read(unit, *) sonic_speed
        read(unit, *) m_body
        read(unit, *) s_ref
        read(unit, *) mac, b_span
        read(unit, *) ix, iy, iz
        read(unit, *) ixy, ixz
        read(unit, *) iyz
        read(unit, *) time_end, delta_t
        read(unit, *) x_start, y_start, z_start
        read(unit, *) u_start, v_start, w_start
        read(unit, *) p_start, q_start, r_start
        read(unit, *) phi_start, tht_start, psi_start
        close(unit)

        p_start = radians * p_start
        q_start = radians * q_start
        r_start = radians * r_start
        phi_start = radians * phi_start
        tht_start = radians * tht_start
        psi_start = radians * psi_start
    end subroutine read_input


    subroutine init_aircraft(ac)
        type(aircraft), intent(inout) :: ac
        ac%x = x_start
        ac%y = y_start
        ac%z = z_start
        ac%u = u_start
        ac%v = v_start
        ac%w = w_start
        ac%p = p_start
        ac%q = q_start
        ac%r = r_start
        ac%phi = phi_start
        ac%theta = tht_start
        ac%psi = psi_start
        ac%elevator = 0

        call ac%set_state
    end subroutine init_aircraft


    subroutine calc_motion(ac_in, ac_out)
        type(aircraft), intent(in) :: ac_in
        type(aircraft), intent(out) :: ac_out

        call rk4(ac_in%u, ac_in%v, ac_in%w, ac_in%p, ac_in%q, ac_in%r, &
                 ac_in%phi, ac_in%theta, ac_in%psi, &
                 ac_out%u, ac_out%v, ac_out%w, ac_out%p, ac_out%q, ac_out%r, &
                 ac_out%phi, ac_out%theta, ac_out%psi)

        call ac_out%set_state

        ac_out%x = ac_in%x + ac_out%u_earth * delta_t
        ac_out%y = ac_in%y + ac_out%v_earth * delta_t
        ac_out%z = ac_in%z + ac_out%w_earth * delta_t
    end subroutine calc_motion


    ! ======================================================
    ! インターフェイス
    ! ======================================================

    subroutine init() bind(c, name="init")
        call read_input
        n_step = time_end / delta_t

        ! call init_krig(num_var=3, pop_size=100, steps=100, scale=[100d0, 100d0, 100d0])
        call load_krig(num_var=3)
    end subroutine init


    subroutine calc(elevator_in, n, d, output) bind(c, name="calc")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in), value :: d
        real(c_double), intent(in) :: elevator_in(n)
        real(c_double), intent(out) :: output(6, n)
        real(8) :: angle_of_attack, elevator, mach_number
        integer :: i

        n_step = n
        delta_t = d

        if (allocated(ac)) then
            if (size(ac) /= n_step) then
                deallocate(ac)
                allocate(ac(n_step))
            end if
        else
            allocate(ac(n_step))
        end if

        call init_aircraft(ac(1))
        ! ac(1)%elevator = elevator_in(1)

        ! print *, 1, ac(1)%x, ac(1)%z
        ! output(1, 1) = ac(1)%x
        ! output(2, 1) = ac(1)%z
        ! output(3, 1) = ac(1)%u
        ! output(4, 1) = ac(1)%mach_number
        ! output(5, 1) = ac(1)%angle_of_attack
        ! output(6, 1) = ac(1)%theta * degrees

        do i = 1, n_step
            ac(i)%elevator = elevator_in(i)

            angle_of_attack = ac(i)%angle_of_attack
            elevator = ac(i)%elevator
            mach_number = ac(i)%mach_number

            cx = estimate(1, angle_of_attack, elevator, mach_number) ! * degrees !?
            cm = estimate(2, angle_of_attack, elevator, mach_number) ! * degrees !?
            cz = estimate(3, angle_of_attack, elevator, mach_number) ! * degrees !?

            ac(i)%cx = cx
            ac(i)%cm = cm
            ac(i)%cz = cz

            print *, cm

            if (i < n_step) then
                call calc_motion(ac(i), ac(i+1))
                call ac(i+1)%set_state
            end if

            ! print *, i, ac(i)%x, ac(i)%z
            output(1, i) = ac(i)%x
            output(2, i) = ac(i)%z
            output(3, i) = ac(i)%u
            output(4, i) = ac(i)%mach_number
            output(5, i) = ac(i)%angle_of_attack
            output(6, i) = ac(i)%cx
        end do
    end subroutine calc
end module dof
