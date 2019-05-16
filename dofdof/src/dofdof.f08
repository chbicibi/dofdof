!===============================================================================
!
!===============================================================================

module global
    implicit none

    ! 定数
    real(8), parameter :: PI = 4 * atan(1d0)
    real(8), parameter :: RADIANS = PI / 180d0
    real(8), parameter :: DEGREES = 180d0 / PI

    ! 運動計算用パラメータ
    real(8) :: dt                        ! 時間刻み幅

    ! 空力・運動計算用パラメータ
    real(8) :: grav                      ! 重力加速度
    real(8) :: rho                       ! 空気密度
    real(8) :: speed_of_sound            ! 音速

    real(8) :: s_ref                     ! 代表面積
    real(8) :: m_body                    ! 機体重量
    real(8) :: mac                       ! 機体MAC長
    real(8) :: b_span                    ! 翼スパン
    real(8) :: cx, cy, cz, cl, cm, cn    ! 空力係数
    real(8) :: ix, iy, iz, ixy, ixz, iyz ! 慣性モーメント
end module global


!===============================================================================
!
!===============================================================================

module mod_calculation
    implicit none

    contains

    subroutine runge_kutta(func, X, h, s, X_next)
        ! 陽的ルンゲ・クッタ法による常微分方程式 dX/dt = f(X) の計算
        ! 入力:
        !   func : 関数f(X) => K
        !   X(m) : 入力変数
        !   h : 刻み幅
        !   s : 次数
        ! 出力:
        !   X_next(m) : 次ステップの変数
        ! 変数:
        !   m : Xの次元数
        !   C : 係数ベクトル

        real(8), intent(in) :: X(:), h
        integer, intent(in) :: s
        real(8), intent(out) :: X_next(:)
        real(8) :: C(s), K(size(X), 0:s)
        integer :: i, m

        select case(s)
        case(1)              ! 1次 (オイラー法)
            C = [1]
        case(2)              ! 2次 (ホイン法)
            C = [1, 1]
        case(4)              ! 4次
            C = [1, 2, 2, 1]
        case default
            print *, "Error: Invalid DEGREES of the Runge–Kutta method."
            call exit(1)
        end select

        m = size(X)
        K = 0

        do i = 1, s
            call func(X + K(:, i-1) * h / C(i), K(:, i))
        end do

        X_next = X + sum(K(:, 1:) * spread(C, 1, m), dim=2) * h / sum(C)
    end subroutine runge_kutta

end module mod_calculation


!===============================================================================
!
!===============================================================================

module mod_aircraft
    use global
    use mod_calculation
    implicit none

    type :: aircraft
        real(8) :: x, y, z       ! 位置
        real(8) :: u, v, w       ! 速度
        real(8) :: p, q, r       ! 角速度
        real(8) :: phi, tht, psi ! 姿勢角

        real(8) :: u_earth, v_earth, w_earth
        real(8) :: velocity
        real(8) :: mach_number
        real(8) :: angle_of_attack

        real(8) :: t             ! 時刻
        real(8) :: elevator
        real(8) :: cx, cm, cz

        contains

        procedure :: set_state
    end type aircraft

    contains

    subroutine set_state(self)
        class(aircraft), intent(inout) :: self
        real(8) :: u, v, w, phi, tht, psi
        real(8) :: u_earth, v_earth, w_earth
        real(8) :: velocity, angle_of_attack, mach_number

        u = self%u
        v = self%v
        w = self%w
        phi = self%phi
        tht = self%tht
        psi = self%psi

        velocity        = sqrt(u ** 2 + v ** 2 + w ** 2)
        angle_of_attack = DEGREES * atan(w / u)
        mach_number     = velocity / speed_of_sound

        u_earth = cos(tht) * cos(psi)                                    * u &
                + (sin(phi) * sin(tht) * cos(psi) - cos(phi) * sin(psi)) * v &
                + (cos(phi) * sin(tht) * cos(psi) + sin(phi) * sin(psi)) * w

        v_earth = cos(tht) * sin(psi)                                    * u &
                + (sin(phi) * sin(tht) * sin(psi) + cos(phi) * cos(psi)) * v &
                + (cos(phi) * sin(tht) * sin(psi) - sin(phi) * cos(psi)) * w

        w_earth = sin(tht)            * u &
                - sin(phi) * cos(tht) * v &
                - cos(phi) * cos(tht) * w

        self%u_earth = u_earth
        self%v_earth = v_earth
        self%w_earth = w_earth
        self%velocity = velocity
        self%mach_number = mach_number
        self%angle_of_attack = angle_of_attack
    end subroutine set_state


    ! ================================================================
    ! function of EOM ( 6 DEGREES of freedom )
    ! ================================================================

    subroutine func(X, K)
        real(8), intent(in) :: X(9)
        real(8), intent(out) :: K(9)
        real(8) :: u, v, w, p, q, r, phi, tht, psi
        real(8) :: fu, fv, fw, fp, fq, fr, fphi, ftht, fpsi
        real(8) :: vel_sq, coef1, coef2, roll, pitch, yaw

        u = X(1)
        v = X(2)
        w = X(3)
        p = X(4)
        q = X(5)
        r = X(6)
        phi = X(7)
        tht = X(8)
        psi = X(9)

        vel_sq = sum(X(1:3) ** 2)    ! == u^2 + v^2 + w^2 (== |V|^2)
        coef1 = rho * vel_sq * s_ref ! == ρ * |V|^2 * S
        coef2 = coef1 / (2 * m_body) ! == ρ * |V|^2 * S / 2m

        roll  = (iy - iz) * q * r + ixz * p * q             + coef1 * b_span * cl / 2
        pitch = (iz - ix) * r * p + ixz * (r ** 2 - p ** 2) + coef1 * mac    * cm / 2
        yaw   = (ix - iy) * p * q + ixz * q * r             + coef1 * b_span * cn / 2

        fu = -q * w + r * v - grav * sin(tht)            + coef2 * cx
        fv = -r * u + p * w + grav * cos(tht) * sin(phi) + coef2 * cy
        fw = -p * v + q * u + grav * cos(tht) * cos(phi) + coef2 * cz

        fp = (roll / ix + ixz / ix * yaw / iz) / (1 - ixz ** 2 / (ix * iz))
        fq = pitch / iy
        fr = (yaw / iz + ixz / iz * roll / ix) / (1 - ixz ** 2 / (ix * iz))

        fphi = (r * cos(phi) + q * sin(phi)) / cos(tht)
        ftht = q * cos(phi) - r * sin(phi)
        fpsi = p + fphi * sin(tht)

        K = [fu, fv, fw, fp, fq, fr, fphi, ftht, fpsi]
    end subroutine func


    subroutine calc_motion(ac, ac_next)
        type(aircraft), intent(in) :: ac
        type(aircraft), intent(out) :: ac_next
        real(8) :: X(9), X_next(9)

        X = [ac%u, ac%v, ac%w, ac%p, ac%q, ac%r, ac%phi, ac%tht, ac%psi]

        call runge_kutta(func, X, dt, 4, X_next)

        ac_next%u = X_next(1)
        ac_next%v = X_next(2)
        ac_next%w = X_next(3)
        ac_next%p = X_next(4)
        ac_next%q = X_next(5)
        ac_next%r = X_next(6)
        ac_next%phi = X_next(7)
        ac_next%tht = X_next(8)
        ac_next%psi = X_next(9)

        call ac_next%set_state

        ac_next%x = ac%x + ac_next%u_earth * dt
        ac_next%y = ac%y + ac_next%v_earth * dt
        ac_next%z = ac%z + ac_next%w_earth * dt
    end subroutine calc_motion

end module mod_aircraft


!===============================================================================
!
!===============================================================================

module dof
    use iso_c_binding
    use global
    use mod_aircraft
    use dof_kriging
    implicit none

    real(8) x_start, y_start, z_start
    real(8) u_start, v_start, w_start
    real(8) p_start, q_start, r_start
    real(8) phi_start, tht_start, psi_start

    contains

    ! ================================================================
    ! 初期化
    ! ================================================================

    subroutine read_input
        integer :: unit

        open(newunit=unit, file='conditions.txt', status="old")
        read(unit, *) grav
        read(unit, *) rho
        read(unit, *) speed_of_sound
        read(unit, *) m_body
        read(unit, *) s_ref
        read(unit, *) mac
        read(unit, *) b_span
        read(unit, *) ix, iy, iz
        read(unit, *) ixy, ixz, iyz
        close(unit)

        open(newunit=unit, file='input.txt', status="old")
        read(unit, *) x_start, y_start, z_start
        read(unit, *) u_start, v_start, w_start
        read(unit, *) p_start, q_start, r_start
        read(unit, *) phi_start, tht_start, psi_start
        close(unit)

        p_start = RADIANS * p_start
        q_start = RADIANS * q_start
        r_start = RADIANS * r_start
        phi_start = RADIANS * phi_start
        tht_start = RADIANS * tht_start
        psi_start = RADIANS * psi_start
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
        ac%tht = tht_start
        ac%psi = psi_start
        ac%elevator = 0

        call ac%set_state
    end subroutine init_aircraft


    logical function check_break(ac) result(res)
        type(aircraft), intent(inout) :: ac

        if (ac%angle_of_attack < -5 .or. ac%angle_of_attack > 10) then
            print "(a,f0.1,a)", 'Error: "Angle of Attack" is out of range! (', ac%t, "s)"
            res = .true.
        else if (ac%elevator < -40 .or. ac%elevator > 40) then
            print "(a,f0.1,a)", 'Error: "Elevator" is out of range! (', ac%t, "s)"
            res = .true.
        else if (ac%mach_number < 0.1d0 .or. ac%mach_number > 0.5d0) then
            print "(a,f0.1,a)", 'Error: "Mach Number" is out of range! (', ac%t, "s)"
            res = .true.
        else
            res = .false.
        end if
    end function check_break


    ! ================================================================
    ! インターフェイス
    ! ================================================================

    subroutine init() bind(c, name="init")
        call read_input

        ! call init_krig(n_var=3, pop_size=100, steps=30, scale=[100d0, 50d0, 100d0])
        call load_krig(n_var=3)
    end subroutine init


    subroutine calc(elevator_in, n, d, output) bind(c, name="calc")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in), value :: d
        real(c_double), intent(in) :: elevator_in(n)
        real(c_double), intent(out) :: output(7, n)
        type(aircraft), allocatable :: ac(:)
        real(8) :: angle_of_attack, elevator, mach_number
        integer :: i

        dt = d
        output = transfer(z'7ff8000000000000', 0d0) ! == NaN

        if (allocated(ac)) then
            if (size(ac) /= n) then
                deallocate(ac)
                allocate(ac(n))
            end if
        else
            allocate(ac(n))
        end if

        call init_aircraft(ac(1))

        do i = 1, n
            ac(i)%t = dt * (i - 1)
            ac(i)%elevator = elevator_in(i)

            angle_of_attack = ac(i)%angle_of_attack
            elevator = ac(i)%elevator
            mach_number = ac(i)%mach_number

            if (check_break(ac(i))) then
                exit
            end if

            ! angle_of_attack = min(max(angle_of_attack, -5d0), 10d0)
            ! elevator = min(max(elevator, -40d0), 40d0)
            ! mach_number = min(max(mach_number, 0.1d0), 0.5d0)

            cx = estimate(1, angle_of_attack, elevator, mach_number)
            cm = estimate(2, angle_of_attack, elevator, mach_number)
            cz = estimate(3, angle_of_attack, elevator, mach_number)

            ac(i)%cx = cx
            ac(i)%cm = cm
            ac(i)%cz = cz

            if (i < n) then
                call calc_motion(ac(i), ac(i+1))
                ! call ac(i+1)%set_state
            end if

            ! print *, i, ac(i)%x, ac(i)%z
            output(1, i) = ac(i)%x
            output(2, i) = ac(i)%z
            output(3, i) = ac(i)%cm
            output(4, i) = ac(i)%q * DEGREES
            output(5, i) = ac(i)%tht * DEGREES
            output(6, i) = ac(i)%angle_of_attack
            output(7, i) = ac(i)%mach_number
        end do
    end subroutine calc

end module dof
