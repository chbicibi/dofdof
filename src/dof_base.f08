!===============================================================================
! 変数管理用モジュール (主に定数と計算パラメータに使用)
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
    ! real(8) :: cx, cy, cz, cl, cm, cn    ! 空力係数
    real(8) :: ix, iy, iz, ixy, ixz, iyz ! 慣性モーメント

end module global


!===============================================================================
! 数値計算用モジュール (常微分方程式を解く)
!===============================================================================

module mod_calculation
    implicit none

    private
    public :: runge_kutta

    contains


    !===========================================================================
    ! 陽的ルンゲ・クッタ法による常微分方程式 dX/dt = f(X) の計算
    !===========================================================================

    subroutine runge_kutta(func, X, h, s, X_next)
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

        interface
            module subroutine func(X, K)
                real(8), intent(in) :: X(9)
                real(8), intent(out) :: K(9)
            end subroutine func
        end interface

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
        case default         ! それ以外はエラー (必要に応じて追加可能)
            print *, "Error: Invalid DEGREES of the Runge–Kutta method."
            call exit(1)
        end select

        m = size(X) ! 求めたい変数の数を代入
        K = 0       ! 初期化

        ! 方程式を解く
        do i = 1, s
            call func(X + K(:, i-1) * h / C(i), K(:, i))
        end do

        ! 結果を代入
        X_next = X + sum(K(:, 1:) * spread(C, 1, m), dim=2) * h / sum(C)
    end subroutine runge_kutta

end module mod_calculation


!===============================================================================
! 運動計算用モジュール (機体の位置と姿勢を計算する)
!===============================================================================

module mod_aircraft
    use global
    use mod_calculation
    implicit none

    private
    public :: aircraft, calc_motion

    !===========================================================================
    ! 機体変数管理用の型 (運動計算に必要な変数をまとめて管理する)
    !===========================================================================

    type :: aircraft
        real(8) :: x, y, z         ! 位置
        real(8) :: u, v, w         ! 速度
        real(8) :: p, q, r         ! 角速度
        real(8) :: phi, tht, psi   ! 姿勢角

        real(8) :: u_earth         ! 速度 (u方向)
        real(8) :: v_earth         ! 速度 (v方向)
        real(8) :: w_earth         ! 速度 (w方向)
        real(8) :: velocity        ! 速度の絶対値
        real(8) :: mach_number     ! マッハ数
        real(8) :: angle_of_attack ! 迎角

        real(8) :: t               ! 時刻
        real(8) :: elevator        ! 舵角
        real(8) :: cx              ! 空力係数 (推力)
        real(8) :: cm              ! 空力係数 (ピッチ角)
        real(8) :: cz              ! 空力係数 (揚力)

        ! 2019.6.17追加
        real(8) :: thrust
        real(8) :: ang_thrust
        integer :: calculated      ! 0=>計算前, 1=>計算済み
        real(8) :: gamt

         ! 2020.6.11追加
        real(8) :: s_ref                     ! 代表面積
        real(8) :: m_body                    ! 機体重量
        real(8) :: mac                       ! 機体MAC長
        real(8) :: b_span                    ! 翼スパン
        real(8) :: ix, iy, iz, ixy, ixz, iyz ! 慣性モーメント

        ! 2020.06.16追加
        real(8) :: forth                     !　分離力
        real(8) :: moment                    !　分離モーメント

        contains

        procedure :: set_state
        procedure :: print_state
    end type aircraft

    ! real(8) :: ref_thrust, ref_ang_thrust ! 運動計算参照用変数 (角度はラジアン)
    ! real(8) :: ref_forth                  ! 分離力
    ! real(8) :: ref_moment                 ! 分離モーメント力


    contains


    !===========================================================================
    ! 運動計算ヘルパールーチン (計算結果から他の値を求めて機体変数管理用の変数に格納する)
    !===========================================================================

    subroutine set_state(this)
        class(aircraft), intent(inout) :: this
        real(8) :: u, v, w, phi, tht, psi
        real(8) :: u_earth, v_earth, w_earth
        real(8) :: velocity, angle_of_attack, mach_number, gamt

        ! ここの値は運動計算で求めてある (u, v, w, phi, tht, psi)
        u = this%u
        v = this%v
        w = this%w
        phi = this%phi
        tht = this%tht
        psi = this%psi

        ! 必要な値を計算
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

        gamt = DEGREES * atan(w_earth / u_earth)

        ! 結果を格納
        this%u_earth = u_earth
        this%v_earth = v_earth
        this%w_earth = w_earth
        this%velocity = velocity
        this%mach_number = mach_number
        this%angle_of_attack = angle_of_attack
        this%gamt = gamt
    end subroutine set_state

    !===========================================================================
    ! 相対位置計算サブルーチン
    !===========================================================================

    ! subroutine calc_position(ac1, ac2, rel_x, rel_z, del_tht)
    !     implicit none
    !     type(aircraft), intent(in) :: ac1, ac2
    !     real(8), intent(out) :: rel_x, rel_z, del_tht
    !     real(8) d_b,d_o
    !     real(8) xo, zo, xb, zb
    !     real(8) bunshi_z, bunbo_z, bunshi_x, bunbo_x

    !     ! rel_x = (ac2%x - ac1%x) / 40.8
    !     ! rel_z = (ac2%z - ac1%z) / 7.8
    !     del_tht = ac2%tht - ac1%tht

    !     d_b = sqrt(4.45d0 * 4.45d0 + 15.7d0 * 15.7d0)
    !     d_o = sqrt(3.8d0 * 3.8d0 + 19.8d0 * 19.8d0)

    !     xo = ac2%x - d_o * cos(ac2%tht + atan(3.8d0/19.8d0))
    !     zo = ac2%z - d_o * sin(ac2%tht + atan(3.8d0/19.8d0))

    !     xb = ac1%x - d_b * cos(abs(ac1%tht - atan(4.45d0/15.7d0)))
    !     zb = ac1%z - d_b * sin(ac1%tht - atan(4.45d0/15.7d0))


    !     bunshi_z = abs(tan(ac1%tht) * (xo - xb) + (zb - zo))
    !     bunbo_z =  8.9d0 * (sqrt(tan(ac1%tht) * tan(ac1%tht) + 1.0d0))

    !     bunshi_x = abs((xb - xo) + tan(ac1%tht) * (zb - zo))
    !     bunbo_x = 40.8d0 * (sqrt(tan(ac1%tht) * tan(ac1%tht) + 1.0d0))


    !     rel_z = bunshi_z / bunbo_z
    !     rel_x = bunshi_x / bunbo_x


    !     ! 相対位置計算式（重心）
    !     ! bunshi_z = abs(tan(ac1%tht) * (ac2%x - ac1%x) + (ac1%z - ac2%z))
    !     ! bunbo_z =  (sqrt(tan(ac1%tht) * tan(ac1%tht) + 1))

    !     ! bunshi_x = abs((ac1%x - ac2%x) + tan(ac1%tht) * (ac1%z - ac2%z))
    !     ! bunbo_x = (sqrt(tan(ac1%tht) * tan(ac1%tht) + 1))

    !     ! print *, "rel_x=", rel_x
    !     ! print *, "rel_z=", rel_z
    !     ! print *, "del_tht=", del_tht

    ! end subroutine calc_position


    ! ==========================================================================
    ! ルンゲクッタ法計算ルーチンに渡す関数
    ! ==========================================================================

    subroutine solve_differential_equations(X, cx, cy, cz, cl, cm, cn, thrust, &
                                            ang_thrust, X_next)
        real(8), intent(in) :: X(:), cx, cy, cz, cl, cm, cn, thrust, ang_thrust
        real(8), intent(out) :: X_next(:)

        call runge_kutta(diff_func, X, dt, 4, X_next)

        contains

        subroutine diff_func(X, K)
            real(8), intent(in) :: X(9)
            real(8), intent(out) :: K(9)
            real(8) :: u, v, w, p, q, r, phi, tht, psi
            real(8) :: fu, fv, fw, fp, fq, fr, fphi, ftht, fpsi
            real(8) :: vel_sq, coef1, coef2, roll, pitch, yaw

            ! 分かりやすいように配列から各変数に代入
            u = X(1)
            v = X(2)
            w = X(3)
            p = X(4)
            q = X(5)
            r = X(6)
            phi = X(7)
            tht = X(8)
            psi = X(9)

            ! 各関数を求める
            vel_sq = sum(X(1:3) ** 2)    ! [速度の二乗] == u^2 + v^2 + w^2 (== |V|^2)
            coef1 = rho * vel_sq * s_ref ! [係数1] == ρ * |V|^2 * S
            coef2 = coef1 / (2 * m_body) ! [係数2] == ρ * |V|^2 * S / 2m

            roll  = (iy - iz) * q * r + ixz * p * q             + coef1 * b_span * cl / 2
            pitch = (iz - ix) * r * p + ixz * (r ** 2 - p ** 2) + coef1 * mac    * cm / 2 ! + ref_moment
            yaw   = (ix - iy) * p * q + ixz * q * r             + coef1 * b_span * cn / 2

            ! 旧: 推力項無し↓
            ! fu = -q * w + r * v - grav * sin(tht)            + coef2 * cx
            ! fv = -r * u + p * w + grav * cos(tht) * sin(phi) + coef2 * cy
            ! fw = -p * v + q * u + grav * cos(tht) * cos(phi) + coef2 * cz

            ! 新: 推力項有り↓
            ! forth項あり
            ! fu = -q * w + r * v - grav * sin(tht)            + coef2 * cx + thrust / m_body * cos(ang_thrust) + forth / m_body * cos(PI/2.0d0 + tht)
            ! fv = -r * u + p * w + grav * cos(tht) * sin(phi) + coef2 * cy
            ! fw = -p * v + q * u + grav * cos(tht) * cos(phi) + coef2 * cz + thrust / m_body * sin(ang_thrust) + forth / m_body * sin(PI/2.0d0 + tht)

            fu = -q * w + r * v - grav * sin(tht)            + coef2 * cx + thrust / m_body * cos(ang_thrust)
            fv = -r * u + p * w + grav * cos(tht) * sin(phi) + coef2 * cy
            fw = -p * v + q * u + grav * cos(tht) * cos(phi) + coef2 * cz - thrust / m_body * sin(ang_thrust) ! - forth / m_body

            fp = (roll / ix + ixz / ix * yaw / iz) / (1 - ixz ** 2 / (ix * iz))
            fq = pitch / iy
            fr = (yaw / iz + ixz / iz * roll / ix) / (1 - ixz ** 2 / (ix * iz))

            fphi = (r * cos(phi) + q * sin(phi)) / cos(tht)
            ftht = q * cos(phi) - r * sin(phi)
            fpsi = p + fphi * sin(tht)

            ! 結果を配列に格納
            K = [fu, fv, fw, fp, fq, fr, fphi, ftht, fpsi]
        end subroutine diff_func
    end subroutine solve_differential_equations


    ! ==========================================================================
    ! 運動計算ルーチン
    ! ==========================================================================

    subroutine calc_motion(ac, ac_next)
        type(aircraft), intent(in) :: ac       ! 管理用変数 (入力)
        type(aircraft), intent(out) :: ac_next ! 管理用変数 (出力)
        real(8) :: X(9), X_next(9), cx, cm, cz, thrust, ang_thrust

        ! 計算しやすいように配列に格納
        X = [ac%u, ac%v, ac%w, ac%p, ac%q, ac%r, ac%phi, ac%tht, ac%psi]

        ! 計算に必要な空力係数等を一時変数に代入
        cx = ac%cx
        cm = ac%cm
        cz = ac%cz
        thrust = ac%thrust
        ang_thrust = ac%ang_thrust

        ! 2020.6.11 追加
        ! s_ref = ac%s_ref          ! 代表面積
        ! m_body = ac%m_body        ! 機体重量
        ! mac = ac%mac              ! 機体MAC長
        ! b_span = ac%b_span        ! 翼スパン
        ! ix = ac%ix
        ! iy = ac%iy
        ! iz = ac%iz
        ! ixy = ac%ixy
        ! ixz = ac%ixz
        ! iyz = ac%iyz              ! 慣性モーメント

        ! ルンゲクッタ法で次ステップの各値を求める
        ! 2019.6.17 推力項追加

        ! 旧: 推力項無し
        ! call runge_kutta(func_normal, X, dt, 4, X_next)

        ! 新: 推力項有り
        ! ref_thrust = ac%thrust         ! 参照用変数に代入
        ! ref_ang_thrust = ac%ang_thrust ! 参照用変数に代入

        ! 2020.06.17　分離力項あり
        ! ref_forth = ac%forth           ! 参照用変数に代入
        ! ref_moment = ac%moment         ! 参照用変数に代入

        call solve_differential_equations(X, cx, 0d0, cz, 0d0, cm, 0d0, thrust, ang_thrust, X_next)


        ! 結果を管理用変数に格納
        ac_next%u = X_next(1)
        ac_next%v = X_next(2)
        ac_next%w = X_next(3)
        ac_next%p = X_next(4)
        ac_next%q = X_next(5)
        ac_next%r = X_next(6)
        ac_next%phi = X_next(7)
        ac_next%tht = X_next(8)
        ac_next%psi = X_next(9)

        ! 追加で求めた値を管理用変数に格納
        call ac_next%set_state
        ac_next%calculated = 1

        ! 機体の座標を計算
        ac_next%x = ac%x + ac_next%u_earth * dt
        ac_next%y = ac%y + ac_next%v_earth * dt
        ac_next%z = ac%z + ac_next%w_earth * dt
    end subroutine calc_motion


    subroutine print_state(this)
        class(aircraft), intent(in) :: this
        print "(1(a9,es12.3))", "t:", this%t
        print "(3(a9,es12.3))", "x:", this%x, "y:", this%y, "z:", this%z
        print "(3(a9,es12.3))", "u:", this%u, "v:", this%v, "w:", this%w
        print "(3(a9,es12.3))", "p:", this%p, "q:", this%q, "r:", this%r
        print "(3(a9,es12.3))", "phi:", this%phi, "tht:", this%tht, "psi:", this%psi
        print "(3(a9,es12.3))", "ue:", this%u_earth, "ve:", this%v_earth, "we:", this%w_earth
        print "(3(a9,es12.3))", "cx:", this%cx, "cm:", this%cm, "cz:", this%cz
        print "(4(a9,es12.3))", "alpha:", this%angle_of_attack, "gamt:", this%gamt, "elv:", this%elevator, "ma:", this%mach_number
        print "(4(a9,es12.3))", "thrust:", this%thrust, "ang_thr:", this%ang_thrust, "forth:", this%forth, "moment:", this%moment
    end subroutine print_state

end module mod_aircraft
