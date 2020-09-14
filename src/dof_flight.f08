!===============================================================================
! 飛行計算モジュール
!===============================================================================
module mod_flight
    use mod_kriging
    use global
    use mod_aircraft
    use mod_wind
    implicit none

    private
    public :: read_input, init_aircraft, flight, flight_with_hook, output_path

    real(8) :: x_start, y_start, z_start
    real(8) :: u_start, v_start, w_start
    real(8) :: p_start, q_start, r_start
    real(8) :: phi_start, tht_start, psi_start

    ! 2020.9.14 追加
    real(8) :: angle_of_attack_lim(2)
    real(8) :: mach_number_lim(2)
    real(8) :: elevator_lim(2)

    contains

    ! ==========================================================================
    ! 初期化
    ! ==========================================================================

    ! [入力ファイル読み込み用サブルーチン]
    ! 入力ファイルの書式を変更する場合はここを書き換えてください
    subroutine read_input(file_conditions, file_input)
        character(*), intent(in) :: file_conditions, file_input
        integer :: unit

        ! 計算条件定義ファイル
        open(newunit=unit, file=file_conditions, status="old")
        read(unit, *) grav
        read(unit, *) rho
        read(unit, *) speed_of_sound
        read(unit, *) m_body
        read(unit, *) s_ref
        read(unit, *) mac
        read(unit, *) b_span
        read(unit, *) ix, iy, iz
        read(unit, *) ixy, ixz, iyz
        read(unit, *) angle_of_attack_lim
        read(unit, *) mach_number_lim
        read(unit, *) elevator_lim
        read(unit, *) wind_length, wind_power, wind_center
        close(unit)

        ! 初期条件定義ファイル
        open(newunit=unit, file=file_input, status="old")
        read(unit, *) x_start, y_start, z_start
        read(unit, *) u_start, v_start, w_start
        read(unit, *) p_start, q_start, r_start
        read(unit, *) phi_start, tht_start, psi_start
        close(unit)

        ! 単位をラジアンに変更
        ! [メモ] 2019.6.17
        ! このプログラムでは角度の単位をラジアンで管理しています (舵角を除く)
        ! 必要に応じて入出力時に変換してください
        p_start = RADIANS * p_start
        q_start = RADIANS * q_start
        r_start = RADIANS * r_start
        phi_start = RADIANS * phi_start
        tht_start = RADIANS * tht_start
        psi_start = RADIANS * psi_start

        is_calc_wind = wind_power > 0
    end subroutine read_input


    ! 機体の状態を初期化
    subroutine init_aircraft(ac, number)
        type(aircraft), intent(inout) :: ac
        character(1), intent(in) :: number
        character(:), allocatable :: file_conditions, file_input

        file_conditions = 'conditions' // number // '.txt'
        file_input = 'input' // number // '.txt'
        call read_input(file_conditions, file_input)

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
        ac%thrust = 0
        ac%ang_thrust = 0
        ac%calculated = 1

        ac%s_ref = s_ref
        ac%m_body = m_body
        ac%mac = mac
        ac%b_span = b_span
        ac%ix = ix
        ac%iy = iy
        ac%iz = iz
        ac%ixy = ixy
        ac%ixz = ixz
        ac%iyz = iyz

        ! 2020.06.16
        ! forth = ac%fotrh
        ac%forth = 0
        ac%moment = 0

        call ac%set_state
    end subroutine init_aircraft


    ! 軌道計算サブルーチン
    subroutine flight(acs, n_last, exit_status)
        type(aircraft), intent(inout) :: acs(:)
        integer, intent(out) :: n_last, exit_status
        integer :: total_steps, i

        total_steps = size(acs)

        do i = 1, total_steps - 1
            acs(i)%t = dt * (i - 1)
            ! if (i > 1) then
            !     acs(i)%elevator = acs(i-1)%elevator
            !     acs(i)%thrust = acs(i-1)%thrust
            !     acs(i)%ang_thrust = acs(i-1)%ang_thrust
            ! end if

            associate (param => [acs(i)%angle_of_attack, &
                                 acs(i)%elevator,        &
                                 acs(i)%mach_number])
                acs(i)%cx = estimate(1, param)
                acs(i)%cm = estimate(2, param)
                acs(i)%cz = estimate(3, param)
            end associate

            call calc_motion(acs(i), acs(i+1))

            n_last = i
            if (check_stop(acs(i+1), exit_status)) exit

            ! if (mod(i, 1000) == 0) print *, i, "step calculating..."
        end do
    end subroutine flight


    ! 軌道計算サブルーチン
    subroutine flight_with_hook(acs, n_last, exit_status, hook)
        interface
            module subroutine hook(ac, idx)
                type(aircraft), intent(in) :: ac
                integer, intent(in) :: idx
            end subroutine hook
        end interface

        type(aircraft), intent(inout) :: acs(:)
        integer, intent(out) :: n_last, exit_status
        integer :: total_steps, i

        total_steps = size(acs)

        do i = 1, total_steps - 1
            acs(i)%t = dt * (i - 1)

            associate (param => [acs(i)%angle_of_attack, &
                                 acs(i)%elevator,        &
                                 acs(i)%mach_number])
                acs(i)%cx = estimate(1, param)
                acs(i)%cm = estimate(2, param)
                acs(i)%cz = estimate(3, param)
            end associate

            call calc_motion(acs(i), acs(i+1))

            n_last = i
            if (check_stop(acs(i+1), exit_status)) exit

            call hook(acs(i), i)

            ! if (mod(i, 1000) == 0) print *, i, "step calculating..."
        end do
    end subroutine flight_with_hook


    ! 計算条件の確認
    logical function check_stop(ac, exit_status) result(is_stop)
        type(aircraft), intent(inout) :: ac
        integer, intent(out) :: exit_status

        if (ac%z <= 0) then
            is_stop = .true.
            exit_status = 1

        else if (ac%angle_of_attack < angle_of_attack_lim(1) .or. &
                 ac%angle_of_attack > angle_of_attack_lim(2)) then
            ! print "(a,f0.1,a)", 'Error: "Angle of Attack" is out of range! (', ac%t, "s)"
            is_stop = .true.
            exit_status = 2

        else if (ac%mach_number < mach_number_lim(1) .or. &
                 ac%mach_number > mach_number_lim(2)) then
            ! print "(a,f0.1,a)", 'Error: "Mach Number" is out of range! (', ac%t, "s)"
            is_stop = .true.
            exit_status = 3

        else if (ac%elevator < elevator_lim(1) .or. &
                 ac%elevator > elevator_lim(2)) then
            ! print "(a,f0.1,a)", 'Error: "Elevator" is out of range! (', ac%t, "s)"
            is_stop = .true.
            exit_status = 4

        else
            is_stop = .false.
            exit_status = 0
        end if

        ! 2020.9.14 追加
        ac%status = exit_status
    end function check_stop


    ! 飛行経路をファイルに出力する (csv形式)
    subroutine output_path(ac, start, end, step, file)
        type(aircraft), intent(in) :: ac(:)
        integer, intent(in) :: start, end, step
        character(*), intent(in) :: file

        character(:), allocatable :: format
        integer :: unit, i, n

        format = "(20(es15.8',')i0)" ! カンマ区切り

        if (start <= 1) then
            open(newunit=unit, file=file, status='replace')
            write(unit, "(a)") "time,theta,alpha,gamt,elv,ma,x,z,ue,we,uu,ww,q,cx,cm,cz,thrust,ang_thrust,forth,moment,status"
        else
            open(newunit=unit, file=file, position='append')
        end if

        do i = start, size(ac), step

            if (i == 0) then
                if (step > 1) then
                    n = 1
                else
                    cycle
                end if
            else
                n = i
            end if

            if (end > 0 .and. n > end) exit
            if (ac(n)%calculated == 0) cycle

            ! forth, momentは i = 1 の値を書き出している
            write(unit, format) ac(n)%t, ac(n)%tht, ac(n)%angle_of_attack, ac(n)%gamt, ac(n)%elevator, &
                ac(n)%mach_number, ac(n)%x, ac(n)%z, ac(n)%u_earth, ac(n)%w_earth, &
                ac(n)%u, ac(n)%w, ac(n)%q, ac(n)%cx, ac(n)%cm, ac(n)%cz, ac(n)%thrust, ac(n)%ang_thrust, &
                ac(1)%forth, ac(1)%moment, ac(n)%status
        end do

        close(unit)
    end subroutine output_path

end module mod_flight