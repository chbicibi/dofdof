
!===============================================================================
! メインモジュール (実行ファイル作成用)
!===============================================================================

program dof_main
    use util, only: elapsed_time
    use global
    use mod_aircraft
    use dof_kriging
    implicit none

    real(8) :: x_start, y_start, z_start
    real(8) :: u_start, v_start, w_start
    real(8) :: p_start, q_start, r_start
    real(8) :: phi_start, tht_start, psi_start
    integer :: is_init_krig
    character(:), allocatable :: tes


    ! Krigingモデルの初期化
    call init_krig("kriging1.txt")
    call flight_test ! [A]


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
        ac%elevator = -25.9136143
        ac%thrust = 17719.3379
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
            write(unit, "(a)") "time,theta,alpha,gamt,elv,ma,x,z,ue,we,uu,ww,q,cx,cm,cz,thrust,ang_thrust,forth,moment"
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
                ac(1)%forth, ac(1)%moment, ac(n)%calculated
        end do

        close(unit)
    end subroutine output_path


    ! ==========================================================================
    ! テスト用
    ! ==========================================================================

    ! [A] 運動計算メイン (飛行のみ)
    subroutine flight_test
        use util, only: str
        type(aircraft), allocatable :: ac(:, :)
        real(8) :: angle_of_attack, elevator, mach_number, rel_x, rel_z, del_tht
        integer :: total_steps, i, phase, j

        dt = 0.005 ! 時間刻み幅を設定
        total_steps = 20000  ! 計算ステップ数

        ! 配列変数の割付
        if (allocated(ac)) then
            if (size(ac) /= total_steps) then
                deallocate(ac)
                allocate(ac(total_steps, 10))
            end if
        else
            allocate(ac(total_steps, 10))
        end if

        ! 機体の状態を初期化
        ! 制御入力が必要な場合はここを書き換える (set_elv, set_thrustを使用)
        call init_aircraft(ac(1, 1), "1")
        ac(1, 2:) = ac(1, 1)
        ! ac1(1)%elevator = 0
        ! ac1(1)%thrust = 0
        ! ac1(1)%ang_thrust = 0
        ! ac2(1)%elevator = 0
        ! ac2(1)%thrust = 1960000
        ! ac2(1)%ang_thrust = ac2(1)%tht
        ! ac1(2:) = ac1(1)
        ! ac2(2:) = ac2(1)

        ! stop

        ! 運動メインループ
        !$omp parallel
        !$omp do
        do j = 1, 10
            do i = 1, total_steps - 1
                ac(i, j)%t = dt * (i - 1)
                if (i > 1) then
                    ac(i, j)%elevator = ac(i-1, j)%elevator
                    ac(i, j)%thrust = ac(i-1, j)%thrust
                    ac(i, j)%ang_thrust = ac(i-1, j)%ang_thrust
                end if

                ac(i, j)%cx = estimate(1, [ac(i, j)%angle_of_attack, ac(i, j)%elevator, ac(i, j)%mach_number])
                ac(i, j)%cm = estimate(2, [ac(i, j)%angle_of_attack, ac(i, j)%elevator, ac(i, j)%mach_number])
                ac(i, j)%cz = estimate(3, [ac(i, j)%angle_of_attack, ac(i, j)%elevator, ac(i, j)%mach_number])

                call calc_motion(ac(i, j), ac(i+1, j))

                ! if (i == 1) then
                !     call ac(1)%print_state
                ! end if

                if (mod(i, 1000) == 0) print *, i, "step calculating...", j
            end do

            ! 結果を出力
            call output_path(ac(:, j), start=1, end=-1, step=50, file="path"//str(j)//".csv")
        end do
        !$omp end do
        !$omp end parallel
    end subroutine flight_test

end program dof_main