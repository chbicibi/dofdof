!===============================================================================
! 最適化モジュール
!===============================================================================
module dof_flight
    use global
    use mod_aircraft
    use dof_kriging
    implicit none

    private
    public :: read_input, init_aircraft, flight, output_path

    real(8) :: x_start, y_start, z_start
    real(8) :: u_start, v_start, w_start
    real(8) :: p_start, q_start, r_start
    real(8) :: phi_start, tht_start, psi_start

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


    ! 計算条件の確認
    logical function check_stop(ac, exit_status) result(is_stop)
        type(aircraft), intent(inout) :: ac
        integer, intent(out) :: exit_status

        if (ac%z <= 0) then
            is_stop = .true.
            exit_status = 1

        else if (ac%angle_of_attack < 0 .or. ac%angle_of_attack > 20) then
            ! print "(a,f0.1,a)", 'Error: "Angle of Attack" is out of range! (', ac%t, "s)"
            is_stop = .true.
            exit_status = 2

        else if (ac%mach_number < 0.1d0 .or. ac%mach_number > 0.5d0) then
            ! print "(a,f0.1,a)", 'Error: "Mach Number" is out of range! (', ac%t, "s)"
            is_stop = .true.
            exit_status = 3

        else if (ac%elevator < -50 .or. ac%elevator > 30) then
            ! print "(a,f0.1,a)", 'Error: "Elevator" is out of range! (', ac%t, "s)"
            is_stop = .true.
            exit_status = 4

        else
            is_stop = .false.
            exit_status = 0
        end if
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
            write(unit, "(a)") "time,theta,alpha,gamt,elv,ma,x,z,ue,we,uu,ww,q,cx,cm,cz,thrust,ang_thrust,forth,moment,calculated"
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

end module dof_flight


!===============================================================================
! 最適化モジュール
!===============================================================================
module dof_ga
    use util, only: str
    use global
    use individual, only: TIndiv
    use dof_kriging, only: init_krig, estimate
    use soga, only: TSOGA, TPopulation
    use nsga2c, only: TNSGA2C
    use mod_aircraft
    use dof_flight
    implicit none

    private
    public :: optimizer, workspace2d, init_ga, optimize, save_result
    public :: n_var_elv, n_var_thrust


    ! 最適化用型
    type, extends(TNSGA2C) :: TOpt
        contains
        procedure :: evaluate_pop => evaluate_pop_para
        ! generic :: evaluate => evaluate_pop_para
        ! procedure :: evaluate_pop_para
    end type TOpt


    ! 最適化用管理変数
    type(TOpt) :: optimizer

    ! 最適化に使用する飛行経路変数
    type(aircraft), allocatable :: workspace2d(:, :)


    ! 設計変数を定義
    real(8) :: elv_lower = -50.0     ! 舵角入力下限値 (度)
    real(8) :: elv_upper = 30.0      ! 舵角入力上限値 (度)
    real(8) :: elv_bound = 2.0       ! 舵角入力変化量最大値 (度/秒)
    real(8) :: elv_start  = -25.9136143        ! 舵角入力初期値 (度/秒)
    integer :: n_var_elv = 400

    real(8) :: thrust_lower =  -140000 ! 推力入力下限値 (N)
    real(8) :: thrust_upper =  140000 ! 推力入力上限値 (N)
    real(8) :: thrust_bound =  1000   ! 推力入力変化量最大値 (N/秒)
    real(8) :: thrust_start  =  17719.3379    ! 推力入力初期値 (N)
    integer :: n_var_thrust = 400


    ! ! 吉田追加
    ! real(8) :: ang_thrust_lower = -5.0 ! 推力方向下限値 (deg.)
    ! real(8) :: ang_thrust_upper =  5.0 ! 推力方向上限値 (deg.)
    ! real(8) :: ang_thrust_bound =  1   ! 推力方向変化量最大値 (deg./秒)
    ! real(8) :: ang_thrust_start  =  0    ! 推力方向初期値 (deg.)

    ! integer :: n_var_ang_thrust, n_var_thrust, n_last

    ! ! 空力データベース範囲
    ! real(8) :: rel_z_low = 0.05
    ! real(8) :: rel_z_upp = 1.0
    ! real(8) :: rel_x_low = 0.00
    ! real(8) :: rel_x_upp = 0.5

    contains

    ! #######################################################################
    ! # 型結合手続きの定義
    ! #######################################################################

    subroutine evaluate_pop_para(this, population)
        class(TOpt), intent(inout) :: this
        type(TPopulation), intent(inout) :: population(:)
        integer :: i

        !$omp parallel
        !$omp do
        do i = 1, size(population)
            call evaluate_indiv(population(i)%indiv, i)
        end do
        !$omp end do
        !$omp end parallel

        contains

        subroutine evaluate_indiv(indiv, idx)
            class(TIndiv), intent(inout) :: indiv
            integer, intent(in) :: idx

            if (indiv%evaluated) return
            ! print "(2(a12,i5),a3,i5,a5,i6)",        &
            !       "generation:", this%current_step, &
            !       "evaluate:", idx,                 &
            !       "/", this%pop_size,               &
            !       "id:", indiv%id

            call evaluate_ga(indiv%dvariables,  &
                             indiv%objectives,  &
                             indiv%constraints, &
                             indiv%feasible,    &
                             idx, indiv%id)
            indiv%evaluated = .true.
        end subroutine evaluate_indiv

    end subroutine evaluate_pop_para


    ! #######################################################################
    ! # 型結合手続きの定義
    ! #######################################################################

    ! GAの初期化
    subroutine init_ga(n_var, n_obj, pop_size)
        integer, intent(in) :: n_var, n_obj, pop_size

        call optimizer%initialize(nx=n_var,        &
                                  m=n_obj,         &
                                  N=pop_size,      &
                                  selection="T",   &
                                  crossover="BLX", &
                                  mutation="PM")

        ! call optimizer%set_fitness_type("VALUE")
        optimizer%elite_preservation = .false.
        optimizer%sharing = .true.
        optimizer%history_preservation = .true.
    end subroutine init_ga


    ! GA実行
    subroutine optimize(n_gen, best)
        integer, intent(in) :: n_gen
        real(8), intent(out), allocatable :: best(:)

        ! call optimizer%set_problem(problem)
        call optimizer%prepare_calculation
        call optimizer%run(n_gen, hook)
        call optimizer%best(best)

        contains

        subroutine hook(ga)
            class(TSOGA), intent(inout) :: ga
            character(:), allocatable :: filename
            integer :: n

            n = ga%current_step

            if (n == n_gen .or. mod(n, 100) == 0) then
                filename = "result/out_" // str(ga%current_step, 5) // ".csv"
                call ga%save_result(filename)
            end if

            if (n == n_gen .or. mod(n, 100) == 0) then
                call ga%save_history("result/history.csv")
            end if
        end subroutine hook
    end subroutine optimize


    ! GA結果ファイル出力
    subroutine save_result

        ! 制約出力なし(仮)
        call optimizer%save_result("out_so_X_1.csv", elite="only")
        call optimizer%save_result("out_so_X_2.csv", elite="all")
        call optimizer%save_history("out_so_X_3.csv", elite="only")
        call optimizer%save_history("out_so_X_4.csv", elite="all")
        call optimizer%save_history("out_so_X_5.csv", elite="all")

        ! call optimizer%save_result("out_so_X_1.csv", elite="only", feasible="all")
        ! call optimizer%save_result("out_so_X_2.csv", elite="all", feasible="all")
        ! call optimizer%save_history("out_so_X_3.csv", elite="only", feasible="all")
        ! call optimizer%save_history("out_so_X_4.csv", elite="all", feasible="all")
        ! call optimizer%save_history("out_so_X_5.csv", elite="all", feasible="only")

        ! call optimizer%save_result("test1.csv")
        ! call optimizer%save_history("test2.csv")
    end subroutine save_result


    ! #######################################################################
    ! # 目的関数を計算する
    ! # このサブルーチンを編集する
    ! #######################################################################

    subroutine evaluate_ga(variables, objectives, constraints, feasible, idx, id)
        real(8), intent(in) :: variables(:)
        real(8), intent(out), allocatable :: objectives(:), constraints(:)
        logical, intent(out) :: feasible
        integer, intent(in) :: idx, id
        real(8) :: ll
        integer :: n_last, end_status, i

        ! type(aircraft), allocatable :: ac(:)
        ! real(8), allocatable :: ll(:)
        ! integer :: n, i
        ! real(8) :: i1, i2

        ! variables: 設計変数の配列 (入力)
        ! objectives: 目的関数の配列 (出力)
        ! constraints: 制約条件の配列 (出力)
        ! feasible: 個体の有効性 (出力) [.true.=>制約条件満足, .false.=>制約違反]

        ! 機体の状態を初期化
        ! do i = 1, size(ac)
        !     call init_aircraft(ac(i))
        ! end do

        ! 設計変数をエレベータと推力に変換し，飛行経路を求める
        ! call set_elv(variables(:n_var_elv), ac_opt)
        ! call set_thrust(variables(n_var_elv+1:), ac_opt)

        associate (acs => workspace2d(:, idx))
            acs(2:)%calculated = 0

            call set_elevator(acs, variables(1:n_var_elv))
            call set_thrust(acs, variables(n_var_elv+1:n_var_elv+n_var_thrust))

            ! 飛行経路を求める
            call flight(acs, n_last, end_status)
            if (mod(id, 1000) == 0) then
                call output_path(acs, start=1, end=-1, step=50, &
                                 file="path/flight_"//str(id)//".csv")
            end if

            ! 評価結果を代入
            allocate(objectives(2))
            allocate(constraints(2))

            ll = 0
            do i = 1, n_last
                ll = ll + 2 * ((acs(i)%z - acs(1)%z) * cos(acs(1)%angle_of_attack)       &
                             - (acs(i)%x - acs(1)%x) * sin(acs(1)%angle_of_attack)) ** 2 &
                        + 0.04 * (acs(i)%elevator - acs(1)%elevator) ** 2
            end do
            objectives(1) = abs(acs(n_last)%gamt) + ll * dt
            objectives(2) = sum(acs(1:n_last)%thrust) * dt

            constraints(1) = max(acs(n_last)%angle_of_attack - 20, &
                                 0.0, &
                                 -acs(n_last)%angle_of_attack)
            constraints(2) = abs(acs(n_last)%x - 4770)

            feasible = end_status == 1 .and. &
                       acs(n_last)%x > 4769.5 .and. &
                       acs(n_last)%x < 4770.5
        end associate
    end subroutine evaluate_ga


    ! #######################################################################
    ! # 飛行制御パラメータ設定
    ! #######################################################################

    subroutine set_elevator(acs, array)
        type(aircraft), intent(inout) :: acs(:)
        real(8), intent(in) :: array(:)
        real(8) :: r, d_elv, elv
        integer :: in_size, out_size, idx, i

        in_size = size(array)
        out_size = size(acs)
        r = dble(in_size) / out_size

        acs(1)%elevator = elv_start
        do i = 2, out_size
            idx = int(r * (i - 1)) + 1
            d_elv = (2 * array(idx) - 1) * elv_bound
            elv = acs(i - 1)%elevator + d_elv * dt
            acs(i)%elevator = min(max(elv, elv_lower), elv_upper)
        end do
    end subroutine set_elevator


    subroutine set_thrust(acs, array)
        type(aircraft), intent(inout) :: acs(:)
        real(8), intent(in) :: array(:)
        real(8) :: r, d_thrust, thrust
        integer :: in_size, out_size, idx, i

        in_size = size(array)
        out_size = size(acs)
        r = dble(in_size) / out_size

        acs(1)%thrust = thrust_start
        do i = 2, out_size
            idx = int(r * (i - 1)) + 1
            d_thrust = (2 * array(idx) - 1) * thrust_bound
            thrust = acs(i - 1)%thrust + d_thrust * dt
            acs(i)%thrust = min(max(thrust, thrust_lower), thrust_upper)
        end do
    end subroutine set_thrust


    ! subroutine set_ang_thrust(variables, acs)
    !     real(8), intent(in) :: variables(:)
    !     type(aircraft), intent(inout) :: acs(:)

    !     ! #############################################
    !     ! # variables配列を推力方向配列に変換する処理を書く
    !     ! # 必要に応じて線形補間, スプライン補間等を使用する
    !     ! # variablesは[0, 1]の範囲を持つ実数配列
    !     ! #############################################

    !     real(8) :: r       ! 入出力比
    !     real(8) :: d_input ! 推力方向変化量 (度/秒)
    !     integer :: in_size, out_size, i

    !     in_size = size(variables)    ! 入力サイズ (設計変数の数)
    !     out_size = size(acs)         ! 出力サイズ (制御入力数)
    !     r = real(in_size) / out_size

    !     do i = 0, out_size - 1
    !         ! [0, 1] => [-1, 1] => [-ang_thrust_bound, ang_thrust_bound] と変換
    !         ! d_input = (2 * variables(int(r * i) + 1) - 1) * ang_thrust_bound

    !         ! if (i == 0) cycle

    !         ! 前ステップの値に加算
    !         ! acs(i+1)%ang_thrust = min(max(acs(i)%ang_thrust + (d_input * dt), ang_thrust_lower), ang_thrust_upper)
    !     end do
    ! end subroutine set_ang_thrust


    ! subroutine set_thrust(variables, acs)
    !     real(8), intent(in) :: variables(:)
    !     type(aircraft), intent(inout) :: acs(:)

    !     ! #############################################
    !     ! # variables配列を推力配列に変換する処理を書く
    !     ! # 必要に応じて線形補間, スプライン補間等を使用する
    !     ! # variablesは[0, 1]の範囲を持つ実数配列
    !     ! #############################################

    !     real(8) :: r       ! 入出力比
    !     real(8) :: d_input ! 推力変化量 (N/秒)
    !     integer :: in_size, out_size, i

    !     in_size = size(variables)    ! 入力サイズ (設計変数の数)
    !     out_size = size(acs)         ! 出力サイズ (制御入力数)
    !     r = real(in_size) / out_size

    !     do i = 0, out_size - 1
    !         ! [0, 1] => [-1, 1] => [-thrust_bound, thrust_bound] と変換
    !         ! d_input = (2 * variables(int(r * i) + 1) - 1) * thrust_bound

    !         ! if (i == 0) cycle

    !         ! 前ステップの値に加算
    !         ! acs(i+1)%thrust = min(max(acs(i)%thrust + (d_input * dt), thrust_lower), thrust_upper)
    !     end do
    ! end subroutine set_thrust

end module dof_ga


!===============================================================================
! メインモジュール (実行ファイル作成用)
!===============================================================================

program dof_main
    use util, only: elapsed_time
    use global
    use mod_aircraft
    use dof_kriging
    use dof_ga
    use dof_flight
    implicit none

    ! integer :: is_init_krig
    ! real(8) :: rel_z_low, rel_z_upp, rel_x_low, rel_x_upp
    ! real(8) :: t1, t2



    ! Krigingモデルの初期化
    call init_krig("kriging1.txt")
    call optimize_test ! [B]


    contains


    ! 運動計算メイン (最適化)
    subroutine optimize_test
        ! type(aircraft), allocatable :: ac(:)
        real(8) :: flight_time
        real(8), allocatable :: best(:)
        integer :: start_time, n_control
        integer :: n_var, n_obj, pop_size, n_gen, n_step, start, end

        dt = 0.005         ! 時間刻み幅を設定
        flight_time = 80   ! 飛行時間
        n_control = 5      ! 操舵回数 (回/秒)
        pop_size = 30      ! GA集団サイズ
        n_gen = 10000      ! GA世代数

        n_step = int(flight_time / dt) ! 計算ステップ数
        ! n_var_elv = flight_time * n_control    ! 設計変数の数 (の定義点数)
        ! n_var_thrust = flight_time * n_control ! 設計変数の数 (推力の定義点数)

        n_var = n_var_elv + n_var_thrust ! 設計変数の数 (舵角と推力の定義点数)
        n_obj = 2                        ! 目的関数の数

        start = 1    ! 最適化開始ステップ数
        end = n_step ! 最適化終了ステップ数


        ! 初期化
        allocate(workspace2d(n_step, pop_size))
        call init_aircraft(workspace2d(1, 1), "1")
        ! block
        !     integer :: i
        !     workspace2d(:, 1)%t = [(dt*(i-1),i = 1, n_step)]
        ! end block
        workspace2d(1, 2:) = workspace2d(1, 1)

        call init_ga(n_var, n_obj, pop_size)

        ! 計算時間計測開始
        call system_clock(start_time)

        ! 最適化
        call optimize(n_gen, best)

        ! 計算時間を表示
        call elapsed_time(start_time)

        ! call save_result

        ! ! 最適結果取得
        ! call set_ang_thrust(best(5:n_var_ang_thrust+5), ac2_opt)
        ! call set_thrust(best(n_var_ang_thrust+6:), ac2_opt)

        ! ! 最適結果で飛行経路を計算
        ! call flight(ac_opt, ac1_opt, ac2_opt, t1, t2)

        ! ! 結果を出力
        ! call output_path(ac_opt, start=0, end=-1, step=50, file="path.csv")
        ! call output_path(ac1_opt, start=0, end=-1, step=50, file="path1.csv")
        ! call output_path(ac2_opt, start=0, end=-1, step=50, file="path2.csv")
    end subroutine optimize_test


end program dof_main
