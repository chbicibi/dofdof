!===============================================================================
! 最適化モジュール
!===============================================================================
module mod_ga
    use util, only: str
    use mod_kriging, only: init_krig, estimate
    use individual, only: TIndiv
    use soga, only: TSOGA, TPopulation
    use nsga2c, only: TNSGA2C
    use global
    use mod_aircraft
    use mod_flight
    implicit none

    private
    public :: optimizer, workspace2d, init_ga, optimize
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
    real(8) :: elv_upper = 10.0      ! 舵角入力上限値 (度)
    real(8) :: elv_bound = 2.0       ! 舵角入力変化量最大値 (度/秒)
    real(8) :: elv_start  = -25.9136143        ! 舵角入力初期値 (度/秒)
    integer :: n_var_elv = 400

    real(8) :: thrust_lower =  -140000 ! 推力入力下限値 (N)
    real(8) :: thrust_upper =  140000 ! 推力入力上限値 (N)
    real(8) :: thrust_bound =  1000   ! 推力入力変化量最大値 (N/秒)
    ! real(8) :: thrust_start  =  17719.3379    ! 推力入力初期値 (N)
    real(8) :: thrust_start  =  0    ! 推力入力初期値 (N)
    integer :: n_var_thrust = 400

    integer :: n_obj_ref, n_con_ref


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
        integer :: n, i!, unit

        n = this%current_step

        ! do while (access("pause", " ") == 0)
        !     call sleep(1)
        ! end do

        ! open(newunit=unit, file="progress.txt", status="replace")
        ! write(unit, "(i0)") n
        ! close(unit)

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
                             idx)
            indiv%evaluated = .true.
        end subroutine evaluate_indiv

    end subroutine evaluate_pop_para


    ! #######################################################################
    ! # 型結合手続きの定義
    ! #######################################################################

    ! GAの初期化
    subroutine init_ga(n_var, n_obj, n_con, pop_size, xover)
        integer, intent(in) :: n_var, n_obj, n_con, pop_size
        character(*), intent(in) :: xover

        call optimizer%initialize(nx=n_var,        &
                                  m=n_obj,         &
                                  N=pop_size,      &
                                  selection="T",   &
                                  crossover=xover, &
                                  mutation="PM")

        n_obj_ref = n_obj
        n_con_ref = n_con
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
            character(:), allocatable :: filename!, dirname
            integer :: n, m

            n = ga%current_step
            m = 10 ** max(int(log10(dble(max(n, 1))))-1, 0)

            if (n == n_gen .or. mod(n, m) == 0) then
                filename = "result/out_" // str(ga%current_step, 5) // ".csv"
                ! [saving current population]
                call system("if not exist result mkdir result")
                call ga%save_result(filename)

                ! [saving history]
                ! call system("if not exist result mkdir result")
                call ga%save_history("result/history.csv")

                ! [saving path]
                if (mod(n, 100) == 0) call save_flight(ga)

                ! dirname = "path/gen_" // str(ga%current_step)
                ! call system('if not exist "' // dirname // '" mkdir "' // dirname // '"')
                ! do i = 1, size(ga%population)
                ! block
                !     associate (acs => workspace2d(:, i),    &
                !                variables => indiv(i)%dvariables, &
                !                id => indiv(i)%id)
                !     class(TIndiv), allocatable :: indiv(:)
                !     allocate(indiv(size(ga%population)), source=ga%population%indiv)
                    ! do i = 1, size(indiv)
                    !     indiv(i) = ga%population(i)%indiv
                    ! end do
                    ! call save_flight(indiv, i, dirname)
                ! end block
                ! end do
            end if
        end subroutine hook
    end subroutine optimize


    subroutine save_flight(ga)
        class(TSOGA), intent(in) :: ga
        character(:), allocatable :: dirname
        integer :: i
        ! [saving path]
        dirname = "path/gen_" // str(ga%current_step)
        call system('if not exist "' // dirname // '" mkdir "' // dirname // '"')

        !$omp parallel
        !$omp do
        do i = 1, size(ga%population)
            call save_flight_inner(ga%population(i)%indiv, i)
        end do
        !$omp end do
        !$omp end parallel

        do i = 1, size(ga%population)
            associate (acs => workspace2d(:, i), &
                       id => ga%population(i)%indiv%id)
                call output_path(acs, start=1, end=-1, step=50, &
                                 file=dirname // "/flight_" // str(i, 4) // "_id=" // str(id) // ".csv")
            end associate
        end do

        contains

        subroutine save_flight_inner(indiv, idx)
            class(TIndiv), intent(in) :: indiv
            integer, intent(in) :: idx
            ! character(*), intent(in) :: dirname
            ! character(:), allocatable :: filename
            integer :: n_last, end_status

            associate (acs => workspace2d(:, idx),    &
                       variables => indiv%dvariables)!, &
                       ! id => indiv%id)
                ! filename = dirname // "/flight_" // str(idx, 4) // "_id=" // str(id) // ".csv"

                acs(2:)%calculated = 0
                call set_elevator(acs, variables(1:n_var_elv))
                call set_thrust(acs, variables(n_var_elv+1:n_var_elv+n_var_thrust))
                call flight(acs, n_last, end_status)

                ! call output_path(acs, start=1, end=-1, step=50, file=filename)
            end associate
        end subroutine save_flight_inner

    end subroutine save_flight


    ! GA結果ファイル出力
    ! subroutine save_result

    !     ! 制約出力なし(仮)
    !     call optimizer%save_result("out_so_X_1.csv", elite="only")
    !     call optimizer%save_result("out_so_X_2.csv", elite="all")
    !     call optimizer%save_history("out_so_X_3.csv", elite="only")
    !     call optimizer%save_history("out_so_X_4.csv", elite="all")
    !     call optimizer%save_history("out_so_X_5.csv", elite="all")

    !     ! call optimizer%save_result("out_so_X_1.csv", elite="only", feasible="all")
    !     ! call optimizer%save_result("out_so_X_2.csv", elite="all", feasible="all")
    !     ! call optimizer%save_history("out_so_X_3.csv", elite="only", feasible="all")
    !     ! call optimizer%save_history("out_so_X_4.csv", elite="all", feasible="all")
    !     ! call optimizer%save_history("out_so_X_5.csv", elite="all", feasible="only")

    !     ! call optimizer%save_result("test1.csv")
    !     ! call optimizer%save_history("test2.csv")
    ! end subroutine save_result


    ! #######################################################################
    ! # 目的関数を計算する
    ! # このサブルーチンを編集する
    ! #######################################################################

    subroutine evaluate_ga(variables, objectives, constraints, feasible, idx)
        real(8), intent(in) :: variables(:)
        real(8), intent(out), allocatable :: objectives(:), constraints(:)
        logical, intent(out) :: feasible
        integer, intent(in) :: idx
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

        associate (acs => workspace2d(:, idx))
            acs(2:)%calculated = 0

            call set_elevator(acs, variables(1:n_var_elv))
            call set_thrust(acs, variables(n_var_elv+1:n_var_elv+n_var_thrust))

            ! 飛行経路を求める
            call flight(acs, n_last, end_status)
            ! if (mod(id, 1000) == 0) then
            !     call output_path(acs, start=1, end=-1, step=50, &
            !                      file="path/flight_"//str(id)//".csv")
            ! end if

            ! 評価結果を代入
            allocate(objectives(n_obj_ref))
            allocate(constraints(n_con_ref))

            associate (aoa0 => acs(1)%angle_of_attack * RADIANS, &
                       elv0 => acs(1)%elevator)
                ll = sum(2 * ((acs(1:n_last)%z - acs(1)%z) * cos(aoa0)       &
                            - (acs(1:n_last)%x - acs(1)%x) * sin(aoa0)) ** 2 &
                     + 0.04 * (acs(1:n_last)%elevator - elv0) ** 2)
            end associate

            do i = 1, size(objectives)
                select case (i)
                case (1)
                    objectives(i) = abs(acs(n_last)%gamt) + ll * dt
                case (2)
                    ! objectives(2) = sum(acs(1:n_last)%thrust) * dt
                    objectives(i) = abs(acs(n_last)%w_earth)
                case (3)
                    objectives(i) = abs(4770 - acs(n_last)%x)
                end select
            end do

            do i = 1, size(constraints)
                select case (i)
                case (1)
                    if (end_status == 1) then
                        constraints(i) = 0
                    else
                        constraints(i) = -1e10
                    end if
                case (2)
                    ! constraints(1) = max(acs(n_last)%angle_of_attack - 20, &
                    !                      0.0, &
                    !                      -acs(n_last)%angle_of_attack)
                    ! constraints(1) = (2 + acs(n_last)%w_earth) * 100
                    constraints(i) = 5 - abs(4770 - acs(n_last)%x)
                end select
            end do

            ! feasible = end_status == 1 .and. &
            !            acs(n_last)%x > 4769.5 .and. &
            !            acs(n_last)%x < 4770.5
            feasible = all(constraints >= 0)
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

        in_size = size(array) - 1
        out_size = size(acs)
        r = dble(in_size) / out_size

        ! acs(1)%elevator = elv_start
        acs(1)%elevator = elv_start + (array(1) * 20 - 10)

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

        acs%thrust = 0
        return

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

end module mod_ga


!===============================================================================
! メインモジュール (実行ファイル作成用)
!===============================================================================

program dof_main
    use util, only: elapsed_time
    use mod_kriging
    use global
    use mod_aircraft
    use mod_flight
    use mod_ga
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
        integer :: n_var, n_obj, n_con, pop_size, n_gen, n_step!, start, end
        integer :: unit
        character(3) :: xover


        open(newunit=unit, file="params.txt", status="old")
        read(unit, *) n_obj
        read(unit, *) n_con
        read(unit, *) pop_size
        read(unit, *) n_gen
        read(unit, *) xover
        close(unit)

        dt = 0.005         ! 時間刻み幅を設定
        flight_time = 80   ! 飛行時間
        n_control = 5      ! 操舵回数 (回/秒)
        ! pop_size = 30      ! GA集団サイズ
        ! n_gen = 2000      ! GA世代数

        n_var_elv = int(flight_time * n_control) + 1    ! 設計変数の数 (の定義点数)
        ! n_var_thrust = int(flight_time * n_control) ! 設計変数の数 (推力の定義点数)
        n_var_thrust = 0 ! 設計変数の数 (推力の定義点数)
        n_var = n_var_elv + n_var_thrust ! 設計変数の数 (舵角と推力の定義点数)
        ! n_obj = 3                        ! 目的関数の数

        n_step = int(flight_time / dt) ! 計算ステップ数
        ! start = 1    ! 最適化開始ステップ数
        ! end = n_step ! 最適化終了ステップ数


        ! 初期化
        allocate(workspace2d(n_step, pop_size))
        call init_aircraft(workspace2d(1, 1), "1")
        ! block
        !     integer :: i
        !     workspace2d(:, 1)%t = [(dt*(i-1),i = 1, n_step)]
        ! end block
        workspace2d(1, 2:) = workspace2d(1, 1)

        call init_ga(n_var, n_obj, n_con, pop_size, xover)

        ! 計算時間計測開始
        call system_clock(start_time)

        ! 最適化
        call optimize(n_gen, best)

        ! 計算時間を表示
        call elapsed_time(start_time)

        ! open(newunit=unit, file="progress.txt", status="replace")
        ! write(unit, "(i0)") -1
        ! close(unit)

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
