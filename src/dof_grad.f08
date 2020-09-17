!===============================================================================
! 最適化モジュール
!===============================================================================
module mod_grad
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
    public :: workspace2d, ws2d_variables, n_var_elv, n_var_thrust
    public :: evaluate_grad, get_grad


    ! 最適化用型
    ! type, extends(TNSGA2C) :: TOpt
    !     contains
    !     procedure :: evaluate_pop => evaluate_pop_para
    !     ! generic :: evaluate => evaluate_pop_para
    !     ! procedure :: evaluate_pop_para
    ! end type TOpt


    ! 最適化用管理変数
    ! type(TOpt) :: optimizer

    ! 最適化に使用する飛行経路変数
    type(aircraft), allocatable :: workspace2d(:, :)
    real(8), allocatable :: ws2d_variables(:, :), variables_prev(:, :)


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

    integer :: n_obj_ref, n_con_ref, n_var_ref


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

    ! subroutine evaluate_pop_para(this, population)
    !     class(TOpt), intent(inout) :: this
    !     type(TPopulation), intent(inout) :: population(:)
    !     integer :: n, i!, unit

    !     n = this%current_step

    !     ! do while (access("pause", " ") == 0)
    !     !     call sleep(1)
    !     ! end do

    !     ! open(newunit=unit, file="progress.txt", status="replace")
    !     ! write(unit, "(i0)") n
    !     ! close(unit)

    !     !$omp parallel
    !     !$omp do
    !     do i = 1, size(population)
    !         call evaluate_indiv(population(i)%indiv, i)
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    !     contains

    !     subroutine evaluate_indiv(indiv, idx)
    !         class(TIndiv), intent(inout) :: indiv
    !         integer, intent(in) :: idx

    !         if (indiv%evaluated) return
    !         ! print "(2(a12,i5),a3,i5,a5,i6)",        &
    !         !       "generation:", this%current_step, &
    !         !       "evaluate:", idx,                 &
    !         !       "/", this%pop_size,               &
    !         !       "id:", indiv%id

    !         call evaluate_ga(indiv%dvariables,  &
    !                          indiv%objectives,  &
    !                          indiv%constraints, &
    !                          indiv%feasible,    &
    !                          idx)
    !         indiv%evaluated = .true.
    !     end subroutine evaluate_indiv

    ! end subroutine evaluate_pop_para


    ! #######################################################################
    ! # 型結合手続きの定義
    ! #######################################################################

    ! GAの初期化
    ! subroutine init_ga(n_var, n_obj, n_con, pop_size, xover)
    !     integer, intent(in) :: n_var, n_obj, n_con, pop_size
    !     character(*), intent(in) :: xover

    !     call optimizer%initialize(nx=n_var,        &
    !                               m=n_obj,         &
    !                               N=pop_size,      &
    !                               selection="T",   &
    !                               crossover=xover, &
    !                               mutation="PM")

    !     n_obj_ref = n_obj
    !     n_con_ref = n_con
    ! end subroutine init_ga


    ! GA実行
    ! subroutine optimize(n_gen, best)
    !     integer, intent(in) :: n_gen
    !     real(8), intent(out), allocatable :: best(:)

    !     ! call optimizer%set_problem(problem)
    !     call optimizer%prepare_calculation
    !     call optimizer%run(n_gen, hook)
    !     call optimizer%best(best)

    !     contains

    !     subroutine hook(ga)
    !         class(TSOGA), intent(inout) :: ga
    !         character(:), allocatable :: filename!, dirname
    !         integer :: n, m

    !         n = ga%current_step
    !         m = 10 ** max(int(log10(dble(max(n, 1))))-1, 0)

    !         if (n == n_gen .or. mod(n, m) == 0) then
    !             filename = "result/out_" // str(ga%current_step, 5) // ".csv"
    !             ! [saving current population]
    !             call system("if not exist result mkdir result")
    !             call ga%save_result(filename)

    !             ! [saving history]
    !             ! call system("if not exist result mkdir result")
    !             call ga%save_history("result/history.csv")

    !             ! [saving path]
    !             if (mod(n, 100) == 0) call save_flight(ga)

    !             ! dirname = "path/gen_" // str(ga%current_step)
    !             ! call system('if not exist "' // dirname // '" mkdir "' // dirname // '"')
    !             ! do i = 1, size(ga%population)
    !             ! block
    !             !     associate (acs => workspace2d(:, i),    &
    !             !                variables => indiv(i)%dvariables, &
    !             !                id => indiv(i)%id)
    !             !     class(TIndiv), allocatable :: indiv(:)
    !             !     allocate(indiv(size(ga%population)), source=ga%population%indiv)
    !                 ! do i = 1, size(indiv)
    !                 !     indiv(i) = ga%population(i)%indiv
    !                 ! end do
    !                 ! call save_flight(indiv, i, dirname)
    !             ! end block
    !             ! end do
    !         end if
    !     end subroutine hook
    ! end subroutine optimize


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

    ! subroutine evaluate_ga(variables, objectives, constraints, feasible, idx)
    !     real(8), intent(in) :: variables(:)
    !     real(8), intent(out), allocatable :: objectives(:), constraints(:)
    !     logical, intent(out) :: feasible
    !     integer, intent(in) :: idx
    !     real(8) :: ll
    !     integer :: n_last, end_status, i

    !     ! type(aircraft), allocatable :: ac(:)
    !     ! real(8), allocatable :: ll(:)
    !     ! integer :: n, i
    !     ! real(8) :: i1, i2

    !     ! variables: 設計変数の配列 (入力)
    !     ! objectives: 目的関数の配列 (出力)
    !     ! constraints: 制約条件の配列 (出力)
    !     ! feasible: 個体の有効性 (出力) [.true.=>制約条件満足, .false.=>制約違反]

    !     associate (acs => workspace2d(:, idx))
    !         acs(2:)%calculated = 0

    !         call set_elevator(acs, variables(1:n_var_elv))
    !         call set_thrust(acs, variables(n_var_elv+1:n_var_elv+n_var_thrust))

    !         ! 飛行経路を求める
    !         call flight(acs, n_last, end_status)
    !         ! if (mod(id, 1000) == 0) then
    !         !     call output_path(acs, start=1, end=-1, step=50, &
    !         !                      file="path/flight_"//str(id)//".csv")
    !         ! end if

    !         ! 評価結果を代入
    !         allocate(objectives(n_obj_ref))
    !         allocate(constraints(n_con_ref))

    !         associate (aoa0 => acs(1)%angle_of_attack * RADIANS, &
    !                    elv0 => acs(1)%elevator)
    !             ll = sum(2 * ((acs(1:n_last)%z - acs(1)%z) * cos(aoa0)       &
    !                         - (acs(1:n_last)%x - acs(1)%x) * sin(aoa0)) ** 2 &
    !                  + 0.04 * (acs(1:n_last)%elevator - elv0) ** 2)
    !         end associate

    !         do i = 1, size(objectives)
    !             select case (i)
    !             case (1)
    !                 objectives(i) = abs(acs(n_last)%gamt) + ll * dt
    !             case (2)
    !                 ! objectives(2) = sum(acs(1:n_last)%thrust) * dt
    !                 objectives(i) = abs(acs(n_last)%w_earth)
    !             case (3)
    !                 objectives(i) = abs(4770 - acs(n_last)%x)
    !             end select
    !         end do

    !         do i = 1, size(constraints)
    !             select case (i)
    !             case (1)
    !                 if (end_status == 1) then
    !                     constraints(i) = 0
    !                 else
    !                     constraints(i) = -1e10
    !                 end if
    !             case (2)
    !                 ! constraints(1) = max(acs(n_last)%angle_of_attack - 20, &
    !                 !                      0.0, &
    !                 !                      -acs(n_last)%angle_of_attack)
    !                 ! constraints(1) = (2 + acs(n_last)%w_earth) * 100
    !                 constraints(i) = 5 - abs(4770 - acs(n_last)%x)
    !             end select
    !         end do

    !         ! feasible = end_status == 1 .and. &
    !         !            acs(n_last)%x > 4769.5 .and. &
    !         !            acs(n_last)%x < 4770.5
    !         feasible = all(constraints >= 0)
    !     end associate
    ! end subroutine evaluate_ga

    subroutine get_grad(variables, objective, grads)
        real(8), intent(in) :: variables(:)
        real(8), intent(out) :: objective
        real(8), intent(out), allocatable :: grads(:)
        real(8) :: h0 = 1d-3
        real(8), allocatable :: objectives(:)
        integer :: n_last, end_status, idx_last, i
        integer :: count = 0
        ! real(8) :: obj_last = 0

        n_var_ref = size(variables)
        allocate(objectives(n_var_ref), source=0d0)
        allocate(grads(n_var_ref), source=0d0)


        ! variables = 0
        call evaluate_grad(variables, objective, n_last, end_status, 0)
        if (end_status /= 1) then
            print *, "end_status:", end_status
            stop 1
        end if

        ! print "(a, i4, es20.5, i8, l3)", "count", count, objective, n_last, all(variables==0.5d0)

        call output_path(workspace2d(:, 0), start=1, end=-1, step=50, &
                         file="path/flight_" // str(count, 4) // ".csv")


        idx_last = n_last * n_var_ref / size(workspace2d, dim=1)
        ws2d_variables = spread(variables, 2, n_var_ref)
        workspace2d(:, 1:) = spread(workspace2d(:, 0), 2, n_var_ref)

        ! print *, idx_last, n_last, n_var_ref, size(workspace2d, dim=1)

        !$omp parallel
        !$omp do
        do i = 1, idx_last
            ! if (i /= modulo(count, idx_last) + 1) cycle
            ! if (i > 1) cycle

            block
                real(8) :: f1, f2
                integer :: n_last_p, end_status_p

                ! if (variables(i) + dxs >= 0 .and. variables(i) + dxs <= 1) then
                !     dx = dxs
                ! else
                !     dx = -dxs
                ! end if
                ws2d_variables(i, i) = variables(i) + h0
                call evaluate_grad(ws2d_variables(:, i), f1, n_last_p, end_status_p, i)
                if (end_status_p /= 1) cycle

                ws2d_variables(i, i) = variables(i) - h0
                call evaluate_grad(ws2d_variables(:, i), f2, n_last_p, end_status_p, i)
                if (end_status_p /= 1) cycle

                grads(i) = (f1 - f2) / (2 * h0)

                ! block
                !     integer :: j

                !     open(100, file="test/test1_"//str(count+1)//".txt")
                !     write(100, *) objectives(i)
                !     do j = 1, size(grads)
                !         write(100, *) j, ws2d_variables(j, i)
                !     end do
                !     close(100)
                ! end block

                ! if (end_status_p == 1) then
                !     grads(i) = (objectives(i) - objective) / dx
                ! end if
                ! print "(a,i5,3es18.5,l5)", "result", i, objective, objectives(i), grads(i), obj_last == objective
                ! if (count > 0 .and. abs(obj_last - objective) > 1) then
                !     print *, obj_last
                !     print *, objective
                !     stop 2
                ! end if
                ! obj_last = objectives(i)
            end block
        end do
        !$omp end do
        !$omp end parallel

        ! block
        !     integer :: j

        !     open(100, file="test/test0_"//str(count)//".txt")
        !     write(100, *) objective
        !     do j = 1, n_var_ref
        !         write(100, *) j, variables(j), grads(j)
        !     end do
        !     close(100)
        ! end block

        count = count + 1

    end subroutine get_grad


    subroutine evaluate_grad(variables, objective, n_last, end_status, idx)
        real(8), intent(in) :: variables(:)
        real(8), intent(out) :: objective
        integer, intent(out) :: n_last, end_status
        integer, intent(in) :: idx
        real(8) :: ll
        integer :: start
        integer :: count = 1

        start = max((idx - 1) * size(workspace2d, dim=1) / n_var_ref, 1)
        ! print *, "start", start, idx, size(workspace2d, dim=1), n_var_ref
        ! if (idx > 1) return

        associate (acs => workspace2d(:, idx))
            acs(start+1:)%calculated = 0
            call set_elevator(acs, variables(1:n_var_elv))
            call set_thrust(acs, variables(n_var_elv+1:n_var_elv+n_var_thrust))
            call flight(acs, n_last, end_status, start=start)

            ! print "(a,4i5)", "flight", count, n_last-start, start, n_last

            associate (aoa0 => acs(1)%angle_of_attack * RADIANS, &
                       elv0 => acs(1)%elevator)
                ll = sum(2 * ((acs(1:n_last)%z - acs(1)%z) * cos(aoa0)       &
                            - (acs(1:n_last)%x - acs(1)%x) * sin(aoa0)) ** 2 &
                     + 0.04 * (acs(1:n_last)%elevator - elv0) ** 2)
            end associate

            objective = abs(acs(n_last)%gamt) + ll * dt
            ! call output_path(acs, start=1, end=-1, step=50, &
            !                  file="path/flight_" // str(i, 4) // ".csv")
        end associate
    end subroutine evaluate_grad

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
        ! acs(1)%elevator = elv_start + (array(1) * 20 - 10)

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

end module mod_grad


!===============================================================================
! メインモジュール (実行ファイル作成用)
!===============================================================================

program dof_main
    use util, only: elapsed_time, str
    use mod_kriging
    use global
    use mod_aircraft
    use mod_flight
    use mod_grad
    implicit none


    ! Krigingモデルの初期化
    call init_krig("kriging1.txt")
    call optimize_test ! [B]


    contains


    ! 運動計算メイン (最適化)
    subroutine optimize_test
        ! type(aircraft), allocatable :: ac(:)
        real(8) :: flight_time, alpha
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
        alpha = 0.05d0

        n_var_elv = int(flight_time * n_control)    ! 設計変数の数 (の定義点数)
        ! n_var_thrust = int(flight_time * n_control) ! 設計変数の数 (推力の定義点数)
        n_var_thrust = 0 ! 設計変数の数 (推力の定義点数)
        n_var = n_var_elv + n_var_thrust ! 設計変数の数 (舵角と推力の定義点数)
        ! n_obj = 3                        ! 目的関数の数

        n_step = int(flight_time / dt) ! 計算ステップ数
        ! start = 1    ! 最適化開始ステップ数
        ! end = n_step ! 最適化終了ステップ数


        ! 初期化
        allocate(workspace2d(n_step, 0:n_var))
        allocate(ws2d_variables(n_var, n_var))
        call init_aircraft(workspace2d(1, 0), "1")
        workspace2d(1, 1:) = workspace2d(1, 0)

        ! call init_ga(n_var, n_obj, n_con, pop_size, xover)

        ! 計算時間計測開始
        call system_clock(start_time)

        ! 最適化
        ! call optimize(n_gen, best)
        block
            real(8) :: variables(n_var), x_diff(n_var)
            real(8), allocatable :: grads(:)
            real(8) :: objective, result(n_gen)
            integer :: i, j

            variables = 0.5d0

            do i = 1, n_gen
                call get_grad(variables, objective, grads)

                result(i) = objective
                x_diff = alpha * grads / result(1)

                print "(a, i5, 3es18.5)", "step", i, objective, &
                                          objective - result(max(i-1,1)), &
                                          objective - result(1)

                ! if (i > 1 .and. objective > result(max(i-1,1))) stop 1

                open(100, file="test.txt")
                write(100, *) objective
                do j = 1, size(grads)
                    write(100, *) j, variables(j), grads(j), x_diff(j)
                end do
                close(100)
                call system("copy test.txt test\\test_" // str(i) // ".txt > nul")

                variables = min(max(variables - x_diff, 0d0), 1d0)
                ! where (x_diff > 0)
                !     variables = variables - 1d-2
                ! else where (x_diff < 0)
                !     variables = variables + 1d-2
                ! end where

                open(200, file="result.csv")
                write(200, "(a)") "step,J"
                do j = 1, i
                    write(200, "(i0,a,es20.5)") j, ",", result(j)
                end do
                close(200)
            end do
            ! print *, objective, n_last

        end block

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
