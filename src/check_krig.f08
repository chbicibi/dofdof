! 実行ファイル作成方法:
! 方法1: バッチファイル使用
! setup.bat krig

! 方法2: 直接コマンドを実行
! gfortran -c -O3 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -Wno-unused-variable -static -ffree-line-length-0 -std=gnu -o dof_kriging.o dof_kriging.f08 -Igalib\include
! gfortran -c -O3 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -Wno-unused-variable -static -ffree-line-length-0 -std=gnu -o check_krig.o check_krig.f08 -Igalib\include
! gfortran -O3 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -Wno-unused-variable -static -ffree-line-length-0 -std=gnu -o check_krig.exe dof_kriging.o check_krig.f08 -Lgalib\lib -lfec -llapack -lrefblas
! cp check_krig.exe ..

! ※ -Wno-*はエラーメッセージ抑制オプション（省略化）
! ※ -ffree-line-length-0は1行の文字数制限解除オプション（省略化）

! ##############################################################################

program check_kriging
    use util
    use dof_kriging
    implicit none
    integer :: n_var

    n_var = 4

    call init_krig("kriging1.txt")

    call call_krig_grid(1, "ck_krig1-12.csv", [1, 2], 0.5d0)
    call call_krig_grid(1, "ck_krig1-13.csv", [1, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig1-14.csv", [1, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig1-23.csv", [2, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig1-24.csv", [2, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig1-34.csv", [3, 4], 0.5d0)

    call call_krig_grid(1, "ck_krig2-12.csv", [1, 2], 0.5d0)
    call call_krig_grid(1, "ck_krig2-13.csv", [1, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig2-14.csv", [1, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig2-23.csv", [2, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig2-24.csv", [2, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig2-34.csv", [3, 4], 0.5d0)

    call call_krig_grid(1, "ck_krig3-12.csv", [1, 2], 0.5d0)
    call call_krig_grid(1, "ck_krig3-13.csv", [1, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig3-14.csv", [1, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig3-23.csv", [2, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig3-24.csv", [2, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig3-34.csv", [3, 4], 0.5d0)

    call init_krig("kriging2.txt")

    call call_krig_grid(1, "ck_krig4-12.csv", [1, 2], 0.5d0)
    call call_krig_grid(1, "ck_krig4-13.csv", [1, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig4-14.csv", [1, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig4-23.csv", [2, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig4-24.csv", [2, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig4-34.csv", [3, 4], 0.5d0)

    call call_krig_grid(1, "ck_krig5-12.csv", [1, 2], 0.5d0)
    call call_krig_grid(1, "ck_krig5-13.csv", [1, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig5-14.csv", [1, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig5-23.csv", [2, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig5-24.csv", [2, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig5-34.csv", [3, 4], 0.5d0)

    call call_krig_grid(1, "ck_krig6-12.csv", [1, 2], 0.5d0)
    call call_krig_grid(1, "ck_krig6-13.csv", [1, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig6-14.csv", [1, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig6-23.csv", [2, 3], 0.5d0)
    call call_krig_grid(1, "ck_krig6-24.csv", [2, 4], 0.5d0)
    call call_krig_grid(1, "ck_krig6-34.csv", [3, 4], 0.5d0)

    !  n_var = 1

    ! call init_krig("kriging3.txt")

    ! call call_krig_grid(1, "ck_krig7.csv", [1, 0], 0.5d0)
    ! call call_krig_grid(1, "ck_krig8.csv", [1, 0], 0.5d0)
    ! call call_krig_grid(1, "ck_krig9.csv", [1, 0], 0.5d0)

    ! call init_krig("kriging4.txt")

    ! call call_krig_grid(1, "ck_krig10.csv", [1, 0], 0.5d0)
    ! call call_krig_grid(1, "ck_krig11.csv", [1, 0], 0.5d0)
    ! call call_krig_grid(1, "ck_krig12.csv", [1, 0], 0.5d0)


    contains


    subroutine call_krig_grid(model_no, output_file, idx, default)
        integer, intent(in) :: model_no, idx(2) ! 変動させる設計変数のインデックス
        character(*), intent(in) :: output_file
        real(8), intent(in) :: default ! 固定させる設計変数 [0-1]
        character(:), allocatable :: format

        real(8) :: dv(n_var), variables(n_var), z
        integer :: unit, i, j, nx, ny

        nx = 100
        ny = 100
        format = "(" // str(n_var) // "(f0.10',')f0.10)"

        open(newunit=unit, file=output_file, status="replace")
        do i = 0, nx
            if (mod(i, 10) == 0) print *, output_file, i

            do j = 0, ny
                dv = default
                dv(idx(1)) = dble(i) / nx
                dv(idx(2)) = dble(j) / ny

                variables = model(model_no)%dec(dv)
                z = model(model_no)%estimate(variables)
                write(unit, format) variables, z
            end do
        end do
        close(unit)
    end subroutine call_krig_grid
end program check_kriging
