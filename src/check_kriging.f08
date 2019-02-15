
program check_kriging
    use dof_kriging
    implicit none

    call chech_kriging

    contains

    subroutine chech_kriging
        real(8) :: aoa, elv, ma, y1, y2, y3
        integer :: unit1, unit2, unit3, i, j

        ! call init_krig(num_var=3, pop_size=50)
        call load_krig(num_var=3)

        open(newunit=unit1, file="test1.csv", status="replace")
        open(newunit=unit2, file="test2.csv", status="replace")
        open(newunit=unit3, file="test3.csv", status="replace")

        ! aoa: [-5, 10]
        ! elv: [-40, 40]
        ! ma: [0.1, 0.5]

        do i = 0, 100
            if (mod(i, 10) == 0) print *, i
            do j = 0, 100
                aoa = 0.01 * i * 15 - 5
                elv = 0.01 * j * 80 - 40
                ma = 0.3d0

                y1 = estimate(1, aoa, elv, ma)
                y2 = estimate(2, aoa, elv, ma)
                y3 = estimate(3, aoa, elv, ma)

                write(unit1, "(3(f0.10','))") aoa, elv, y1
                write(unit2, "(3(f0.10','))") aoa, elv, y2
                write(unit3, "(3(f0.10','))") aoa, elv, y3
            end do
        end do
        close(unit3)
        close(unit2)
        close(unit1)
    end subroutine chech_kriging
end program check_kriging
