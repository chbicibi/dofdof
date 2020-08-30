module module1
    implicit none

    integer :: i, a(3, 10)

    contains

    subroutine sub(n)
        integer, intent(in) :: n
        ! integer :: a(3)
        integer :: i

        a(:, n) = 1

        do i = 1, 10
            call p
        end do

        contains

        subroutine p
            print *, n, a(:, n)
            a(:, n) = a(:, n) + 1
            call sleep(1)
        end subroutine p
    end subroutine sub
end module module1


program main
    use module1
    implicit none

    !$omp parallel
    !$omp do
    do i = 1, 10
        call sub(i)
    end do
    !$omp end do
    !$omp end parallel

end program main