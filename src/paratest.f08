!


program main
    use util
    implicit none

    integer, parameter :: N = 9
    integer i, j, a(N)

    j = 1

    !$omp parallel
    !$omp do
    do i=1,N
        a(i) = i
        call printn(i)
    end do
    !$omp end do
    !$omp end parallel

    print *, a

    contains

    subroutine printn(a)
        implicit none
        integer, intent(in) :: a
        real(8) :: r

        r = random() * 10

        call sleep(int(r))

        print *, j, a, int(r)
        j = j + 1

    end subroutine printn
end program main