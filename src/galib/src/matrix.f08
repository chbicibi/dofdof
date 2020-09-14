module matrix
    use util, only: str
    implicit none

    private
    public :: inverse, sym_inverse, inverse_with_logdet, sym_inverse_with_logdet
    public :: blas_dot, blas_norm, blas_matmul, blas_sym_matmul
    public :: copy_l2u

    real(8), allocatable :: work(:)
    integer :: lwork = 0

    contains

    function inverse(a, det, flag) result (a_inv)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(out), optional :: det
        logical, intent(out), optional :: flag
        real(8), allocatable :: a_inv(:, :), lu(:, :)
        integer, allocatable :: ipiv(:)
        integer :: info, i

        call lapack_inverse(a, a_inv, lu, ipiv, info)
        if (present(det)) then
            if (info /= 0) then
                det = 0d0
            else
                det = abs(product([(lu(ipiv(i), i), i = 1, size(a, dim=1))]))
            end if
        end if
        if (present(flag)) flag = info /= 0
    end function inverse

    function sym_inverse(a, det, flag) result (a_inv)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(out), optional :: det
        logical, intent(out), optional :: flag
        real(8), allocatable :: a_inv(:, :), ld(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n, info, i

        n = size(a, dim=1)
        call lapack_sym_inverse(a, a_inv, ld, ipiv, info)
        if (present(det)) then
            if (info /= 0) then
                det = 0d0
            else
                det = 1d0
                i = 1
                do while (i <= n)
                    if (ipiv(i) > 0) then
                        det = det * ld(i, i)
                        i = i + 1
                    else
                        det = det * (ld(i, i) * ld(i+1, i+1) - ld(i+1, i) ** 2)
                        i = i + 2
                    end if
                end do
                det = abs(det)
            end if
        end if
        if (present(flag)) flag = info /= 0
    end function sym_inverse

    function inverse_with_logdet(a, det, flag) result (a_inv)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(out) :: det
        logical, intent(out), optional :: flag
        real(8), allocatable :: a_inv(:, :), lu(:, :)
        integer, allocatable :: ipiv(:)
        integer :: info, i

        ! call naive_inverse(a, a_inv, lu, ipiv, info)
        call lapack_inverse(a, a_inv, lu, ipiv, info)
        if (info /= 0) then
            det = -huge(0d0)
        else
            det = sum([(log(abs(lu(i, i))), i = 1, size(a, dim=1))])
        end if
        if (present(flag)) flag = info /= 0
    end function inverse_with_logdet

    function sym_inverse_with_logdet(a, det, flag) result (a_inv)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(out) :: det
        logical, intent(out), optional :: flag
        real(8), allocatable :: a_inv(:, :), ld(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n, info, i

        n = size(a, dim=1)
        call lapack_sym_inverse(a, a_inv, ld, ipiv, info)
        if (info /= 0) then
            det = -huge(0d0)
        else
            det = 0d0
            i = 1
            do while (i <= n)
                if (ipiv(i) > 0) then
                    det = det + log(abs(ld(i, i)))
                    i = i + 1
                else
                    det = det + log(abs(ld(i, i) * ld(i+1, i+1) - ld(i+1, i) ** 2))
                    i = i + 2
                end if
            end do
        end if
        if (present(flag)) flag = info /= 0
    end function sym_inverse_with_logdet

    ! ============================================================================

    subroutine naive_inverse(a, a_inv, lu, ipiv, info)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(out), allocatable :: a_inv(:, :), lu(:, :)
        integer, intent(out), allocatable :: ipiv(:)
        integer, intent(out) :: info
        real(8), allocatable :: b(:)
        integer :: n, i

        n = size(a, dim=1)
        lu = a
        allocate(a_inv(n, n), b(n))
        call naive_ludecomposition(lu, ipiv, info)
        if (info /= 0) return
        do i = 1, n
            b = 0d0
            b(i) = 1d0
            call naive_lusolution(lu, b, ipiv, a_inv(:, i), info)
        end do
    end subroutine naive_inverse

    subroutine naive_ludecomposition(a, ipiv, info)
        implicit none
        real(8), intent(inout) :: a(:, :)
        integer, intent(out), allocatable :: ipiv(:)
        integer, intent(out) :: info
        real(8) :: dd
        integer :: imax, sign, n, i, j, k

        n = size(a, dim=1)
        ipiv = [(i, i = 1, n)]
        sign = 1
        info = 0
        do j = 1, n - 1
            imax = maxloc([(abs(a(ipiv(i), j)), i = j, n)], dim=1) + j - 1
            if (imax /= j) then
                k = ipiv(j)
                ipiv(j) = ipiv(imax)
                ipiv(imax) = k
                sign = -sign
            end if
            do i = j + 1, n
                if (a(ipiv(j), j) == 0d0) then
                    info = j
                    return
                end if
                dd = a(ipiv(i), j) / a(ipiv(j), j)
                a(ipiv(i), j+1:n) = a(ipiv(i), j+1:n) - dd * a(ipiv(j), j+1:n)
                a(ipiv(i), j) = dd
            end do
        end do
    end subroutine naive_ludecomposition

    subroutine naive_lusolution(a, b, ipiv, x, info)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(inout) :: b(:)
        integer, intent(in) :: ipiv(:)
        real(8), intent(out) :: x(:)
        integer, intent(out) :: info
        real(8) :: dd
        integer :: n, i, j

        n = size(a, dim=1)
        info = 0
        do j = 1, n - 1
            do i = j + 1, n
                b(ipiv(i)) = b(ipiv(i)) - a(ipiv(i), j) * b(ipiv(j))
            end do
        end do
        x(n) = b(ipiv(n)) / a(ipiv(n), n)
        do i = n - 1, 1, -1
            dd = b(ipiv(i))
            do j = i + 1, n
                dd = dd - a(ipiv(i), j) * x(j)
            end do
            x(i) = dd / a(ipiv(i), i)
        end do
    end subroutine naive_lusolution


    ! ============================================================================
    ! LAPACK / BLAS
    ! ============================================================================

    ! ******************************************************************
    ! LAPACK
    ! ******************************************************************

    subroutine lapack_inverse(a, a_inv, lu, ipiv, info)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(out), allocatable :: a_inv(:, :), lu(:, :)
        integer, intent(out), allocatable :: ipiv(:)
        integer, intent(out) :: info
        integer :: n

        n = size(a, dim=1)
        allocate(ipiv(n), a_inv(n, n))
        lu = a

        call lapack_alloc_workspace(n)
        call dgetrf(n, n, lu, n, ipiv, info)
        if (info /= 0) return

        a_inv = lu
        call dgetri(n, a_inv, n, ipiv, work, lwork, info)
        lwork = size(work)
    end subroutine lapack_inverse

    subroutine lapack_sym_inverse(a, a_inv, ld, ipiv, info)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(out), allocatable :: a_inv(:, :), ld(:, :)
        integer, intent(out), allocatable :: ipiv(:)
        integer, intent(out) :: info
        integer :: n

        n = size(a, dim=1)
        allocate(ipiv(n), a_inv(n, n))
        ld = a

        call lapack_alloc_workspace(n)
        call dsytrf('L', n, ld, n, ipiv, work, lwork, info)
        if (info /= 0) return
        a_inv = ld
        call dsytri('L', n, a_inv, n, ipiv, work, lwork, info)
        lwork = size(work)
    end subroutine lapack_sym_inverse

    subroutine lapack_alloc_workspace(n)
        implicit none
        integer, intent(in) :: n

        if (lwork < 64 * n) then
            if (allocated(work)) deallocate(work)
            lwork = 64 * n
            allocate(work(lwork))
        end if
    end subroutine lapack_alloc_workspace


    ! ******************************************************************
    ! BLAS
    ! ******************************************************************

    real(8) function blas_dot(a, b) result(res)
        implicit none
        real(8), intent(in) :: a(:), b(:)

        call ddot(3, a, 1, b, 1, res)
    end function blas_dot

    real(8) function blas_norm(a) result(res)
        implicit none
        real(8), intent(in) :: a(:)

        call dnrm2(size(a), a, 1, res)
    end function blas_norm

    function blas_matmul(a, b) result(c)
        implicit none
        real(8), intent(in) :: a(:, :), b(:, :)
        real(8), allocatable :: c(:, :)
        integer :: m, n, k

        m = size(a, dim=1)
        n = size(b, dim=2)
        k = size(b, dim=1)
        allocate(c(m, n))
        call dgemm('N', 'N', m, n, k, 1d0, a, m, b, k, 0d0, c, m)
    end function blas_matmul

    function blas_sym_matmul(a, b) result(res)
        implicit none
        real(8), intent(in) :: a(:, :), b(:)
        real(8), allocatable :: res(:)
        integer :: n

        n = size(a, dim=1)
        allocate(res(n))
        call dsymv('L', n, 1d0, a, n, b, 1, 0d0, res, 1)
    end function blas_sym_matmul

    ! ============================================================================

    subroutine copy_l2u(a)
        implicit none
        real(8) :: a(:, :)
        integer :: i

        do i = 1, size(a, dim=1) - 1
            a(i, i+1:) = a(i+1:, i)
        end do
    end subroutine copy_l2u

    subroutine dump_mat(a, file)
        use util
        implicit none
        real(8), intent(in) :: a(:, :)
        character(*), intent(in) :: file
        character(:), allocatable :: format
        integer :: unit, i

        format = "(" // str(size(a, dim=2) - 1) // "(g0.5',')g0.5)"
        open(newunit=unit, file=file, status='replace')
            do i = 1, size(a, dim=1)
                write(unit, format) a(i, :)
            end do
        close(unit)
    end subroutine dump_mat
end module matrix
