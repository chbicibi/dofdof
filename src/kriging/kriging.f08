module mod_kriging
    use util
    use matrix
    implicit none

    private
    public kriging

    type :: kriging
        real(8), allocatable :: var(:, :)       ! 設計変数(変数番号, 個体番号)
        real(8), allocatable :: obj(:)          ! 目的関数(個体番号)
        real(8), allocatable :: theta(:)        ! 重み係数θ(変数番号)

        real(8), allocatable :: r_mat(:, :)     ! 相関行列R(個体番号, 個体番号)
        real(8), allocatable :: r_inv(:, :)     ! Rの逆行列
        real(8), allocatable :: dev_obj(:)      ! 目的関数の偏差(個体番号) (== f-μ)
        real(8), allocatable :: dev_mul(:)      ! (R^-1)・(f-μ)
        integer :: n_var                        ! 設計変数の数
        integer :: n_samples                    ! サンプル数

        real(8) :: r_det                        ! Rの行列式の対数値
        real(8) :: mu                           ! 平均μ
        real(8) :: sigma                        ! 分散σ^2
        real(8) :: ln                           ! 尤度Ln

        real(8), allocatable :: var_in(:, :)
        real(8), allocatable :: bounds(:, :), scale(:)
        logical :: feasible = .false.

        contains

        generic :: set_bounds => set_bounds1, set_bounds2
        procedure :: set_bounds1, set_bounds2

        generic :: set_scale => set_scale1, set_scale2
        procedure :: set_scale1, set_scale2

        procedure :: initialize
        procedure :: set_theta
        procedure :: calc_likelihood
        procedure :: estimate
        procedure :: enc, dec
        procedure :: error
        procedure :: check_error
        procedure :: show
        final :: destroy_instance
    end type kriging

    ! 定数
    ! real(8), parameter :: PI  = 4 * atan(1d0)

    contains


    ! ==========================================================================
    ! 初期化 (設計変数, 目的関数設定)
    ! ==========================================================================
    subroutine initialize(this, var, obj)
        class(kriging), intent(inout) :: this
        real(8), intent(in) :: var(:, :), obj(:)

        call destroy_instance(this)

        this%n_var = size(var, dim=1)
        this%n_samples = size(var, dim=2)

        this%var_in = var
        this%var = var
        this%obj = obj

        allocate(this%r_mat(this%n_samples, this%n_samples))
        allocate(this%r_inv(this%n_samples, this%n_samples))
        allocate(this%dev_obj(this%n_samples))
        allocate(this%dev_mul(this%n_samples))

        call this%set_bounds
        this%feasible = .false.
    end subroutine initialize


    ! ==========================================================================
    ! 設計範囲設定
    ! ==========================================================================
    subroutine set_bounds1(this)
        implicit none
        class(kriging), intent(inout) :: this

        call this%set_bounds(spread([0d0, 1d0], 2, this%n_var))
    end subroutine set_bounds1

    subroutine set_bounds2(this, bounds)
        implicit none
        class(kriging), intent(inout) :: this
        real(8), intent(in) :: bounds(:, :)
        integer :: i

        if (size(bounds, dim=1) /= 2) then
            print *, "設計範囲が不正です"
            call exit(1)
        end if

        this%bounds = bounds

        do i = 1, this%n_samples
            this%var(:, i) = this%enc(this%var_in(:, i))
        end do
    end subroutine set_bounds2


    ! ==========================================================================
    ! スケール設定
    ! ==========================================================================
    subroutine set_scale1(this, scale)
        implicit none
        class(kriging), intent(inout) :: this
        real(8), intent(in) :: scale

        this%scale = spread(scale, 1, this%n_var)
    end subroutine set_scale1

    subroutine set_scale2(this, scale)
        implicit none
        class(kriging), intent(inout) :: this
        real(8), intent(in) :: scale(:)

        this%scale = scale
    end subroutine set_scale2


    ! ==========================================================================
    ! 重み係数θ設定
    ! ==========================================================================
    subroutine set_theta(this, theta)
        implicit none
        class(kriging), intent(inout) :: this
        real(8), intent(in) :: theta(:)

        this%theta = theta * this%scale
        this%feasible = .false.
    end subroutine set_theta


    ! ==========================================================================
    ! 尤度計算
    ! ==========================================================================
    subroutine calc_likelihood(this, theta, show)
        implicit none
        class(kriging), intent(inout) :: this
        real(8), intent(in), optional :: theta(:)
        logical, intent(in), optional :: show
        character(:), allocatable :: format
        logical :: flag
        integer :: i, j

        if (present(theta)) call this%set_theta(theta)

        do i = 1, this%n_samples
            this%r_mat(i:, i) = [1d0, (corr(this%theta, this%var(:, i), &
                                            this%var(:, j)),            &
                                       j = i + 1, this%n_samples)]
        end do

        this%r_inv = sym_inverse_with_logdet(this%r_mat, this%r_det, flag)

        if (flag) then
            print "(a)", "Ln: infeasible (singular matrix)"
            this%ln = NaN
            this%feasible = .false.
            return
        end if

        call copy_l2u(this%r_inv)
        this%mu = dot_product(sum(this%r_inv, dim=1), this%obj) &
                            / sum(this%r_inv)
        this%dev_obj = this%obj - this%mu
        this%dev_mul = blas_sym_matmul(this%r_inv, this%dev_obj)
        this%sigma = dot_product(this%dev_obj, this%dev_mul) / this%n_samples

        if (this%sigma <= 0d0 .or. this%sigma > 10d0) then
            print "(a)", "Ln: infeasible (sigma <= 0)"
            this%ln = NaN
            this%feasible = .false.
            return
        end if

        this%ln = -0.5d0 * (this%n_samples * log(this%sigma) + this%r_det)
        this%feasible = .true.

        if (present(show)) then
            format = "('theta='" // str(this%n_var) // &
                     "f8.3' Ln='f8.2' det='f8.2' mu='es9.3e1' sig='es9.3e1)"
            print format, this%theta, this%ln, this%r_det, this%mu, this%sigma
        end if
    end subroutine calc_likelihood


    real(8) function estimate(this, x) result(res)
        implicit none
        class(kriging), intent(inout) :: this
        real(8), intent(in)  :: x(:)
        ! real(8) :: rmse, xi, phi, chi
        real(8), allocatable :: r_new(:)
        integer :: i

        if (.not. this%feasible) call this%calc_likelihood
        if (.not. this%feasible) then
            write(0, "(a)") "kriging%estimate: infeasible krigind model"
            call exit(1)
        end if

        r_new = [(corr(this%theta, this%enc(x), this%var(:, i)), &
                         i = 1, this%n_samples)]
        res = this%mu + dot_product(r_new, this%dev_mul)

        ! rmse = sqrt(sigma * (1d0 - dot_product(matmul(r_new, r_inv), r_new) &
        !        + (1d0 - dot_product(sum(r_inv, 1), r_new)) ** 2 / sum(r_inv)))
        ! xi = (minval(obj) - res) / rmse
        ! phi = 0.5d0 * erfc(-xi / sqrt(2d0))
        ! chi = 1d0 / sqrt(2d0 * PI) * exp(-0.5d0 * xi ** 2)
        ! ei = rmse * (xi * phi + chi)

        ! 確認
        ! print "('rmse =', f15.5)", rmse
        ! print "('xi   =', f15.5)", xi
        ! print "('phi  =', f15.5)", phi
        ! print "('chi  =', f15.5)", chi
        ! print "('ei   =', f15.5/)", ei
    end function estimate


    ! ==========================================================================
    ! 正規化
    ! ==========================================================================
    function enc(this, x) result(res)
        implicit none
        class(kriging), intent(in) :: this
        real(8), intent(in)  :: x(:)
        real(8), allocatable :: bounds_range(:), res(:)

        allocate(bounds_range(this%n_var), &
                 source=this%bounds(2, :) - this%bounds(1, :))
        where (bounds_range == 0) bounds_range = 1
        res = (x - this%bounds(1, :)) / bounds_range
    end function enc


    ! ==========================================================================
    ! 逆正規化
    ! ==========================================================================
    function dec(this, x) result(res)
        implicit none
        class(kriging), intent(in) :: this
        real(8), intent(in)  :: x(:)
        real(8), allocatable :: res(:)

        res = (this%bounds(2, :) - this%bounds(1, :)) * x + this%bounds(1, :)
    end function dec

    ! subroutine estimate(x, new_obj, ei)
    !   real(8), intent(in)  :: x(:)
    !   real(8), intent(out) :: new_obj, ei
    !   real(8) :: rmse, xi, phi, chi
    !   real(8), allocatable :: r_newvar(:)
    !   integer :: i

    !   r_newvar = [(corr(x, var(:, i)), i = 1, n_samples)]

    !   new_obj = mu + dot_product(matmul(r_newvar, r_inv), dev_obj)

    !   rmse = sqrt(sigma * (1d0 - dot_product(matmul(r_newvar, r_inv), r_newvar) &
    !          + (1d0 - dot_product(sum(r_inv, 1), r_newvar)) ** 2 / sum(r_inv)))
    !   xi = (minval(obj) - new_obj) / rmse
    !   phi = 0.5d0 * erfc(-xi / sqrt(2d0))
    !   chi = 1d0 / sqrt(2d0 * PI) * exp(-0.5d0 * xi ** 2)
    !   ei = rmse * (xi * phi + chi)

    !   ! 確認
    !   print "('rmse =', f15.5)", rmse
    !   print "('xi   =', f15.5)", xi
    !   print "('phi  =', f15.5)", phi
    !   print "('chi  =', f15.5)", chi
    !   print "('ei   =', f15.5/)", ei
    ! end subroutine estimate


    ! ==========================================================================
    ! inspection
    ! ==========================================================================

    real(8) function error(this)
        implicit none
        class(kriging), intent(inout) :: this
        integer :: i

        error = sum([((this%estimate(this%var_in(:, i)) - this%obj(i)) ** 2, &
                                    i = 1, this%n_samples)])
    end function error

    subroutine check_error(this)
        implicit none
        class(kriging), intent(inout) :: this
        real(8) :: range, obj, est, err, err_total(2)
        integer :: i

        err_total = 0d0
        range = maxval(this%obj) - minval(this%obj)
        do i = 1, this%n_samples
            obj = this%obj(i)
            est = this%estimate(this%var_in(:, i))
            err = est - obj
            err_total(1) = err_total(1) + abs(err)
            err_total(2) = err_total(2) + err ** 2
            ! print "(i5,3es13.3,es13.3' %')", i, obj, est, err, 100d0 * err / range
        end do
        print "('total'2es13.3,2(es13.3' %'))", err_total, &
                        100d0 * err_total(1) / this%n_samples / range, &
                        100d0 * sqrt(err_total(2) / this%n_samples) / range
    end subroutine check_error

    subroutine show(this)
        implicit none
        class(kriging), intent(inout) :: this
        character(:), allocatable :: format

        if (.not. this%feasible) call this%calc_likelihood

        format = "('theta='"//str(this%n_var)//"f10.5)"

        print "(/a)", "Kriging Model Data:"
        print format, this%theta
        print "('R=')"
        print "(5(5e12.3/))", transpose(this%r_mat(1:5, 1:5))
        print "('R^-1=')"
        print "(5(5e12.3/))", transpose(this%r_inv(1:5, 1:5))

        print "('detR  = ', g0.10)", this%r_det
        ! print "('lndetR= ', g0.10)", log(this%r_det)
        print "('mu    = ', g0.10)", this%mu
        print "('sigma = ', g0.10)", this%sigma
        print "('Ln    = ', g0.10/)", this%ln
    end subroutine show


    ! ==========================================================================
    ! destructor
    ! ==========================================================================

    subroutine destroy_instance(this)
        implicit none
        type(kriging), intent(inout) :: this

        if (allocated(this%var)) deallocate(this%var)
        if (allocated(this%obj)) deallocate(this%obj)
        if (allocated(this%theta)) deallocate(this%theta)
        if (allocated(this%r_mat)) deallocate(this%r_mat)
        if (allocated(this%r_inv)) deallocate(this%r_inv)
        if (allocated(this%dev_obj)) deallocate(this%dev_obj)
        if (allocated(this%dev_mul)) deallocate(this%dev_mul)
        if (allocated(this%var_in)) deallocate(this%var_in)
        if (allocated(this%bounds)) deallocate(this%bounds)
        if (allocated(this%scale)) deallocate(this%scale)
    end subroutine destroy_instance


    ! ==========================================================================
    ! 相関関数
    ! ==========================================================================
    real(8) function corr(theta, x1, x2)
        implicit none
        real(8), intent(in) :: theta(:), x1(:), x2(:)
        real(8) :: a

        a = dot_product(theta, (x1 - x2) ** 2)
        ! a = sum(abs(x1 - x2))
        if (a < 30) then
            corr = exp(-a)
            ! corr = 2d0 / (1d0 + exp(a))
            ! corr = max(1d0 - a * 0.05d0, 0d0)
        else
            corr = 0d0
        end if
    end function corr

end module mod_kriging
