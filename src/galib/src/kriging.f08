module kriging
  use util
  use matrix
  use problem
  implicit none

  private
  public TKriging

  type, extends(TProblem) :: TKriging
    real(8), allocatable :: variables(:, :) ! 設計変数(変数番号, 個体番号)
    real(8), allocatable :: objectives(:)   ! 目的関数(個体番号)
    real(8), allocatable :: theta(:)        ! 重み係数θ(変数番号)

    real(8), allocatable :: r_mat(:, :)     ! 相関行列R(個体番号, 個体番号)
    real(8), allocatable :: r_inv(:, :)     ! Rの逆行列
    real(8), allocatable :: dev_obj(:)      ! 目的関数の偏差(個体番号) (== f-μ)
    real(8), allocatable :: dev_mul(:)      ! (R^-1)・(f-μ)
    integer :: num_var                      ! 設計変数の数
    integer :: pop_size                     ! 個体数

    real(8) :: r_det                        ! Rの行列式の対数値
    real(8) :: mu                           ! 平均μ
    real(8) :: sigma                        ! 分散σ^2
    real(8) :: ln                           ! 尤度Ln

    real(8), allocatable :: var_raw(:, :)
    real(8), allocatable :: low_bounds(:), upp_bounds(:), scale(:)
    logical :: feasible = .false.
    ! logical :: ready = .false.

    contains

    generic :: initialize => initialize_krig
    generic :: set_bounds => set_bounds1, set_bounds2
    generic :: set_scale => set_scale1, set_scale2

    procedure :: initialize_krig
    procedure :: set_bounds1, set_bounds2
    procedure :: set_scale1, set_scale2
    procedure :: set_theta
    procedure :: call_func_so_con
    procedure :: calc_likelihood
    procedure :: estimate
    procedure :: enc, dec
    procedure :: error
    procedure :: check_error
    procedure :: show
    final :: destroy_instance
  end type TKriging

  real(8), parameter :: PI  = 4d0 * atan(1d0)

  contains

  subroutine initialize_krig(this, variables, objectives)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in) :: variables(:, :), objectives(:)

    call destroy_instance(this)

    this%num_var = size(variables, dim=1)
    this%pop_size = size(variables, dim=2)

    this%var_raw = variables
    this%variables = variables
    this%objectives = objectives

    allocate(this%r_mat(this%pop_size, this%pop_size))
    allocate(this%r_inv(this%pop_size, this%pop_size))
    allocate(this%dev_obj(this%pop_size))
    allocate(this%dev_mul(this%pop_size))

    call this%set_bounds
    this%feasible = .false.
  end subroutine initialize_krig

  ! ============================================================================

  subroutine set_bounds1(this)
    implicit none
    class(TKriging), intent(inout) :: this

    this%low_bounds = filled(this%num_var, 0d0)
    this%upp_bounds = filled(this%num_var, 1d0)
  end subroutine set_bounds1

  subroutine set_bounds2(this, lower, upper)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in) :: lower(:), upper(:)
    integer :: i

    this%low_bounds = lower
    this%upp_bounds = upper

    do i = 1, this%pop_size
      this%variables(:, i) = this%enc(this%var_raw(:, i))
    end do
  end subroutine set_bounds2

  subroutine set_scale1(this, scale)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in) :: scale

    this%scale = filled(this%num_var, scale)
  end subroutine set_scale1

  subroutine set_scale2(this, scale)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in) :: scale(:)

    this%scale = scale
  end subroutine set_scale2

  subroutine set_theta(this, theta)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in) :: theta(:)

    this%theta = theta * this%scale
    this%feasible = .false.
  end subroutine set_theta


  ! ============================================================================
  ! target functions for optimization (minimize)
  ! ============================================================================

  ! subroutine call_func_so(this, variables, objective)
  !   implicit none
  !   class(TKriging) :: this
  !   real(8), intent(in) :: variables(:)
  !   real(8), intent(out) :: objective

  !   this%theta = variables * this%scale
  !   call this%calc_likelihood
  !   objective = -this%ln
  ! end subroutine call_func_so

  subroutine call_func_so_con(this, variables, objective, feasible)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out) :: objective
    logical, intent(out) :: feasible

    call this%calc_likelihood(variables)

    feasible = this%feasible
    ! if (feasible) then
    !   objective = 1d0 / (1d0 + 0.01d0 * this%ln)
    ! else
    !   objective = 1d0
    ! end if
    objective = -this%ln
    ! if (feasible) then
    !   objective = this%error()
    ! else
    !   objective = 1d0
    ! end if
  end subroutine call_func_so_con

  ! ============================================================================

  subroutine calc_likelihood(this, theta)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in), optional :: theta(:)
    character(:), allocatable :: format
    logical :: flag
    integer :: i, j

    if (present(theta)) call this%set_theta(theta)

    do i = 1, this%pop_size
      this%r_mat(i:, i) = [1d0, (corr(this%theta, this%variables(:, i),  &
                                                  this%variables(:, j)), &
                                 j = i + 1, this%pop_size)]
    end do
    ! call copy_l2u(this%r_mat)

    ! this%r_inv = inverse_with_logdet(this%r_mat, this%r_det, flag)
    this%r_inv = sym_inverse_with_logdet(this%r_mat, this%r_det, flag)

    if (flag) then
      print "(a)", "Ln: infeasible (singular matrix)"
      this%ln = -1e30
      this%feasible = .false.
      return
    end if

    call copy_l2u(this%r_inv)
    this%mu = dot_product(sum(this%r_inv, dim=1), this%objectives) &
              / sum(this%r_inv)
    this%dev_obj = this%objectives - this%mu
    this%dev_mul = blas_sym_matmul(this%r_inv, this%dev_obj)
    this%sigma = dot_product(this%dev_obj, this%dev_mul) / this%pop_size

    if (this%sigma <= 0d0 .or. this%sigma > 10d0) then
      print "(a)", "Ln: infeasible (sigma <= 0)"
      this%ln = -1e30
      this%feasible = .false.
      return
    end if

    this%ln = -0.5d0 * (this%pop_size * log(this%sigma) + this%r_det)
    this%feasible = .true.
    format = "('theta='" // str(this%num_var) // &
             "f8.3' Ln='f8.2' det='f8.2' mu='es9.3e1' sig='es9.3e1)"
    print format, this%theta, this%ln, this%r_det, this%mu, this%sigma
  end subroutine calc_likelihood

  real(8) function estimate(this, x) result(res)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8), intent(in)  :: x(:)
    ! real(8) :: rmse, xi, phi, chi
    real(8), allocatable :: r_new(:)
    integer :: i

    if (.not. this%feasible) call this%calc_likelihood
    if (.not. this%feasible) then
      write(0, "(a)") "TKriging%estimate: infeasible krigind model"
      call exit(1)
    end if

    r_new = [(corr(this%theta, this%enc(x), this%variables(:, i)), &
             i = 1, this%pop_size)]
    res = this%mu + dot_product(r_new, this%dev_mul)

    ! rmse = sqrt(sigma * (1d0 - dot_product(matmul(r_new, r_inv), r_new) &
    !        + (1d0 - dot_product(sum(r_inv, 1), r_new)) ** 2 / sum(r_inv)))
    ! xi = (minval(objectives) - res) / rmse
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


  ! ============================================================================

  function enc(this, x) result(res)
    implicit none
    class(TKriging), intent(in) :: this
    real(8), intent(in)  :: x(:)
    real(8), allocatable :: res(:)

    res = (x - this%low_bounds) / (this%upp_bounds - this%low_bounds)
  end function enc

  function dec(this, x) result(res)
    implicit none
    class(TKriging), intent(in) :: this
    real(8), intent(in)  :: x(:)
    real(8), allocatable :: res(:)

    res = x * (this%upp_bounds - this%low_bounds) + this%low_bounds
  end function dec

  ! subroutine estimate(x, new_obj, ei)
  !   real(8), intent(in)  :: x(:)
  !   real(8), intent(out) :: new_obj, ei
  !   real(8) :: rmse, xi, phi, chi
  !   real(8), allocatable :: r_newvar(:)
  !   integer :: i

  !   r_newvar = [(corr(x, variables(:, i)), i = 1, pop_size)]

  !   new_obj = mu + dot_product(matmul(r_newvar, r_inv), dev_obj)

  !   rmse = sqrt(sigma * (1d0 - dot_product(matmul(r_newvar, r_inv), r_newvar) &
  !          + (1d0 - dot_product(sum(r_inv, 1), r_newvar)) ** 2 / sum(r_inv)))
  !   xi = (minval(objectives) - new_obj) / rmse
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


  ! ============================================================================
  ! inspection
  ! ============================================================================

  real(8) function error(this)
    implicit none
    class(TKriging), intent(inout) :: this
    integer :: i

    error = sum([((this%estimate(this%var_raw(:, i)) - this%objectives(i)) ** 2, &
                  i = 1, this%pop_size)])
  end function error

  subroutine check_error(this)
    implicit none
    class(TKriging), intent(inout) :: this
    real(8) :: range, obj, est, err, err_total(2)
    integer :: i

    err_total = 0d0
    range = maxval(this%objectives) - minval(this%objectives)
    do i = 1, this%pop_size
      obj = this%objectives(i)
      est = this%estimate(this%var_raw(:, i))
      err = est - obj
      err_total(1) = err_total(1) + abs(err)
      err_total(2) = err_total(2) + err ** 2
      ! print "(i5,3es13.3,es13.3' %')", i, obj, est, err, 100d0 * err / range
    end do
    print "('total'2es13.3,2(es13.3' %'))", err_total, &
            100d0 * err_total(1) / this%pop_size / range, &
            100d0 * sqrt(err_total(2) / this%pop_size) / range
  end subroutine check_error

  subroutine show(this)
    implicit none
    class(TKriging), intent(inout) :: this
    character(:), allocatable :: format

    if (.not. this%feasible) call this%calc_likelihood

    format = "('theta='"//str(this%num_var)//"f10.5)"

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


  ! ============================================================================
  ! destructor
  ! ============================================================================

  subroutine destroy_instance(this)
    implicit none
    type(TKriging), intent(inout) :: this

    if (allocated(this%variables)) deallocate(this%variables)
    if (allocated(this%objectives)) deallocate(this%objectives)
    if (allocated(this%theta)) deallocate(this%theta)
    if (allocated(this%r_mat)) deallocate(this%r_mat)
    if (allocated(this%r_inv)) deallocate(this%r_inv)
    if (allocated(this%dev_obj)) deallocate(this%dev_obj)
    if (allocated(this%dev_mul)) deallocate(this%dev_mul)
    if (allocated(this%var_raw)) deallocate(this%var_raw)
    if (allocated(this%low_bounds)) deallocate(this%low_bounds)
    if (allocated(this%upp_bounds)) deallocate(this%upp_bounds)
    if (allocated(this%scale)) deallocate(this%scale)
  end subroutine destroy_instance


  ! ============================================================================
  ! auxiliary functions
  ! ============================================================================

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
end module kriging
