module moead_unit
  use interface
  use util
  implicit none

  private
  public :: get_scalar_func, bi_theta

  real(8) :: bi_theta = 5d0

  contains


  function get_scalar_func(ftype) result(res)
    implicit none
    character(*), intent(in) :: ftype
    procedure(func_1d_1d_1d_d), pointer :: res

    select case (ftype)
    case ("WEIGHTEDSUM", "WS")
      res => weighted_sum
    case ("CHEBYSHEV", "TCHEBYCHEFF", "CH", "TCH")
      res => chebyshev1
    case ("BOUNDARYINTERSECTION", "BI")
      res => boundaryintersection
    case default
      write(0, "(3a/a)") "Error: unknown keyword '", ftype, "'", &
                         "in function 'get_scalar_func'"
      call exit(1)
    end select
  end function get_scalar_func


  ! ============================================================================
  ! private
  ! ============================================================================

  real(8) function weighted_sum(objectives, weight, reference_point) result(res)
    implicit none
    real(8), intent(in) :: objectives(:), weight(:), reference_point(:)

    res = sum(weight * abs(objectives - reference_point))
  end function weighted_sum

  real(8) function chebyshev1(objectives, weight, reference_point) result(res)
    implicit none
    real(8), intent(in) :: objectives(:), weight(:), reference_point(:)

    res = maxval(weight * abs(objectives - reference_point))
  end function chebyshev1

  real(8) function boundaryintersection(objectives, weight, reference_point) result(res)
    implicit none
    real(8), intent(in) :: objectives(:), weight(:), reference_point(:)
    real(8) :: d1, d2

    d1 = abs(dot_product((objectives - reference_point), normalize(weight)))
    d2 = norm(objectives - (reference_point - d1 * weight))

    res = d1 + bi_theta * d2
  end function boundaryintersection
end module moead_unit
