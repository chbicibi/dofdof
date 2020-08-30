module interface
  implicit none

  ! func_<n><t><o>_...
  ! <n> => dimension
  ! - none => scalar
  ! - n > 0 => array
  ! <t> => type
  ! - i => integer
  ! - d => real(8) (== double precision)
  ! <o> => intent option
  ! - o => intent(inout)
  ! - none => intent(in)

  abstract interface
    function func_i_d_d(a, b) result(res)
      integer, intent(in) :: a
      real(8), intent(in) :: b
      real(8), allocatable :: res(:)
    end function func_i_d_d

    function func_1i_1i(a) result(res)
      real(8), intent(in) :: a(:)
      real(8), allocatable :: res(:)
    end function func_1i_1i

    function func_1d_1d(a) result(res)
      real(8), intent(in) :: a(:)
      real(8), allocatable :: res(:)
    end function func_1d_1d

    logical function func_2d_l(a) result(res)
      real(8), intent(in) :: a(:, :)
    end function func_2d_l

    function func_d_1i_1i(a, b) result(res)
      real(8), intent(in) :: a
      integer, intent(in) :: b(:)
      integer, allocatable :: res(:)
    end function func_d_1i_1i

    function func_d_1d_1d(a, b) result(res)
      real(8), intent(in) :: a, b(:)
      real(8), allocatable :: res(:)
    end function func_d_1d_1d

    function func_d_2i_2i(a, b) result(res)
      real(8), intent(in) :: a
      integer, intent(in) :: b(:, :)
      integer, allocatable :: res(:, :)
    end function func_d_2i_2i

    function func_d_2d_2d(a, b) result(res)
      real(8), intent(in) :: a, b(:, :)
      real(8), allocatable :: res(:, :)
    end function func_d_2d_2d

    logical function func_1d_1d_l(a, b) result(res)
      real(8), intent(in) :: a(:), b(:)
    end function func_1d_1d_l

    function func_1d_1d_1d(a, b) result(res)
      real(8), intent(in) :: a(:), b(:)
      real(8), allocatable :: res(:)
    end function func_1d_1d_1d

    integer function func_1d_1l_i(a, b)
      real(8), intent(in) :: a(:)
      logical, intent(in) :: b(:)
    end function func_1d_1l_i

    real(8) function func_1d_1d_1d_d(a, b, c)
      real(8), intent(in) :: a(:), b(:), c(:)
    end function func_1d_1d_1d_d

    logical function func_2d_1d_1l_l(a, b, c) result(res)
      real(8), intent(in) :: a(:, :), b(:)
      logical, intent(in) :: c(:)
    end function func_2d_1d_1l_l

    logical function func_2d_2d_1l_l(a, b, c) result(res)
      real(8), intent(in) :: a(:, :), b(:, :)
      logical, intent(in) :: c(:)
    end function func_2d_2d_1l_l
  end interface
end module interface
