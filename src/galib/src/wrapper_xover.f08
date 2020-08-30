module wrapper
use iso_c_binding
use util
implicit none

contains

subroutine order_crossover(iargs1, iargs2, n1, result1, result2) bind(c, name="order_crossover")
    implicit none
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(in) :: iargs1(n1), iargs2(n1)
    integer(c_int), intent(out) :: result1(n1), result2(n1)
    integer :: s, n

    if (is_disjoint(iargs1, iargs2)) then
      call crossover_1pt(iargs1, iargs2, result1, result2)
      return
    end if

    s = n1
    n = random(s-3) + 1

    result1 = iargs1
    result2 = iargs2

    call rearrange(result1(1:n), iargs2)
    call rearrange(result2(n+1:s), iargs1)

    ! pair = shuffle(s-1, 2) + random(2) - 1
    ! pair = pair(sort(pair))
    ! children = parents
    ! call rearrange(children(pair(1):pair(2), 1), parents(:, 2))
    ! call rearrange(children(pair(1):pair(2), 2), parents(:, 1))
end subroutine order_crossover

  logical function is_disjoint(array1, array2) result(flag)
    implicit none
    integer, intent(in) :: array1(:), array2(:)
    integer :: i, j

    flag = .false.
    do j = 1, size(array1)
      do i = 1, size(array2)
        if (array1(i) == array2(j)) return
      end do
    end do
    flag = .true.
  end function is_disjoint

  subroutine rearrange(array, other)
    implicit none
    integer, intent(inout) :: array(:)
    integer, intent(in) :: other(:)
    integer, allocatable :: index(:), order(:)
    integer :: s1, s2, i

    s1 = size(array)
    s2 = size(other)
    index = pack(integers(s1), [(any(array(i) == other), i = 1, s1)])
    order = pack(other, [(any(other(i) == array(index)), i = 1, s2)])

    if (size(index) /= size(order)) then
      write(0, "(a,2i5)") "Error: size(index) /= size(order)", size(index), size(order)
      call exit(1)
    end if

    array(index) = order
  end subroutine rearrange

  subroutine crossover_1pt(parent1, parent2, child1, child2)
    implicit none
    integer, intent(in) :: parent1(:), parent2(:)
    integer, intent(out) :: child1(:), child2(:)
    integer :: s, n

    ! print *, "crossover_1pt"

    s = size(parent1)
    n = random(s-1)

    child1 = parent1
    child2 = parent2
    child1(1:n) = parent2(1:n)
    child2(1:n) = parent1(1:n)
    ! print *, "end crossover_1pt"
  end subroutine crossover_1pt

end module wrapper