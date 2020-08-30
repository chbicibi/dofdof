submodule (ga_unit) ga_unit_crossover
  use interface
  ! use util
  implicit none

  contains


  ! ============================================================================
  ! constructor
  ! ============================================================================

  subroutine init_crossover(this, ftype, rate, dparam)
    implicit none
    class(TCrossover), intent(inout) :: this
    character(*), intent(in) :: ftype
    real(8), intent(in) :: rate, dparam

    select case (ftype)
    case ("SBX")
      this%dproc => crossover_sbx
    case ("BLX")
      this%dproc => crossover_blx
    case ("ORDER")
      this%iproc => crossover_order
    case default
      write(0, "(3a/a)") "Error: unknown keyword '", ftype, "'", &
                         "in subroutine 'init_crossover'"
      call exit(1)
    end select

    call this%set_params(ftype, rate, dparam)
  end subroutine init_crossover


  ! ============================================================================
  ! call procedure
  ! ============================================================================

  subroutine call_crossover1(this, parents, children)
    implicit none
    class(TCrossover), intent(in) :: this
    real(8), intent(in) :: parents(:, :)
    real(8), intent(out), allocatable :: children(:, :)

    if (random() <= this%rate) then
      children = this%dproc(this%dparam, parents)
    else
      children = parents
    end if
  end subroutine call_crossover1

  subroutine call_crossover2(this, parents, child)
    implicit none
    class(TCrossover), intent(in) :: this
    real(8), intent(in) :: parents(:, :)
    real(8), intent(out), allocatable :: child(:)

    if (random() <= this%rate) then
      associate (children => this%dproc(this%dparam, parents))
        child = children(:, 1)
      end associate
    else
      child = parents(:, 1)
    end if
  end subroutine call_crossover2

  subroutine call_icrossover1(this, parents, children)
    implicit none
    class(TCrossover), intent(in) :: this
    integer, intent(in) :: parents(:, :)
    integer, intent(out), allocatable :: children(:, :)

    if (random() <= this%rate) then
      children = this%iproc(this%dparam, parents)
    else
      children = parents
    end if
  end subroutine call_icrossover1

  subroutine call_icrossover2(this, parents, child)
    implicit none
    class(TCrossover), intent(in) :: this
    integer, intent(in) :: parents(:, :)
    integer, intent(out), allocatable :: child(:)

    if (random() <= this%rate) then
      associate (children => this%iproc(this%dparam, parents))
        child = children(:, 1)
      end associate
    else
      child = parents(:, 1)
    end if
  end subroutine call_icrossover2


  ! ============================================================================
  ! body
  ! ============================================================================

  function crossover_sbx(eta_c, parents) result(children)
    implicit none
    real(8), intent(in) :: eta_c, parents(:, :)
    real(8), allocatable :: children(:, :)
    real(8), allocatable :: x1(:), x2(:), sum(:), diff(:), u(:)
    real(8), allocatable :: alpha(:, :), beta(:, :), betaq(:, :)
    integer :: s

    s = size(parents, dim=1)
    allocate(x1(s), source=min(parents(:, 1), parents(:, 2)))
    allocate(x2(s), source=max(parents(:, 1), parents(:, 2)))
    allocate(u(s), source=random_array(s))
    allocate(sum(s), source= x2 + x1)
    allocate(diff(s), source= x2 - x1)

    where (diff == 0d0) diff = 1d0
    beta = 1d0 + 2d0 / reshape([x1 / diff, (1d0 - x2) / diff], [s, 2])
    alpha = 2d0 - beta ** (-eta_c - 1d0)
    betaq = spread(u, 2, 2)
    where (betaq <= 1d0 / alpha)
      betaq = (betaq * alpha) ** (1d0 / (eta_c + 1d0))
    else where
      betaq = (2d0 - betaq * alpha) ** (-1d0 / (eta_c + 1d0))
    end where
    children = reshape(0.5d0 * [sum - betaq(:, 1) * diff, sum + betaq(:, 2) * diff], [s, 2])
  end function crossover_sbx

  function crossover_blx(alpha, parents) result(children)
    implicit none
    real(8), intent(in) :: alpha, parents(:, :)
    real(8), allocatable :: children(:, :)
    integer :: len

    len = size(parents, dim=1)
    children = reshape(parents(:, 1) + ((1d0 + 2d0 * alpha) * random_array(len) - alpha) * &
                                       (parents(:, 2) - parents(:, 1)), [len, 1])
  end function crossover_blx

  function crossover_order(ddum, parents) result(children)
    implicit none
    real(8), intent(in) :: ddum
    integer, intent(in) :: parents(:, :)
    integer, allocatable :: children(:, :)
    integer :: s, n, pair(2)

    if (is_disjoint(parents(:, 1), parents(:, 2))) then
      children = crossover_1pt(ddum, parents)
      return
    end if

    s = size(parents, dim=1)

    n = random(s-3) + 1
    allocate(children(s, 2))
    children(:, 1) = parents(:, 1)
    children(:, 2) = parents(:, 1)
    call rearrange(children(1:n, 1), parents(:, 2))
    call rearrange(children(n+1:s, 2), parents(:, 2))
    return

    pair = shuffle(s-1, 2) + random(2) - 1
    pair = pair(sort(pair))
    children = parents
    call rearrange(children(pair(1):pair(2), 1), parents(:, 2))
    call rearrange(children(pair(1):pair(2), 2), parents(:, 1))

    ! print *, pair
    ! print "(a, 10i3)", "P1", parents(:, 1)
    ! print "(a, 10i3)", "P2", parents(:, 2)
    ! print "(a, 10i3)", "C1", children(:, 1)
    ! print "(a, 10i3)", "C2", children(:, 2)
    ! stop
  end function crossover_order

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

  function crossover_1pt(ddum, parents) result(children)
    implicit none
    real(8), intent(in) :: ddum
    integer, intent(in) :: parents(:, :)
    integer, allocatable :: children(:, :)
    integer :: s, n

    ! print *, "crossover_1pt"

    s = size(parents, dim=1)
    n = random(s-1)

    children = parents
    children(1:n, 1) = parents(1:n, 2)
    children(1:n, 2) = parents(1:n, 1)
    ! print *, "end crossover_1pt"
  end function crossover_1pt


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_crossover(this)
    implicit none
    type(TCrossover), intent(inout) :: this

    if (associated(this%dproc)) nullify(this%dproc)
    if (associated(this%iproc)) nullify(this%iproc)
  end subroutine destroy_crossover
end submodule ga_unit_crossover
