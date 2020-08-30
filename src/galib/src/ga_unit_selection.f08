submodule (ga_unit) ga_unit_selection
  use interface
  ! use util
  implicit none

  contains


  ! ============================================================================
  ! constructor
  ! ============================================================================

  subroutine init_selection1(this, ftype, nsel, iparam)
    implicit none
    class(TSelection), intent(inout) :: this
    character(*), intent(in) :: ftype
    integer, intent(in) :: nsel
    integer, intent(in) :: iparam

    this%iparam = iparam
    call this%initialize(ftype, nsel)
  end subroutine init_selection1

  subroutine init_selection2(this, ftype, nsel, dparam)
    implicit none
    class(TSelection), intent(inout) :: this
    character(*), intent(in) :: ftype
    integer, intent(in) :: nsel
    real(8), intent(in) :: dparam

    this%dparam = dparam
    call this%initialize(ftype, nsel)
  end subroutine init_selection2

  subroutine init_selection3(this, ftype, nsel)
    implicit none
    class(TSelection), intent(inout) :: this
    character(*), intent(in) :: ftype
    integer, intent(in) :: nsel

    select case (ftype)
    case ("ROULETTE", "R")
      this%proc => selection_roulette
      this%iparam = 0
      this%dparam = 0d0
    case ("TOURNAMENT", "T")
      this%proc => selection_tournament
    case default
      write(0, "(3a/ a)") "Error: unknown keyword '", ftype, "'", &
                          "in subroutine 'init_selection'"
      call exit(1)
    end select

    this%num_selection = nsel
    this%ftype = ftype
  end subroutine init_selection3


  ! ============================================================================
  ! call procedure
  ! ============================================================================

  subroutine call_selection1x(this, probability, res)
    implicit none
    class(TSelection), intent(in) :: this
    real(8), intent(in) :: probability(:)
    integer, intent(out) :: res
    logical, allocatable :: rest(:)
    integer :: s, n

    s = size(probability)
    n = max(this%iparam, int(this%dparam * s))
    allocate(rest(s), source=.true.)

    if (n > 0 .and. n < s) then
      res = this%proc(max(probability, 0d0), shuffle(rest, n))
    else
      res = this%proc(max(probability, 0d0), rest)
    end if
  end subroutine call_selection1x

  subroutine call_selection1(this, probability, rest, res)
    implicit none
    class(TSelection), intent(in) :: this
    real(8), intent(in) :: probability(:)
    logical, intent(inout) :: rest(:)
    integer, intent(out) :: res
    integer :: s, n

    s = size(probability)
    n = max(this%iparam, int(this%dparam * s))
    if (n > 0 .and. n < s) then
      res = this%proc(max(probability, 0d0), shuffle(rest, n))
    else
      res = this%proc(max(probability, 0d0), rest)
    end if
    if (this%uniq) rest(res) = .false.
  end subroutine call_selection1

  subroutine call_selection2(this, probability, res)
    implicit none
    class(TSelection), intent(in) :: this
    real(8), intent(in) :: probability(:)
    integer, intent(out), allocatable :: res(:)
    logical, allocatable :: rest(:)
    integer :: s, n, i

    s = size(probability)
    n = max(this%iparam, int(this%dparam * s))
    allocate(rest(s), source=.true.)
    allocate(res(this%num_selection))
    do i = 1, this%num_selection
      if (n > 0 .and. n < s) then
        res(i) = this%proc(max(probability, 0d0), shuffle(rest, n))
      else
        res(i) = this%proc(max(probability, 0d0), rest)
      end if
      if (this%uniq) rest(res(i)) = .false.
    end do
  end subroutine call_selection2


  ! ============================================================================
  ! selection body
  ! ============================================================================

  integer function selection_roulette(probability, mask) result(res)
    implicit none
    real(8), intent(in) :: probability(:)
    logical, intent(in) :: mask(:)
    real(8) :: rand
    integer :: i

    if (all(probability <= 0d0 .or. .not. mask)) then
      res = selection_random(probability, mask)
      return
    end if

    rand = random() * sum(probability, mask=mask)
    do i = 1, size(probability)
      if (.not. mask(i)) cycle
      rand = rand - probability(i)
      if (rand < 0) then
        res = i
        return
      end if
    end do
    write(0, "(a/a)") "Error: selection failure", &
                      "in subroutine 'selection_roulette'"
    call exit(1)
    res = 0
  end function selection_roulette

  integer function selection_tournament(probability, mask) result(res)
    implicit none
    real(8), intent(in) :: probability(:)
    logical, intent(in) :: mask(:)

    res = maxloc(probability, dim=1, mask=mask)
  end function selection_tournament

  integer function selection_random(probability, mask) result(res)
    implicit none
    real(8), intent(in) :: probability(:)
    logical, intent(in) :: mask(:)
    integer :: rand, i

    rand = random(count(mask))
    do i = 1, size(probability)
      if (.not. mask(i)) cycle
      rand = rand - 1
      if (rand < 1) then
        res = i
        return
      end if
    end do
    write(0, "(a/a)") "Error: selection failure", &
                      "in subroutine 'selection_random'"
    call exit(1)
    res = 0
  end function selection_random


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_selection(this)
    implicit none
    type(TSelection), intent(inout) :: this

    if (associated(this%proc)) nullify(this%proc)
  end subroutine destroy_selection
end submodule ga_unit_selection
