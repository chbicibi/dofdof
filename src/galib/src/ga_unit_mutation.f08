submodule (ga_unit) ga_unit_mutation
  use interface
  implicit none

  contains


  ! ============================================================================
  ! constructor
  ! ============================================================================

  subroutine init_mutation(this, ftype, rate, dparam)
    implicit none
    class(TMutation), intent(inout) :: this
    character(*), intent(in) :: ftype
    real(8), intent(in) :: rate, dparam

    select case (ftype)
    case ("POLYNOMIAL", "PM")
      this%dproc => mutation_pm
    case ("SWAP")
      this%iproc => mutation_swap
    case ("REVERSE")
      this%iproc => mutation_reverse
    case ("SHUFFLE")
      this%iproc => mutation_shuffle
    case default
      write(0, "(3a/a)") "Error: unknown keyword '", ftype, "'", &
                         "in subroutine 'init_mutation'"
      call exit(1)
    end select

    call this%set_params(ftype, rate, dparam)
  end subroutine init_mutation


  ! ============================================================================
  ! call procedure
  ! ============================================================================

  subroutine call_mutation(this, variables)
    implicit none
    class(TMutation), intent(in) :: this
    real(8), intent(inout) :: variables(:)

    if (random() <= this%rate) variables = this%dproc(this%dparam, variables)
  end subroutine call_mutation

  subroutine call_imutation(this, variables)
    implicit none
    class(TMutation), intent(in) :: this
    integer, intent(inout) :: variables(:)

    if (random() <= this%rate) variables = this%iproc(this%dparam, variables)
  end subroutine call_imutation


  ! ============================================================================
  ! mutation body
  ! ============================================================================

  function mutation_pm(eta_m, indiv) result(new_indiv)
    implicit none
    real(8), intent(in) :: eta_m, indiv(:)
    real(8), allocatable :: new_indiv(:)
    real(8), allocatable :: u(:), delta(:)
    integer :: s

    s = size(indiv)
    allocate(u(s), source=random_array(s))
    allocate(delta(s), source=u)
    where (u <= 0.5d0)
      delta = ((2d0 * u) ** (1d0 / (1d0 + eta_m)) - 1d0) * indiv
    else where
      delta = ((2d0 - 2d0 * u) ** (1d0 / (1d0 + eta_m)) - 1d0) * (indiv - 1d0)
    end where
    new_indiv = indiv + delta
  end function mutation_pm

  function mutation_swap(eta_m, indiv) result(new_indiv)
    implicit none
    real(8), intent(in) :: eta_m
    integer, intent(in) :: indiv(:)
    integer, allocatable :: new_indiv(:)
    integer :: index(2)

    index = shuffle(size(indiv), 2)
    new_indiv = indiv
    new_indiv(index) = reverse(indiv(index))
  end function mutation_swap

  function mutation_reverse(eta_m, indiv) result(new_indiv)
    implicit none
    real(8), intent(in) :: eta_m
    integer, intent(in) :: indiv(:)
    integer, allocatable :: new_indiv(:)
    integer :: index(2)

    index = shuffle(size(indiv), 2)
    new_indiv = indiv
    new_indiv(index(1):index(2)) = reverse(indiv(index(1):index(2)))
  end function mutation_reverse

  function mutation_shuffle(eta_m, indiv) result(new_indiv)
    implicit none
    real(8), intent(in) :: eta_m
    integer, intent(in) :: indiv(:)
    integer, allocatable :: new_indiv(:)

    new_indiv = shuffle(indiv)
  end function mutation_shuffle


  ! ============================================================================
  ! destructor
  ! ============================================================================

  elemental subroutine destroy_mutation(this)
    implicit none
    type(TMutation), intent(inout) :: this

    if (associated(this%dproc)) nullify(this%dproc)
    if (associated(this%iproc)) nullify(this%iproc)
  end subroutine destroy_mutation
end submodule ga_unit_mutation
