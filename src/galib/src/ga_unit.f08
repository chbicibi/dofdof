module ga_unit
  use interface
  use util
  implicit none

  private
  public :: TSelection, TCrossover, TMutation
  public :: rank_pareto, crowding_distance

  type :: TOperator
    character(:), allocatable :: ftype
    real(8) :: rate = 1d0
    integer :: iparam = 0
    real(8) :: dparam = 0d0

    contains

    generic :: set_params => set_params1, set_params2
    procedure :: set_params1
    procedure :: set_params2
  end type TOperator

  type, extends(TOperator) :: TSelection
    integer :: num_selection
    logical :: uniq = .true.
    procedure(func_1d_1l_i), nopass, pointer :: proc

    contains

    generic :: initialize => init_selection1, init_selection2, init_selection3
    generic :: call => call_selection1, call_selection2
    procedure :: init_selection1, init_selection2, init_selection3
    procedure :: call_selection1, call_selection2
    final :: destroy_selection
  end type TSelection

  type, extends(TOperator) :: TCrossover
    procedure(func_d_2i_2i), nopass, pointer :: iproc
    procedure(func_d_2d_2d), nopass, pointer :: dproc

    contains

    generic :: call => call_crossover1, call_crossover2, call_icrossover1, call_icrossover2
    procedure :: initialize => init_crossover
    procedure :: call_crossover1, call_crossover2, call_icrossover1, call_icrossover2
    final :: destroy_crossover
  end type TCrossover

  type, extends(TOperator) :: TMutation
    procedure(func_d_1i_1i), nopass, pointer :: iproc
    procedure(func_d_1d_1d), nopass, pointer :: dproc

    contains

    generic :: call => call_mutation, call_imutation
    procedure :: initialize => init_mutation
    procedure :: call_mutation, call_imutation
    final :: destroy_mutation
  end type TMutation

  ! type, extends(TCrossover) :: TICrossover
  !   procedure(func_d_2d_2d), nopass, pointer :: iproc

  !   contains

  !   generic :: call => call_icrossover1, call_icrossover2
  !   ! procedure :: initialize => init_crossover
  !   procedure :: call_icrossover1, call_icrossover2
  !   ! final :: destroy_crossover
  ! end type TICrossover

  ! type, extends(TMutation) :: TIMutation
  !   procedure(func_d_1d_1d), nopass, pointer :: iproc

  !   contains

  !   ! procedure :: initialize => init_mutation
  !   procedure :: call => call_imutation
  !   ! final :: destroy_mutation
  ! end type TIMutation


  ! ******************************************************************
  ! workspace
  ! ******************************************************************

  ! logical, allocatable :: ws_dominated(:, :)

  ! interface check_uniq
  !   procedure :: dcheck_uniq, icheck_uniq
  ! end interface check_uniq

  interface rank_pareto
    procedure :: rank_pareto_base, rank_pareto_proc, rank_pareto_div
    procedure :: rank_pareto_con1_proc, rank_pareto_con2_proc
  end interface rank_pareto

  interface
    module subroutine init_selection1(this, ftype, nsel, iparam)
      class(TSelection), intent(inout) :: this
      character(*), intent(in) :: ftype
      integer, intent(in) :: nsel
      integer, intent(in) :: iparam
    end subroutine init_selection1

    module subroutine init_selection2(this, ftype, nsel, dparam)
      class(TSelection), intent(inout) :: this
      character(*), intent(in) :: ftype
      integer, intent(in) :: nsel
      real(8), intent(in) :: dparam
    end subroutine init_selection2

    module subroutine init_selection3(this, ftype, nsel)
      class(TSelection), intent(inout) :: this
      character(*), intent(in) :: ftype
      integer, intent(in) :: nsel
    end subroutine init_selection3

    module subroutine call_selection1(this, probability, rest, res)
      class(TSelection), intent(in) :: this
      real(8), intent(in) :: probability(:)
      logical, intent(inout) :: rest(:)
      integer, intent(out) :: res
    end subroutine call_selection1

    module subroutine call_selection2(this, probability, res)
      class(TSelection), intent(in) :: this
      real(8), intent(in) :: probability(:)
      integer, intent(out), allocatable :: res(:)
    end subroutine call_selection2

    elemental module subroutine destroy_selection(this)
      type(TSelection), intent(inout) :: this
    end subroutine destroy_selection
  end interface

  interface
    module subroutine init_crossover(this, ftype, rate, dparam)
      class(TCrossover), intent(inout) :: this
      character(*), intent(in) :: ftype
      real(8), intent(in) :: rate, dparam
    end subroutine init_crossover

    module subroutine call_crossover1(this, parents, children)
      class(TCrossover), intent(in) :: this
      real(8), intent(in) :: parents(:, :)
      real(8), intent(out), allocatable :: children(:, :)
    end subroutine call_crossover1

    module subroutine call_crossover2(this, parents, child)
      class(TCrossover), intent(in) :: this
      real(8), intent(in) :: parents(:, :)
      real(8), intent(out), allocatable :: child(:)
    end subroutine call_crossover2

    module subroutine call_icrossover1(this, parents, children)
      class(TCrossover), intent(in) :: this
      integer, intent(in) :: parents(:, :)
      integer, intent(out), allocatable :: children(:, :)
    end subroutine call_icrossover1

    module subroutine call_icrossover2(this, parents, child)
      class(TCrossover), intent(in) :: this
      integer, intent(in) :: parents(:, :)
      integer, intent(out), allocatable :: child(:)
    end subroutine call_icrossover2

    elemental module subroutine destroy_crossover(this)
      type(TCrossover), intent(inout) :: this
    end subroutine destroy_crossover
  end interface

  interface
    module subroutine init_mutation(this, ftype, rate, dparam)
      class(TMutation), intent(inout) :: this
      character(*), intent(in) :: ftype
      real(8), intent(in) :: rate, dparam
    end subroutine init_mutation

    module subroutine call_mutation(this, variables)
      class(TMutation), intent(in) :: this
      real(8), intent(inout) :: variables(:)
    end subroutine call_mutation

    module subroutine call_imutation(this, variables)
      class(TMutation), intent(in) :: this
      integer, intent(inout) :: variables(:)
    end subroutine call_imutation

    elemental module subroutine destroy_mutation(this)
      type(TMutation), intent(inout) :: this
    end subroutine destroy_mutation
  end interface

  ! interface
  !   function crossover_order(ddum, parents) result(children)
  !     real(8), intent(in) :: ddum
  !     integer, intent(in) :: parents(:, :)
  !   end function crossover_order
  ! end interface

  contains


  ! ============================================================================
  ! setter
  ! ============================================================================

  subroutine set_params1(this, ftype, dparam)
    implicit none
    class(TOperator), intent(inout) :: this
    character(*), intent(in) :: ftype
    real(8), intent(in) :: dparam

    this%ftype = ftype
    this%dparam = dparam
  end subroutine set_params1

  subroutine set_params2(this, ftype, rate, dparam)
    implicit none
    class(TOperator), intent(inout) :: this
    character(*), intent(in) :: ftype
    real(8), intent(in) :: rate, dparam

    this%ftype = ftype
    this%rate = rate
    this%dparam = dparam
  end subroutine set_params2


  ! ============================================================================
  ! check uniq
  ! ============================================================================

  ! logical function dchech_uniq(var, sample) result(flag)
  !   implicit none
  !   real(8), intent(in) :: var(:), sample(:, :)

  !   do i = 1, size(sample, dim=2)
  !     if (all(var == sample(:, i))) then
  !       flag = .false.
  !       return
  !     end if
  !   end do
  !   flag = .true.
  ! end function dchech_uniq

  ! logical function ichech_uniq(var, sample) result(flag)
  !   implicit none
  !   integer, intent(in) :: var(:), sample(:, :)

  !   do i = 1, size(sample, dim=2)
  !     if (all(var == sample(:, i))) then
  !       flag = .false.
  !       return
  !     end if
  !   end do
  !   flag = .true.
  ! end function ichech_uniq


  ! ============================================================================
  ! ranking method, etc.
  ! ============================================================================

  function rank_pareto_base(val, mode) result(rank)
    implicit none
    real(8), intent(in) :: val(:, :)
    integer, intent(in), optional :: mode
    integer, allocatable :: rank(:)
    integer :: n, i, j

    n = size(val, dim=2)
    if (present(mode) .and. mode > 0) then
      rank = rank_pareto_core([((all(val(:, i) > val(:, j)), j = 1, n), i = 1, n)], n)
    else
      rank = rank_pareto_core([((all(val(:, i) >= val(:, j)) .and. &
                                 any(val(:, i) /= val(:, j)), j = 1, n), i = 1, n)], n)
    end if
  end function rank_pareto_base

  function rank_pareto_proc(val, proc) result(rank)
    implicit none
    real(8), intent(in) :: val(:, :)
    procedure(func_1d_1d_l) :: proc
    integer, allocatable :: rank(:)
    integer :: n, i, j

    n = size(val, dim=2)
    rank = rank_pareto_core([((proc(val(:, i), val(:, j)), j = 1, n), i = 1, n)], n)
  end function rank_pareto_proc

  function rank_pareto_div(val, div, mode) result(rank)
    implicit none
    real(8), intent(in) :: val(:, :)
    integer, intent(in) :: div(:)
    integer, intent(in), optional :: mode
    integer, allocatable :: rank(:), index(:)
    logical, allocatable :: mask(:)
    integer :: i

    allocate(rank(size(val, dim=2)))
    do i = 1, maxval(div)
      mask = div == i
      index = mask_index(mask)
      rank(index) = rank_pareto(val(:, index), mode)
    end do
  end function rank_pareto_div

  function rank_pareto_con1_proc(val, con, mask, proc) result(rank)
    implicit none
    real(8), intent(in) :: val(:, :), con(:)
    logical, intent(in) :: mask(:)
    procedure(func_2d_1d_1l_l) :: proc
    integer, allocatable :: rank(:)
    integer :: n, i, j

    n = size(val, dim=2)
    rank = rank_pareto_core([((proc(val(:, [i, j]), con([i, j]), mask([i, j])), &
                               j = 1, n), i = 1, n)], n)
  end function rank_pareto_con1_proc

  function rank_pareto_con2_proc(val, con, mask, proc) result(rank)
    implicit none
    real(8), intent(in) :: val(:, :), con(:, :)
    logical, intent(in) :: mask(:)
    procedure(func_2d_2d_1l_l) :: proc
    integer, allocatable :: rank(:)
    integer :: n, i, j

    n = size(val, dim=2)
    rank = rank_pareto_core([((proc(val(:, [i, j]), con(:, [i, j]), mask([i, j])), &
                               j = 1, n), i = 1, n)], n)
  end function rank_pareto_con2_proc

  function rank_pareto_core(dominated, n) result(rank)
    implicit none
    logical, intent(in) :: dominated(n, n)
    integer, intent(in) :: n
    integer, allocatable :: rank(:)
    integer, allocatable :: num_dominated(:)
    logical, allocatable :: front_mask(:)
    integer :: rank_no, i

    rank = filled(n, 0)
    ! num_dominated = [(count(dominated(n*i+1:n*(i+1))), i = 0, n - 1)]
    allocate(num_dominated(n), source=[(count(dominated(:, i)), i = 1, n)])

    do rank_no = 1, n
      front_mask = rank == 0 .and. num_dominated == 0
      where (front_mask) rank = rank_no
      if (all(rank > 0)) exit
      do i = 1, n
        ! if (rank(i) == 0) num_dominated(i) = num_dominated(i) - &
        !                     count(front_mask .and. dominated(n*(i-1)+1:n*i))
        if (rank(i) == 0) num_dominated(i) = num_dominated(i) - &
                            count(front_mask .and. dominated(:, i))
      end do
    end do
  end function rank_pareto_core

  ! ============================================================================

  function crowding_distance(val, rank) result(distance)
    implicit none
    real(8), intent(in) :: val(:, :)
    integer, intent(in) :: rank(:)
    real(8), allocatable :: distance(:)
    integer, allocatable :: order(:)
    integer :: len, s, i, j

    len = size(rank)
    allocate(distance(len))
    do i = 1, maxval(rank)
      order = sort(val(1, :), pack(integers(len), rank == i))
      s = size(order)
      select case (s)
      case (1)
        distance(order) = huge(0d0)
      case (2:)
        distance(order) = [huge(0d0),                                           &
                           ((val(1, order(j + 1)) - val(1, order(j - 1)) +      &
                             val(2, order(j - 1)) - val(2, order(j + 1))) *     &
                            safe_inv(abs(val(1, order(j)) - val(2, order(j)))), &
                            j = 2, s - 1), huge(0d0)]
      end select
    end do
  end function crowding_distance
end module ga_unit
