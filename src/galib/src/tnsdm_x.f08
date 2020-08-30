!!! NOT use !!!

submodule (tnsdm) tnsdm_1
  implicit none

  contains


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine save_result(this, filename, mode1, mode2)
    implicit none
    class(TTNSDM), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: mode1, mode2
    character(:), allocatable :: format
    type(string), allocatable :: s(:)
    integer :: unit, i

    format = "(" // str(this%num_obj+this%num_var+this%num_con) // "(es15.8',')i0,2(','es15.8))"

    open(newunit=unit, file=filename)
      s = [(string("obj" // str(i)), i = 1, this%num_obj), &
           (string("var" // str(i)), i = 1, this%num_var), &
           (string("con" // str(i)), i = 1, this%num_con), &
           string("rank"), string("fitness"), string("crowding-dist")]
      write(unit, "(a)") join(s, ",")

      do i = 1, this%pop_size
        if ((mode1 == 0 .or. xor(mode1 == 2, this%population(i)%rank == 1)) .and. &
            (mode2 == 0 .or. xor(mode2 == 2, this%population(i)%indiv%feasible))) &
              write(unit, format) this%population(i)%indiv%objectives,            &
                                  this%dec(this%population(i)%indiv%dvariables),   &
                                  this%population(i)%indiv%constraints,           &
                                  this%population(i)%rank,                        &
                                  this%population(i)%fitness,                     &
                                  this%population(i)%crowding
      end do
    close(unit)
  end subroutine save_result

  subroutine save_history(this, filename, mode1, mode2)
    implicit none
    class(TTNSDM), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(in) :: mode1, mode2
    character(:), allocatable :: format
    type(string), allocatable :: s(:)
    integer :: unit, i, j

    format = "(i0," // str(this%num_obj+this%num_var+this%num_con) // "(','es15.8)"

    open(newunit=unit, file=filename)
      s = [string("step"), &
           (string("obj" // str(i)), i = 1, this%num_obj), &
           (string("var" // str(i)), i = 1, this%num_var), &
           (string("con" // str(i)), i = 1, this%num_con), &
           string("rank"), string("fitness"), string("crowding-dist")]
      write(unit, "(a)") join(s, ",")

      do j = 1, size(this%history, dim=2)
        do i = 1, this%pop_size
          if ((mode1 == 0 .or. xor(mode1 == 2, this%history(i, j)%rank == 1)) .and. &
              (mode2 == 0 .or. xor(mode2 == 2, this%history(i, j)%indiv%feasible))) &
                write(unit, format) j, &
                                    this%history(i, j)%indiv%objectives,            &
                                    this%dec(this%history(i, j)%indiv%dvariables),  &
                                    this%history(i, j)%indiv%constraints,           &
                                    this%history(i, j)%rank,                        &
                                    this%history(i, j)%fitness,                     &
                                    this%history(i, j)%crowding
        end do
      end do
    close(unit)
  end subroutine save_history
end submodule tnsdm_1
