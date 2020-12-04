!##############################################################################
module domain_decomposition_mod
  implicit none
  type pencil_t
    integer :: isz
    integer :: jsz
    integer :: ksz
    integer :: irange(2)
    integer :: jrange(2)
    integer :: krange(2)
  contains
    private
    procedure :: Print_debug
    generic :: Print => Print_debug
    generic :: write(formatted) => Print_debug
  end type pencil_t

  type(pencil_t) :: xpencil

  private
  public :: Initialize_domain_decompsition

contains

  !--------------------------------
  subroutine Print_debug(this, unit, iotype, v_list, iostat, iomsg)
    use iso_fortran_env, only : error_unit
    use mpi_mod
    class(pencil_t), intent(in) :: this
    integer, intent(in) :: unit
    character(len = *), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(len = *), intent(inout) :: iomsg
    integer :: i_pass

    iostat = 0
    iomsg = ""
    
    this_block: do i_pass = 1, 1
      !write(unit, *, iostat = iostat, iomsg = iomsg) 'thermalProperty'
      !if(iostat /= 0) exit this_block
      if(iotype(1:2) == 'DT' .and. len(iotype) > 2) &
        write(unit, *, iostat = iostat, iomsg = iomsg) iotype(3:)
      if(iostat /= 0) exit this_block

      write(unit, '(I6,A,3I7)', iostat = iostat, iomsg = iomsg) &
      myid, ":x-index", this%irange(1:2), this%isz
      write(unit, '(I6,A,3I7)', iostat = iostat, iomsg = iomsg) &
      myid, ":y-index", this%jrange(1:2), this%jsz
      write(unit, '(I6,A,3I7)', iostat = iostat, iomsg = iomsg) &
      myid, ":z-index", this%krange(1:2), this%ksz
      
      if(iostat /= 0) exit this_block
    end do this_block

    if(iostat /= 0) then
      write (error_unit, "(A)") "print error : " // trim(iomsg)
      write (error_unit, "(A, I0)") "  iostat : ", iostat
    end if
  end subroutine Print_debug

  subroutine Initialize_domain_decompsition ()
    use mpi_mod
    use input_general_mod
    use decomp_2d

    call decomp_2d_init( npx, npy, npz, p_row, p_col, is_periodic(:) )
    
    xpencil%irange(1) = xstart(1)
    xpencil%jrange(1) = xstart(2)
    xpencil%krange(1) = xstart(3)

    xpencil%irange(2) = xend(1)
    xpencil%jrange(2) = xend(2)
    xpencil%krange(2) = xend(3)

    xpencil%isz = xsize(1)
    xpencil%jsz = xsize(2)
    xpencil%ksz = xsize(3)

    write(*, '(dt)') xpencil

  end subroutine Initialize_domain_decompsition

end module domain_decomposition_mod



!