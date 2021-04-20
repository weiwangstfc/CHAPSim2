!##############################################################################
module domain_decomposition_mod
  implicit none
  type pencil_t
    integer :: ipsz
    integer :: jpsz
    integer :: kpsz
    integer :: iprange(2)
    integer :: jprange(2)
    integer :: kprange(2)
  contains
    private
    procedure :: Print_debug
    generic :: Print => Print_debug
    generic :: write(formatted) => Print_debug
  end type pencil_t
  
  type(pencil_t), save :: local0_xpencil
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

  subroutine Initialize_domain_decompsition (d)
    use mpi_mod, only: nrow, ncol
    use input_general_mod
    use decomp_2d
    implicit none
    type(t_domain), intent(in)   :: d

    call decomp_2d_init( d%np(1), d%np(2), d%np(3), nrow, ncol, is_periodic(:) )
    
    local0_xpencil%iprange(1) = xstart(1)
    local0_xpencil%jprange(1) = xstart(2)
    local0_xpencil%kprange(1) = xstart(3)

    local0_xpencil%iprange(2) = xend(1)
    local0_xpencil%jprange(2) = xend(2)
    local0_xpencil%kprange(2) = xend(3)

    local0_xpencil%ipsz = xsize(1)
    local0_xpencil%jpsz = xsize(2)
    local0_xpencil%kpsz = xsize(3)

    write(*, '(dt)') local_xpencil

  end subroutine Initialize_domain_decompsition

end module domain_decomposition_mod