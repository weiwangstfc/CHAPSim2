module boundary_conditions_mod


  public :: map_inlet_uprofile
  public :: Apply_BC_velocity
  public :: Apply_convective_outlet

contains

  
!==========================================================================================================
!> \brief Apply b.c. conditions 
!---------------------------------------------------------------------------------------------------------- 
!> Scope:  mpi    called-freq    xdomain     module
!>         all    once           specified   public
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[in]     d             domain
!> \param[out]    f             flow
!==========================================================================================================
  subroutine Apply_BC_velocity (dm, ux, uy, uz)
    use parameters_constant_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in )   :: dm
    real(WP), dimension(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3)), intent(inout) :: ux
    real(WP), dimension(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3)), intent(inout) :: uy
    real(WP), dimension(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3)), intent(inout) :: uz
    integer :: m, s
    type(DECOMP_INFO) :: dtmp


    ! if(dm%ibcx(1, 1) /= IBC_DIRICHLET .and. &
    !    dm%ibcx(2, 1) /= IBC_DIRICHLET .and. &
    !    dm%ibcy(1, 2) /= IBC_DIRICHLET .and. &
    !    dm%ibcy(2, 2) /= IBC_DIRICHLET .and. &
    !    dm%ibcz(1, 3) /= IBC_DIRICHLET .and. &
    !    dm%ibcz(2, 3) /= IBC_DIRICHLET ) return

!----------------------------------------------------------------------------------------------------------
!   ux at x-direction. BC of uy and uz at x-direction are given in C2P.
!----------------------------------------------------------------------------------------------------------
    m = 1
    dtmp = dm%dpcc
    ! constant value bc.
    do s = 1, 2
      if(dm%ibcx(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          ux(1, :, :) = dm%fbcx(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          ux(dtmp%xsz(m), :, :) = dm%fbcx(s, m)
        end if
      end if
    end do

    ! variation value bc
    if(dm%ibcx(1, 1) == IBC_UPROFILE .and. dtmp%xst(1) == 1) then
      do j = 1, dtmp%xsz(2)
        jj = dtmp%xst(2) + j - 1
        do k = 1, dtmp%xsz(3)
          kk = dtmp%xsz(3) + k - 1
          ux(1, j, k) = dm%fbcxinlet(jj, kk, 1)
        end do
      end do
    end if

!----------------------------------------------------------------------------------------------------------
!   uy at y-direction. BC of ux and uz at y-direction are given in C2P.
!----------------------------------------------------------------------------------------------------------
    m = 2
    dtmp = dm%dcpc
    do s = 1, 2
      if(dm%ibcy(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          uy(:, 1, :) = dm%fbcy(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          uy(:, dtmp%xsz(m), :) = dm%fbcy(s, m)
        end if
      end if
    end do

!----------------------------------------------------------------------------------------------------------
!   uz at z-direction. BC of ux and uy at z-direction are given in C2P.
!----------------------------------------------------------------------------------------------------------
    m = 3
    dtmp = dm%dccp
    do s = 1, 2
      if(dm%ibcz(s, m) == IBC_DIRICHLET) then
        if(dtmp%xst(m) == 1) then
          uz(:, :, 1) = dm%fbcz(s, m)
        end if
        if(dtmp%xen(m) == dm%np(m)) then
          uz(:, :, dtmp%xsz(m)) = dm%fbcz(s, m)
        end if
      end if
    end do

    return
  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine  map_inlet_uprofile(dm)
    use parameters_constant_mod 
    use udf_type_mod
    implicit none 
    type(t_domain), intent(inout) :: dm

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    character(len = 80) :: str
    real(WP) :: rtmp
    integer :: i, n
    real(WP), allocatable :: uprofile(:)
    real(WP), allocatable :: yy(:)
    real(WP), allocatable :: uxz(:, :)

    !----------------------------------------------------------------------------------------------------------
    ! to read given u-velocity profile, dimensionless, x/H, U/Umean
    !----------------------------------------------------------------------------------------------------------

    if(dm%ibcx(1, 1) /= IBC_UPROFILE) return

    open ( newunit = inputUnit,     &
           file    = "inlet_uprofile.dat", &
           status  = 'old',         &
           action  = 'read',        &
           iostat  = ioerr,         &
           iomsg   = iotxt)
    if(ioerr /= 0) then
      write (*, *) 'Problem openning : inlet_profile.dat.'
      write (*, *) 'Message: ', trim (iotxt)
      stop 4
    end if

    n = 0
    read(inputUnit, *, iostat = ioerr) str

    do
      read(inputUnit, *, iostat = ioerr) rtmp, rtmp
      if(ioerr /= 0) exit
      n = n + 1
    end do
    rewind(inputUnit)
    !----------------------------------------------------------------------------------------------------------
    ! to read given u-velocity profile, dimensionless, x/H, U/Umean
    !----------------------------------------------------------------------------------------------------------
    allocate ( uprofile (n) )
    allocate ( yy (n) )
    allocate ( uxz (dm%nc(2)) )

    read(inputUnit, *, iostat = ioerr) str
    do i = 1, n
      read(inputUnit, *, iostat = ioerr) yy(i), uprofile(i)
    end do
    close(inputUnit)
    !----------------------------------------------------------------------------------------------------------
    ! to map read profile to our case
    !----------------------------------------------------------------------------------------------------------
    allocate(dm%fbcxinlet(dm%nc(2), dm%nc(3), 5))
    do i = 1, 5
      dm%fbcxinlet(:, :, i) = dm%fbcx(1, i)
    end do

    call map_1d_profile_to_case(n, yy, uprofile, dm%nc(2), dm%yc(:), uxz(:))

    do k = 1, dm%nc(3)
      do j = 1, dm%nc(2)
        dm%fbcxinlet(j, k, 1) = uxz(j)
      end do
    end do

    deallocate(uprofile)
    deallocate(yy)
    deallocate(uxz)

  end subroutine

!==========================================================================================================
!==========================================================================================================
  subroutine Apply_convective_outlet


  end subroutine

end module
