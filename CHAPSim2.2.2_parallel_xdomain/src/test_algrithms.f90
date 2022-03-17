module test_algorithms_mod
  use operations

  public   :: Test_schemes
  private  :: Test_interpolation
  !private  :: Test_1st_derivative
  !private  :: Test_2nd_derivative
  
contains
!===============================================================================
!===============================================================================
!> \brief In-code independent test code for algorithms and schemes
!>
!> This subroutine is only called in the main program for testing.
!> Please select the test options which you are interested in.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
subroutine Test_schemes()
  use vars_df_mod
  implicit none

  ! Please use below information for input file
  ! x = 0, 2pi
  ! y = 0, 2pi
  ! z = 0, 2pi

  !call Test_TDMA_cyclic
  !call Test_TDMA_noncyclic
  call Test_interpolation (domain(1))
  !call Test_1st_derivative(flow, domain)
  !call Test_2nd_derivative(flow, domain)
  return 
end subroutine 

!===============================================================================
!===============================================================================
!> \brief To test this subroutine for mid-point interpolation.
!>
!> This subroutine is called in \ref Test_schemes. Define the logicals to choose
!> which test section is required. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
  subroutine Test_interpolation(dm)
    use operations
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: dm
    integer :: i, j, k
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    real(WP) :: err, err_Linf, err_L2
    integer  :: wrt_unit

    real(WP) :: den_ccc (dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) )
    real(WP) :: den_ppp (dm%dppp%xsz(1), dm%dppp%xsz(2), dm%dppp%xsz(3) )

    real(WP) :: den_xpcc(dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) )
    real(WP) :: den_ycpc(dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) )
    real(WP) :: den_zccp(dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) )

    real(WP) :: den_xcpp(dm%dcpp%xsz(1), dm%dcpp%xsz(2), dm%dcpp%xsz(3) )
    real(WP) :: den_ypcp(dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) )
    real(WP) :: den_zppc(dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) )

    open (newunit = wrt_unit, file = 'check_test_algorithms.dat', position="append")

    ! initialise a 3d variable on cell-centre
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, dm%nc(2)
        yc = dm%yc(j)
        do i = 1, dm%nc(1)
          xc = dm%h(1) * (real(i - 1, WP) + HALF)
          den_ccc(i, j, k) = sin_wp ( xc / FOUR) + sin_wp(yc / FOUR) + sin_wp(zc / FOUR)
        end do
      end do
    end do

    ! c2p in x
    call Get_x_midp_C2P_3D(den_ccc, den_xpcc, dm, dm%ibcx(:, 5), dm%fbcx(:, 5) )
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, dm%nc(2)
        yc = dm%yc(j)
        do i = 1, dm%np(1)
          xp = dm%h(1) * real(i - 1, WP)
          ref = sin_wp ( xp  / FOUR) + sin_wp(yc / FOUR) + sin_wp(zc / FOUR)
          err = dabs(den_xpcc(i, j, k) - ref)
          if(err > err_Linf) err_Linf = err
          err_L2 = err_L2 + err**2
        end do
      end do
    end do
    err_L2 = dsqrt ( err_L2 / real(dm%nc(3) * dm%nc(2) * dm%np(1), WP) )
    write (wrt_unit,'(A, 1I5.1, 2ES13.5)') '# interp_c2p_x', dm%np(1), err_Linf, err_L2

    ! c2p in y
    call Get_y_midp_C2P_3D(den_ccc, den_ycpc, dm, dm%ibcy(:, 5), dm%fbcy(:, 5) )
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, dm%np(2)
        yp = dm%yp(j)
        do i = 1, dm%nc(1)
          xc = dm%h(1) * (real(i - 1, WP) + HALF)
          ref = sin_wp ( xc / FOUR ) + sin_wp(yp / FOUR) + sin_wp(zc / FOUR)
          err = dabs(den_ycpc(i, j, k) - ref)
          if(err > err_Linf) err_Linf = err
          err_L2 = err_L2 + err**2
        end do
      end do
    end do
    err_L2 = dsqrt ( err_L2 / real(dm%nc(3) * dm%np(2) * dm%nc(1), WP) )
    write (wrt_unit,'(A, 1I5.1, 2ES13.5)') '# interp_c2p_y', dm%np(2), err_Linf, err_L2

    ! c2p in z
    call Get_z_midp_C2P_3D(den_ccc, den_zccp, dm, dm%ibcz(:, 5), dm%fbcz(:, 5) )
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      do j = 1, dm%nc(2)
        yc = dm%yc(j)
        do i = 1, dm%nc(1)
          xc = dm%h(1) * (real(i - 1, WP) + HALF)
          ref = sin_wp ( xc / FOUR ) + sin_wp(yc / FOUR) + sin_wp(zp / FOUR)
          err = dabs(den_zccp(i, j, k) - ref)
          if(err > err_Linf) err_Linf = err
          err_L2 = err_L2 + err**2
        end do
      end do
    end do
    err_L2 = dsqrt ( err_L2 / real(dm%np(3) * dm%nc(2) * dm%nc(1), WP) )
    write (wrt_unit,'(A, 1I5.1, 2ES13.5)') '# interp_c2p_z', dm%np(3), err_Linf, err_L2

    ! initialise a 3d variable on point based
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      do j = 1, dm%np(2)
        yp = dm%yp(j)
        do i = 1, dm%np(1)
          xp = dm%h(1) * real(i - 1, WP)
          den_ppp(i, j, k) = sin_wp ( xp / FOUR ) + sin_wp(yp / FOUR) + sin_wp(zp / FOUR)
        end do
      end do
    end do

    ! p2c in x
    call Get_x_midp_P2C_3D(den_ppp, den_xcpp, dm, dm%ibcx(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      do j = 1, dm%np(2)
        yp = dm%yp(j)
        do i = 1, dm%nc(1)
          xc = dm%h(1) * (real(i - 1, WP) + HALF)
          ref = sin_wp ( xc / FOUR ) + sin_wp(yp / FOUR) + sin_wp(zp / FOUR)
          err = dabs(den_xcpp(i, j, k) - ref)
          if(err > err_Linf) err_Linf = err
          err_L2 = err_L2 + err**2
        end do
      end do
    end do
    err_L2 = dsqrt ( err_L2 / real(dm%np(3) * dm%np(2) * dm%nc(1), WP) )
    write (wrt_unit,'(A, 1I5.1, 2ES13.5)') '# interp_p2c_x', dm%np(1), err_Linf, err_L2

    ! p2c in y
    call Get_y_midp_P2C_3D(den_ppp, den_ypcp, dm, dm%ibcy(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%np(3)
      zp = dm%h(3) * real(k - 1, WP)
      do j = 1, dm%nc(2)
        yc = dm%yc(j)
        do i = 1, dm%np(1)
          xp = dm%h(1) * real(i - 1, WP)
          ref = sin_wp ( xp / FOUR ) + sin_wp(yc / FOUR) + sin_wp(zp / FOUR)
          err = dabs(den_ypcp(i, j, k) - ref)
          if(err > err_Linf) err_Linf = err
          err_L2 = err_L2 + err**2
        end do
      end do
    end do
    err_L2 = dsqrt ( err_L2 / real(dm%np(3) * dm%nc(2) * dm%np(1), WP) )
    write (wrt_unit,'(A, 1I5.1, 2ES13.5)') '# interp_p2c_y', dm%np(2), err_Linf, err_L2

    ! p2c in z
    call Get_z_midp_P2C_3D(den_ppp, den_zppc, dm, dm%ibcz(:, 5))
    err_Linf = ZERO
    err_L2   = ZERO
    do k = 1, dm%nc(3)
      zc = dm%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, dm%np(2)
        yp = dm%yp(j)
        do i = 1, dm%np(1)
          xp = dm%h(1) * real(i - 1, WP)
          ref = sin_wp ( xp / FOUR ) + sin_wp(yp / FOUR) + sin_wp(zc / FOUR)
          err = dabs(den_zppc(i, j, k) - ref)
          if(err > err_Linf) err_Linf = err
          err_L2 = err_L2 + err**2
        end do
      end do
    end do
    err_L2 = dsqrt ( err_L2 / real(dm%nc(3) * dm%np(2) * dm%np(1), WP) )
    write (wrt_unit,'(A, 1I5.1, 2ES13.5)') '# interp_p2c_z', dm%np(3), err_Linf, err_L2

    close(wrt_unit)

    return 
  end subroutine
!===============================================================================
!===============================================================================
!> \brief To test this subroutine for 1st derivative.
!>
!> This subroutine is called in \ref Test_schemes. Define the logicals to choose
!> which test section is required. 
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     d             domain
!_______________________________________________________________________________
!   subroutine Test_1st_derivative(f, d)
!     use parameters_constant_mod
!     use udf_type_mod
!     use math_mod
!     implicit none
!     type(t_flow)  , intent(in) :: f
!     type(t_domain), intent(in) :: d
!     real(WP), allocatable :: fi(:), fo(:)
!     real(WP) :: xc, yc, zc
!     real(WP) :: xp, yp, zp
!     real(WP) :: ref
!     integer :: i, j, k
!     real(WP) :: err(3), errmax
!     logical :: dbg = .false.

!     logical :: dudx_P2C = .true.
!     logical :: dudx_P2P = .true.
!     logical :: dvdx_C2P = .true.
!     logical :: dvdx_C2C = .true.

!     logical :: dudy_C2P = .true.
!     logical :: dudy_C2C = .true.
!     logical :: dvdy_P2P = .true.
!     logical :: dvdy_P2C = .true.

!     if(dudx_P2C) then
!       ! du / dx, P2C
!       ! (i', j, k) --> (i, j, k)
!       allocate ( fi( dm%np(1) ) ); fi = ZERO
!       allocate ( fo( dm%nc(1) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# du/dx : P2C'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do j = 1, dm%nc(2)
!           yc = dm%yc(j)
!           fi(:) = f%qx(:, j, k)
!           call Get_1st_derivative_1D('x', 'P2C', d, fi(:), fo(:))
!           do i = 4, dm%nc(1)-3
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = cos_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = 1, 3
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = cos_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = dm%nc(1)-2, dm%nc(1)
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = cos_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if

!     if(dudx_P2P) then
!     ! du / dx, P2P
!     ! (i', j, k) --> (i', j, k)
!       allocate ( fi( dm%np(1) ) ); fi = ZERO
!       allocate ( fo( dm%np(1) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# du/dx : P2P'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do j = 1, dm%nc(2)
!           yc = dm%yc(j)
!           fi(:) = f%qx(:, j, k)
!           call Get_1st_derivative_1D('x', 'P2P', d, fi(:), fo(:))
!           do i = 4, dm%np(1)-3
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = cos_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = 1, 3
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = cos_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = dm%np(1)-2, dm%np(1)
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = cos_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if

!     if(dudy_C2P) then
!       ! du / dy, C2P
!       ! (i', j, k) --> (i', j', k)
!       allocate ( fi( dm%nc(2) ) ); fi = ZERO
!       allocate ( fo( dm%np(2) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# du/dy : C2P'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do i = 1, dm%np(1)
!           xp = dm%h(1) * real(i - 1, WP)
!           fi(:) = f%qx(i, :, k)
!           call Get_1st_derivative_1D('y', 'C2P', d, fi(:), fo(:))
!           do j = 4, dm%np(2)-3
!             yp = dm%yp(j)
!             ref = cos_wp(yp)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = 1, 3
!             yp = dm%yp(j)
!             ref = cos_wp(yp)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = dm%np(2)-2, dm%np(2)
!             yp = dm%yp(j)
!             ref = cos_wp(yp)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)   
!     end if

!     if(dudy_C2C) then
!       ! du / dy, C2C
!       ! (i', j, k) --> (i, j, k)
!       allocate ( fi( dm%nc(2) ) ); fi = ZERO
!       allocate ( fo( dm%nc(2) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# du/dy : C2C'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do i = 1, dm%np(1)
!           xp = dm%h(1) * real(i - 1, WP)
!           fi(:) = f%qx(i, :, k)
!           call Get_1st_derivative_1D('y', 'C2C', d, fi(:), fo(:))
!           do j = 4, dm%nc(2)-3
!             yc = dm%yc(j)
!             ref = cos_wp(yc)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = 1, 3
!             yc = dm%yc(j)
!             ref = cos_wp(yc)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = dm%nc(2)-2, dm%nc(2)
!             yc = dm%yc(j)
!             ref = cos_wp(yc)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if

!     if(dvdy_P2C) then
!       ! dv / dy, P2C
!       ! (i, j', k) --> (i, j, k)
!       allocate ( fi( dm%np(2) ) ); fi = ZERO
!       allocate ( fo( dm%nc(2) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# dv/dy : P2C'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do i = 1, dm%nc(1)
!           xc = dm%h(1) * (real(i - 1, WP) + HALF)
!           fi(:) = f%qy(i, :, k)
!           call Get_1st_derivative_1D('y', 'P2C', d, fi(:), fo(:))
!           do j = 4, dm%nc(2)-3
!             yc = dm%yc(j)
!             ref = cos_wp ( yc )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = 1, 3
!             yc = dm%yc(j)
!             ref = cos_wp ( yc )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = dm%nc(2)-2, dm%nc(2)
!             yc = dm%yc(j)
!             ref = cos_wp ( yc )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3) 
!     end if 
    
!     if(dvdy_P2P) then
!       ! dv / dy, P2P
!       ! (i, j', k) --> (i, j', k)
!       allocate ( fi( dm%np(2) ) ); fi = ZERO
!       allocate ( fo( dm%np(2) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# dv/dy : P2P'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do i = 1, dm%nc(1)
!           xc = dm%h(1) * (real(i - 1, WP) + HALF)
!           fi(:) = f%qy(i, :, k)
!           call Get_1st_derivative_1D('y', 'P2P', d, fi(:), fo(:))
!           do j = 4, dm%np(2)-3
!             yp = dm%yp(j)
!             ref = cos_wp ( yp )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = 1, 3
!             yp = dm%yp(j)
!             ref = cos_wp ( yp )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = dm%np(2)-2, dm%np(2)
!             yp = dm%yp(j)
!             ref = cos_wp ( yp )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3) 
!     end if

!     if(dvdx_C2C) then
!       ! du / dx, P2C
!       ! (i', j, k) --> (i, j, k)
!       allocate ( fi( dm%nc(1) ) ); fi = ZERO
!       allocate ( fo( dm%nc(1) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# du/dx : P2C'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do j = 1, dm%nc(2)
!           yc = dm%yc(j)
!           fi(:) = f%qy(:, j, k)
!           call Get_1st_derivative_1D('x', 'C2C', d, fi(:), fo(:))
!           do i = 4, dm%nc(1)-3
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = cos_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = 1, 3
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = cos_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = dm%nc(1)-2, dm%nc(1)
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = cos_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if

!     if(dvdx_C2P) then
!     ! du / dx, P2P
!     ! (i', j, k) --> (i', j, k)
!       allocate ( fi( dm%np(1) ) ); fi = ZERO
!       allocate ( fo( dm%np(1) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# du/dx : P2P'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do j = 1, dm%nc(2)
!           yc = dm%yc(j)
!           fi(:) = f%qx(:, j, k)
!           call Get_1st_derivative_1D('x', 'P2P', d, fi(:), fo(:))
!           do i = 4, dm%np(1)-3
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = cos_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = 1, 3
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = cos_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = dm%np(1)-2, dm%np(1)
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = cos_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if
    
    
!     return 
!   end subroutine


!   !===============================================================================
! !===============================================================================
! !> \brief To test this subroutine for 1st derivative.
! !>
! !> This subroutine is called in \ref Test_schemes. Define the logicals to choose
! !> which test section is required. 
! !>
! !-------------------------------------------------------------------------------
! ! Arguments
! !______________________________________________________________________________.
! !  mode           name          role                                           !
! !______________________________________________________________________________!
! !> \param[in]     d             domain
! !_______________________________________________________________________________
!   subroutine Test_2nd_derivative(f, d)
!     use parameters_constant_mod
!     use udf_type_mod
!     use math_mod
!     implicit none
!     type(t_flow),   intent(in) :: f
!     type(t_domain), intent(in) :: d
!     real(WP), allocatable :: fi(:), fo(:)
!     real(WP) :: xc, yc, zc
!     real(WP) :: xp, yp, zp
!     real(WP) :: ref
!     integer :: i, j, k
!     real(WP) :: err(3), errmax
!     logical :: dbg = .false.

!     logical :: d2udx2_P2P = .true.
!     logical :: d2vdy2_P2P = .true.
    
!     logical :: d2udy2_C2C = .true.
!     logical :: d2vdx2_C2C = .true.

    
!     if(d2udx2_P2P) then
!     ! d2u / dx2, P2P
!     ! (i', j, k) --> (i', j, k)
!       allocate ( fi( dm%np(1) ) ); fi = ZERO
!       allocate ( fo( dm%np(1) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# d2u/dx2 : P2P'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do j = 1, dm%nc(2)
!           yc = dm%yc(j)
!           fi(:) = f%qx(:, j, k)
!           call Get_2nd_derivative_1D('x', 'P2P', d, fi(:), fo(:))
!           do i = 4, dm%np(1)-3
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = -sin_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = 1, 3
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = -sin_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = dm%np(1)-2, dm%np(1)
!             xp = dm%h(1) * real(i - 1, WP)
!             ref = -sin_wp ( xp )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if

!     if(d2udy2_C2C) then
!       ! C2C
!       ! (i', j, k) --> (i, j, k)
!       allocate ( fi( dm%nc(2) ) ); fi = ZERO
!       allocate ( fo( dm%nc(2) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '# d2u/dy2 : C2C'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do i = 1, dm%np(1)
!           xp = dm%h(1) * real(i - 1, WP)
!           fi(:) = f%qx(i, :, k)
!           call Get_2nd_derivative_1D('y', 'C2C', d, fi(:), fo(:))
!           do j = 4, dm%nc(2)-3
!             yc = dm%yc(j)
!             ref = -sin_wp(yc)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = 1, 3
!             yc = dm%yc(j)
!             ref = -sin_wp(yc)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = dm%nc(2)-2, dm%nc(2)
!             yc = dm%yc(j)
!             ref = -sin_wp(yc)
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if
    
!     if(d2vdy2_P2P) then
!       ! dv / dy, P2P
!       ! (i, j', k) --> (i, j', k)
!       allocate ( fi( dm%np(2) ) ); fi = ZERO
!       allocate ( fo( dm%np(2) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '#  d2v/dy2 : P2P'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do i = 1, dm%nc(1)
!           xc = dm%h(1) * (real(i - 1, WP) + HALF)
!           fi(:) = f%qy(i, :, k)
!           call Get_2nd_derivative_1D('y', 'P2P', d, fi(:), fo(:))
!           do j = 4, dm%np(2)-3
!             yp = dm%yp(j)
!             ref = -sin_wp ( yp )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = 1, 3
!             yp = dm%yp(j)
!             ref = -sin_wp ( yp )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!           do j = dm%np(2)-2, dm%np(2)
!             yp = dm%yp(j)
!             ref = -sin_wp ( yp )
!             errmax = dabs(fo(j)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3) 
!     end if

!     if(d2vdx2_C2C) then
!       ! du / dx, P2C
!       ! (i', j, k) --> (i, j, k)
!       allocate ( fi( dm%nc(1) ) ); fi = ZERO
!       allocate ( fo( dm%nc(1) ) ); fo = ZERO
!       xc = ZERO; yc = ZERO; zc = ZERO
!       xp = ZERO; yp = ZERO; zp = ZERO
!       err = ZERO
!       write (wrt_unit,'(A)') '  '
!       write (wrt_unit,'(A)') '#  d2v/dx2 : C2C'
!       do k = 1, dm%nc(3)
!         zc = dm%h(3) * (real(k - 1, WP) + HALF)
!         do j = 1, dm%nc(2)
!           yc = dm%yc(j)
!           fi(:) = f%qy(:, j, k)
!           call Get_2nd_derivative_1D('x', 'C2C', d, fi(:), fo(:))
!           do i = 4, dm%nc(1)-3
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = -sin_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(2)) err(2) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = 1, 3
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = -sin_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(1)) err(1) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!           do i = dm%nc(1)-2, dm%nc(1)
!             xc = dm%h(1) * (real(i - 1, WP) + HALF)
!             ref = -sin_wp ( xc )
!             errmax = dabs(fo(i)-ref)
!             if (errmax > err(3)) err(3) = errmax
!             if(dbg .and. errmax>0.1_WP) & 
!             write (wrt_unit,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
!           end do
!         end do
!       end do
!       deallocate (fi)
!       deallocate (fo)
!       write (wrt_unit, '(3ES15.7)') err(1:3)
!     end if

!     return 
!   end subroutine
end module