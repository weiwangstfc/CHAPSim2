module test_algrithms_mod
  use operations

  public   :: Test_schemes
  private  :: Test_interpolation
  private  :: Test_1st_derivative
  private  :: Test_2nd_derivative
  
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
  use type_vars_mod, only : flow, domain
  implicit none

  !call Test_TDMA_cyclic
  !call Test_TDMA_noncyclic
  call Test_interpolation (flow, domain)
  call Test_1st_derivative(flow, domain)
  call Test_2nd_derivative(flow, domain)
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
  subroutine Test_interpolation(f, d)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: d
    type(t_flow),   intent(in) :: f
    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    integer(4) :: i, j, k
    real(WP) :: err(3), errmax
    logical :: dbg = .false.
    logical :: uix_p2c = .true.
    logical :: vix_c2p = .true.
    logical :: uiy_c2p = .true.
    logical :: viy_p2c = .true.

    if(uix_p2c) then
      !test interpolation. u in x, P2C
      !(i', j, k) --> (i, j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%nc(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# ave(u)_x : P2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = f%qx(:, j, k)
          call Get_midp_interpolation_1D('x', 'P2C', d, fi(:), fo(:))

          do i = 4, d%nc(1)-3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = 1, 3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = d%nc(1)-2, d%nc(1)
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(uiy_c2p) then
      ! test interpolation. u in y, C2P 
      ! (i', j, k) --> (i', j', k)
      allocate ( fi( d%nc(2) ) ); fi = ZERO
      allocate ( fo( d%np(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# ave(u)_y : C2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          fi(:) = f%qx(i, :, k)
          call Get_midp_interpolation_1D('y', 'C2P', d, fi(:), fo(:))
          do j = 4, d%np(2)-3
            yp = d%yp(j)
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yp = d%yp(j)
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%np(2)-2, d%np(2)
            yp = d%yp(j)
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(vix_c2p) then
      !test interpolation. v in x, C2P
      !(i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(1) ) ); fi = ZERO
      allocate ( fo( d%np(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# ave(v)_x : C2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yp = d%yp(j)
          fi(:) = f%qy(:, j, k)
          call Get_midp_interpolation_1D('x', 'C2P', d, fi(:), fo(:))
          do i = 4, d%np(1)-3
            xp = d%h(1) * (real(i - 1, WP))
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = 1, 3
            xp = d%h(1) * (real(i - 1, WP))
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

          do i = d%np(1)-2, d%np(1)
            xp = d%h(1) * (real(i - 1, WP))
            ref = sin_wp ( xp ) + sin_wp(yp) + sin_wp(zc)
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do

        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(viy_p2c) then
      ! test interpolation. v in y, P2C 
      ! (i', j, k) --> (i', j', k)
      allocate ( fi( d%np(2) ) ); fi = ZERO
      allocate ( fo( d%nc(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# ave(v)_y : P2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%nc(1)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          fi(:) = f%qy(i, :, k)
          call Get_midp_interpolation_1D('y', 'P2C', d, fi(:), fo(:))
          do j = 4, d%nc(2)-3
            yc = d%yc(j)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yc = d%yc(j)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%nc(2)-2, d%nc(2)
            yc = d%yc(j)
            ref = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) &
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

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
  subroutine Test_1st_derivative(f, d)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_flow)  , intent(in) :: f
    type(t_domain), intent(in) :: d
    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    integer(4) :: i, j, k
    real(WP) :: err(3), errmax
    logical :: dbg = .false.

    logical :: dudx_P2C = .true.
    logical :: dudx_P2P = .true.
    logical :: dvdx_C2P = .true.
    logical :: dvdx_C2C = .true.

    logical :: dudy_C2P = .true.
    logical :: dudy_C2C = .true.
    logical :: dvdy_P2P = .true.
    logical :: dvdy_P2C = .true.

    if(dudx_P2C) then
      ! du / dx, P2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%nc(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# du/dx : P2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = f%qx(:, j, k)
          call Get_1st_derivative_1D('x', 'P2C', d, fi(:), fo(:))
          do i = 4, d%nc(1)-3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%nc(1)-2, d%nc(1)
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(dudx_P2P) then
    ! du / dx, P2P
    ! (i', j, k) --> (i', j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%np(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# du/dx : P2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = f%qx(:, j, k)
          call Get_1st_derivative_1D('x', 'P2P', d, fi(:), fo(:))
          do i = 4, d%np(1)-3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%np(1)-2, d%np(1)
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(dudy_C2P) then
      ! du / dy, C2P
      ! (i', j, k) --> (i', j', k)
      allocate ( fi( d%nc(2) ) ); fi = ZERO
      allocate ( fo( d%np(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# du/dy : C2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          fi(:) = f%qx(i, :, k)
          call Get_1st_derivative_1D('y', 'C2P', d, fi(:), fo(:))
          do j = 4, d%np(2)-3
            yp = d%yp(j)
            ref = cos_wp(yp)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yp = d%yp(j)
            ref = cos_wp(yp)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%np(2)-2, d%np(2)
            yp = d%yp(j)
            ref = cos_wp(yp)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)   
    end if

    if(dudy_C2C) then
      ! du / dy, C2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(2) ) ); fi = ZERO
      allocate ( fo( d%nc(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# du/dy : C2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          fi(:) = f%qx(i, :, k)
          call Get_1st_derivative_1D('y', 'C2C', d, fi(:), fo(:))
          do j = 4, d%nc(2)-3
            yc = d%yc(j)
            ref = cos_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yc = d%yc(j)
            ref = cos_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%nc(2)-2, d%nc(2)
            yc = d%yc(j)
            ref = cos_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(dvdy_P2C) then
      ! dv / dy, P2C
      ! (i, j', k) --> (i, j, k)
      allocate ( fi( d%np(2) ) ); fi = ZERO
      allocate ( fo( d%nc(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# dv/dy : P2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%nc(1)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          fi(:) = f%qy(i, :, k)
          call Get_1st_derivative_1D('y', 'P2C', d, fi(:), fo(:))
          do j = 4, d%nc(2)-3
            yc = d%yc(j)
            ref = cos_wp ( yc )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yc = d%yc(j)
            ref = cos_wp ( yc )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%nc(2)-2, d%nc(2)
            yc = d%yc(j)
            ref = cos_wp ( yc )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3) 
    end if 
    
    if(dvdy_P2P) then
      ! dv / dy, P2P
      ! (i, j', k) --> (i, j', k)
      allocate ( fi( d%np(2) ) ); fi = ZERO
      allocate ( fo( d%np(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# dv/dy : P2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%nc(1)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          fi(:) = f%qy(i, :, k)
          call Get_1st_derivative_1D('y', 'P2P', d, fi(:), fo(:))
          do j = 4, d%np(2)-3
            yp = d%yp(j)
            ref = cos_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yp = d%yp(j)
            ref = cos_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%np(2)-2, d%np(2)
            yp = d%yp(j)
            ref = cos_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3) 
    end if

    if(dvdx_C2C) then
      ! du / dx, P2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(1) ) ); fi = ZERO
      allocate ( fo( d%nc(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# du/dx : P2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = f%qy(:, j, k)
          call Get_1st_derivative_1D('x', 'C2C', d, fi(:), fo(:))
          do i = 4, d%nc(1)-3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%nc(1)-2, d%nc(1)
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = cos_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(dvdx_C2P) then
    ! du / dx, P2P
    ! (i', j, k) --> (i', j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%np(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# du/dx : P2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = f%qx(:, j, k)
          call Get_1st_derivative_1D('x', 'P2P', d, fi(:), fo(:))
          do i = 4, d%np(1)-3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%np(1)-2, d%np(1)
            xp = d%h(1) * real(i - 1, WP)
            ref = cos_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if
    
    
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
  subroutine Test_2nd_derivative(f, d)
    use parameters_constant_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_flow),   intent(in) :: f
    type(t_domain), intent(in) :: d
    real(WP), allocatable :: fi(:), fo(:)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    real(WP) :: ref
    integer(4) :: i, j, k
    real(WP) :: err(3), errmax
    logical :: dbg = .false.

    logical :: d2udx2_P2P = .true.
    logical :: d2vdy2_P2P = .true.
    
    logical :: d2udy2_C2C = .true.
    logical :: d2vdx2_C2C = .true.

    
    if(d2udx2_P2P) then
    ! d2u / dx2, P2P
    ! (i', j, k) --> (i', j, k)
      allocate ( fi( d%np(1) ) ); fi = ZERO
      allocate ( fo( d%np(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# d2u/dx2 : P2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = f%qx(:, j, k)
          call Get_2nd_derivative_1D('x', 'P2P', d, fi(:), fo(:))
          do i = 4, d%np(1)-3
            xp = d%h(1) * real(i - 1, WP)
            ref = -sin_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xp = d%h(1) * real(i - 1, WP)
            ref = -sin_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%np(1)-2, d%np(1)
            xp = d%h(1) * real(i - 1, WP)
            ref = -sin_wp ( xp )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    if(d2udy2_C2C) then
      ! C2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(2) ) ); fi = ZERO
      allocate ( fo( d%nc(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '# d2u/dy2 : C2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          fi(:) = f%qx(i, :, k)
          call Get_2nd_derivative_1D('y', 'C2C', d, fi(:), fo(:))
          do j = 4, d%nc(2)-3
            yc = d%yc(j)
            ref = -sin_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yc = d%yc(j)
            ref = -sin_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%nc(2)-2, d%nc(2)
            yc = d%yc(j)
            ref = -sin_wp(yc)
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, i, j, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if
    
    if(d2vdy2_P2P) then
      ! dv / dy, P2P
      ! (i, j', k) --> (i, j', k)
      allocate ( fi( d%np(2) ) ); fi = ZERO
      allocate ( fo( d%np(2) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '#  d2v/dy2 : P2P'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do i = 1, d%nc(1)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          fi(:) = f%qy(i, :, k)
          call Get_2nd_derivative_1D('y', 'P2P', d, fi(:), fo(:))
          do j = 4, d%np(2)-3
            yp = d%yp(j)
            ref = -sin_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = 1, 3
            yp = d%yp(j)
            ref = -sin_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
          do j = d%np(2)-2, d%np(2)
            yp = d%yp(j)
            ref = -sin_wp ( yp )
            errmax = dabs(fo(j)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(j), ref, dabs(ref-fo(j))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3) 
    end if

    if(d2vdx2_C2C) then
      ! du / dx, P2C
      ! (i', j, k) --> (i, j, k)
      allocate ( fi( d%nc(1) ) ); fi = ZERO
      allocate ( fo( d%nc(1) ) ); fo = ZERO
      xc = ZERO; yc = ZERO; zc = ZERO
      xp = ZERO; yp = ZERO; zp = ZERO
      err = ZERO
      write(*,'(A)') '  '
      write(*,'(A)') '#  d2v/dx2 : C2C'
      do k = 1, d%nc(3)
        zc = d%h(3) * (real(k - 1, WP) + HALF)
        do j = 1, d%nc(2)
          yc = d%yc(j)
          fi(:) = f%qy(:, j, k)
          call Get_2nd_derivative_1D('x', 'C2C', d, fi(:), fo(:))
          do i = 4, d%nc(1)-3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = -sin_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(2)) err(2) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = 1, 3
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = -sin_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(1)) err(1) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
          do i = d%nc(1)-2, d%nc(1)
            xc = d%h(1) * (real(i - 1, WP) + HALF)
            ref = -sin_wp ( xc )
            errmax = dabs(fo(i)-ref)
            if (errmax > err(3)) err(3) = errmax
            if(dbg .and. errmax>0.1_WP) & 
            write(*,'(3I5, 2F8.4, 1ES15.7)') k, j, i, fo(i), ref, dabs(ref-fo(i))
          end do
        end do
      end do
      deallocate (fi)
      deallocate (fo)
      write(*, '(3ES15.7)') err(1:3)
    end if

    return 
  end subroutine
end module