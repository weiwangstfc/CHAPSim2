module continuity_eq_mod

  private :: Calculate_drhodt
  public  :: Get_divergence
contains

  subroutine Calculate_drhodt(dDens, dDensm1, dDensm2, drhodt, itime)
    use precision_mod
    use udf_type_mod, only : t_flow, t_domain
    use input_general_mod, only: dt, iTimeScheme, ITIME_RK3, ITIME_RK3_CN, ITIME_AB2, &
         nsubitr, tAlpha
    use parameters_constant_mod, only: ONEPFIVE, TWO, HALF
    implicit none

    real(WP), dimension(:, :, :), intent ( in  ) :: dDens, dDensm1, dDensm2
    real(WP), dimension(:, :, :), intent ( out ) :: drhodt
    integer(4),     intent( in  ) :: itime
    
    integer(4) :: i

    if(itime == 0 .or. itime == 1 ) then
      drhodt(:, :, :) = dDens(:, :, :) - dDensm1(:, :, :)
      drhodt(:, :, :) = drhodt(:, :, :) * dt
    else 

      if(iTimeScheme == ITIME_AB2) then

        drhodt(:, :, :) = ONEPFIVE * dDens  (:, :, :) - &
                          TWO      * dDensm1(:, :, :) + &
                          HALF     * dDensm2(:, :, :)
        drhodt(:, :, :) = drhodt(:, :, :) * dt
      else if (iTimeScheme == ITIME_RK3 .or. iTimeScheme == ITIME_RK3_CN) then
        ! to check this part!
        drhodt(:, :, :) = dDens  (:, :, :)
        do i = 1, nsubitr
          drhodt(:, :, :) = drhodt(:, :, :) + tAlpha(i) * &
                            (dDensm1(:, :, :) - dDensm2(:, :, :))
        end do
      else  
        ! default, Euler 1st order 
        drhodt(:, :, :) = dDens(:, :, :) - dDensm1(:, :, :)
        drhodt(:, :, :) = drhodt(:, :, :) * dt
      end if
    end if

    return
  end subroutine Calculate_drhodt


  subroutine Get_divergence(ux, uy, uz, div, d)
    use precision_mod
    use udf_type_mod, only: t_domain, t_flow
    use parameters_constant_mod, only: ZERO
    use operations, only: Get_1st_derivative
    implicit none

    type(t_domain), intent( in ) :: d
    real(WP), dimension(:, :, :), intent (in   ) :: ux, uy, uz
    real(WP), dimension(:, :, :), intent (inout) :: div

    real(WP), allocatable :: fi(:), fo(:)

    integer(4) :: i, j, k

!-------------------------------------------------------------------------------
! du/dx at cell centre
!_______________________________________________________________________________
    allocate ( fi( d%np(1) ) ); fi = ZERO
    allocate ( fo( d%nc(1) ) ); fo = ZERO
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        fi(:) = ux(:, j, k)
        call Get_1st_derivative('x', 'P2C', d, fi(:), fo(:))
        div(:, j, k) = div(:, j, k) + fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! dv/dy at cell centre
!_______________________________________________________________________________
    allocate ( fi( d%np(2) ) ); fi = ZERO
    allocate ( fo( d%nc(2) ) ); fo = ZERO
    do k = 1, d%nc(3)
      do i = 1, d%nc(1)
        fi(:) = uy(i, :, k)
        call Get_1st_derivative('y', 'P2C', d, fi(:), fo(:))
        div(i, :, k) = div(i, :, k) + fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)
!-------------------------------------------------------------------------------
! $dw/dz$ at cell centre
!_______________________________________________________________________________
    allocate ( fi( d%np(3) ) ); fi = ZERO
    allocate ( fo( d%nc(3) ) ); fo = ZERO
    do j = 1, d%nc(2)
      do i = 1, d%nc(1)
        fi(:) = uz(i, j, :)
        call Get_1st_derivative('z', 'P2C', d, fi(:), fo(:))
        div(i, j, :) = div(i, j, :) + fo(:)
      end do
    end do
    deallocate (fi)
    deallocate (fo)

    return
  end subroutine


  subroutine Calculate_continuity_constrains(f, d, itime)
    use precision_mod
    use udf_type_mod, only: t_domain, t_flow
    use input_general_mod, only: ithermo
    use parameters_constant_mod, only: ZERO
    implicit none

    type(t_domain), intent( in ) :: d
    type(t_flow),   intent( in ) :: f                    
    integer(4),     intent( in ) :: itime


    real(WP), allocatable :: div(:, :, :)

    allocate ( div( d%nc(1), d%nc(2), d%nc(3) ) ); div = ZERO

!-------------------------------------------------------------------------------
! $d\rho / dt$ at cell centre
!_______________________________________________________________________________
    if (ithermo == 1) then
      call Calculate_drhodt(f%dDens, f%dDensm1, f%dDensm2, div, itime)
    end if
!-------------------------------------------------------------------------------
! $d(\rho u_i)) / dx_i $ at cell centre
!_______________________________________________________________________________
    if (ithermo == 1) then
      call Get_divergence(f%gx, f%gy, f%gz, div, d)
    else
      call Get_divergence(f%qx, f%qy, f%qz, div, d)
    end if

    write(*,*) "-------------------------------------------------------------------------------"
    write(*,*) "The maximum divergence :"
    write(*,'(12X, 1ES13.5)') MAXVAL(div)
    write(*,*) "-------------------------------------------------------------------------------"
   

    deallocate (div)

    return
  end subroutine


end module continuity_eq_mod
