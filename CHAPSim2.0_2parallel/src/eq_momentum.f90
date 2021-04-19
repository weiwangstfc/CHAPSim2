module eq_momentum_mod

  private :: Calculate_momentum_driven_source
  private :: Calculate_momentum_fractional_step
  private :: Compute_momentum_rhs

contains
  subroutine Calculate_momentum_driven_source(rhs, d, isub)
    use input_general_mod, only: idriven, drvf, &
          IDRVF_NO, IDRVF_MASSFLUX, IDRVF_SKINFRIC, IDRVF_PRESLOSS, &
          tAlpha, dt
    use solver_tools_mod, only: Calculate_y_bulk
    implicit none
    real(WP),   intent(inout) :: rhs(:, :, :)
    integer(4), intent(in   ) :: isub

    real(WP) :: rhs_bulk

    rhs_bulk = ZERO
    if(idriven == IDRVF_MASSFLUX) then

      call Calculate_y_bulk(rhs, d, rhs_bulk)
      
    else if (idriven == IDRVF_SKINFRIC) then
      rhs_bulk = - HALF * drvf * tAlpha(isub) * dt
    else if (idriven == IDRVF_PRESLOSS ) then
      ! to check this part
      rhs_bulk = - HALF * drvf * tAlpha(isub) * dt
    else 
      return
    end if

    rhs(:, :, :) = rhs(:, :, :) - rhs_bulk

    return
  end subroutine 

  subroutine Calculate_momentum_fractional_step(rhs0, rhs1, rhs1_semi, isub)
    use input_general_mod, only: tGamma, tZeta, tAlpha, dt
    implicit none
    real(WP), dimension(:, :, :), intent(in) :: rhs0, rhs1
    integer(4), intent(in) :: isub
    integer(4) :: n(3)

    n(1:3) = shape(rhs1)

    allocate( rhs_dummy (n(1), n(2), n(3)) )

    ! add explicit terms
    rhs_dummy(:, :, :) = rhs1(:, :, :)
    rhs1(:, :, :) = tGamma(isub) * rhs1(:, :, :) + &
                    tZeta (isub) * rhs0(:, :, :)
    rhs0(:, :, :) = rhs_dummy(:, :, :)

    ! add semi-implicit
    rhs1(:, :, :) = rhs1_expl(:, :, :) + &
                    tAlpha(isub) * rhs1_semi(:, :, :)

    ! times the time step 
    rhs1(:, :, :) = dt * rhs1(:, :, :)

    deallocate (rhs_dummy)

    return
  end subroutine

  subroutine Compute_momentum_rhs(f, d, isub)
    use input_general_mod, only: ithermo, iviscous, igravity, idriven, IDRVF_NO
    use udf_type_mod, only: t_flow, t_domain
    implicit none
    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer(4),     intent(in   ) :: isub
    
    ! ithermo == 0, vars
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: qxxc
    real(WP), dimension( d%np(1), d%np(2), d%nc(3) ) :: qyxp
    real(WP), dimension( d%np(1), d%nc(2), d%np(3) ) :: qzxp

    real(WP), dimension( d%np(1), d%np(2), d%nc(3) ) :: qxyp
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: qyyc
    real(WP), dimension( d%nc(1), d%np(2), d%np(3) ) :: qzyp

    real(WP), dimension( d%np(1), d%nc(2), d%np(3) ) :: qxzp
    real(WP), dimension( d%nc(1), d%np(2), d%np(3) ) :: qyzp
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: qzzc

    ! ithermo == 1, vars
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: gxxc
    real(WP), dimension( d%np(1), d%np(2), d%nc(3) ) :: gyxp
    real(WP), dimension( d%np(1), d%nc(2), d%np(3) ) :: gzxp

    real(WP), dimension( d%np(1), d%np(2), d%nc(3) ) :: gxyp
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: gyyc
    real(WP), dimension( d%nc(1), d%np(2), d%np(3) ) :: gzyp

    real(WP), dimension( d%np(1), d%nc(2), d%np(3) ) :: gxzp
    real(WP), dimension( d%nc(1), d%np(2), d%np(3) ) :: gyzp
    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: gzzc

    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: mxp
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: dmdy_xp
    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: dmdz_xp

    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: myp
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: dmdx_yp
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: dmdz_yp

    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: mzp
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: dmdx_zp
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: dmdy_zp

    real(WP), dimension( d%np(1), d%nc(2), d%nc(3) ) :: m1_rhs_semimplt
    real(WP), dimension( d%nc(1), d%np(2), d%nc(3) ) :: m2_rhs_semimplt
    real(WP), dimension( d%nc(1), d%nc(2), d%np(3) ) :: m3_rhs_semimplt

    real(WP), dimension( d%nc(1), d%nc(2), d%nc(3) ) :: div
    

    real(WP), dimension( d%np(1) ) :: foxm, foxd
    real(WP), dimension( d%np(2) ) :: foym, foyd
    real(WP), dimension( d%np(3) ) :: fozm, fozd
    
    ! common vars
    real(WP), dimension( d%np(1) ) :: fix, fox
    real(WP), dimension( d%np(2) ) :: fiy, foy
    real(WP), dimension( d%np(3) ) :: fiz, foz

    integer(4), parameter :: II = 1, JJ = 2, KK = 3
    real(WP) :: one_third_rre, two_rre

!===============================================================================
! Initilisation
!===============================================================================
    fix (:) = ZERO
    fox (:) = ZERO
    fiy (:) = ZERO
    foy (:) = ZERO
    fiz (:) = ZERO
    foz (:) = ZERO

    f%m1_rhs(:, :, :) = ZERO
    f%m2_rhs(:, :, :) = ZERO
    f%m3_rhs(:, :, :) = ZERO

    m1_rhs_semimplt(:, :, :) =  ZERO
    m2_rhs_semimplt(:, :, :) =  ZERO
    m3_rhs_semimplt(:, :, :) =  ZERO

    if(ithermo == 1) then
      foxm(:) = ZERO
      foxd(:) = ZERO
      foym(:) = ZERO
      foyd(:) = ZERO
      fozd(:) = ZERO
      fozm(:) = ZERO
      div(:, :, :) = ZERO
      one_third_rre = ONE_THIRD * f%rre
      two_rre       = TWO * f%rre
    end if
      
!===============================================================================
! interpolation
!===============================================================================
!_______________________________________________________________________________
! interpolation operation in x direction
!_______________________________________________________________________________
    do k = 1, d%np(KK)
      do j = 1, d%np(JJ)

        if( k <= d%nc(KK) .and. j <= d%nc(JJ) ) then
          ! qx, x average
          fix(:) = f%qx(:, j, k) !qx(i', j, k)
          call Get_midp_interpolation( 'x', 'P2C', d, fix(:), fox( 1 : d%nc(II) ) ) ! qx(i, j, k)
          qxxc(:, j, k) = fox( 1 : d%nc(II) )
        end if

        if( k <= d%nc(KK) ) then
          ! qy, x average
          fix( 1 : d%nc(II) ) = f%qy(:, j, k) !qy(i, j', k)
          call Get_midp_interpolation('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:) ) ! qy(i', j', k)
          qyxp(:, j, k) = fox(:)
        end if
        
        if( j <= d%nc(JJ) ) then
          ! qz, x average
          fix( 1 : d%nc(II) ) = f%qz(:, j, k)  !qz(i, j, k')
          call Get_midp_interpolation('x', 'C2P', d, fix( 1 : d%nc(II) ), fo(:) ) ! qz(i', j, k')
          qzxp(:, j, k) = fox(:)
        end if

      end do
    end do

    if(ithermo == 1) then
      do k = 1, d%np(KK)
        do j = 1, d%np(JJ)
          
          if( k <= d%nc(KK) .and. j <= d%nc(JJ) ) then
            ! gx, x average
            fix(:) = f%gx(:, j, k) !gx(i', j, k)
            call Get_midp_interpolation('x', 'P2C', d, fix(:), fox( 1 : d%nc(II) ) ) ! gx(i, j, k)
            gxxc(:, j, k) = fox( 1 : d%nc(II) )

            ! mu, x average
            fix( 1 : d%nc(II) ) = f%mVisc(:, j, k) !mu(i, j, k)
            call Get_midp_interpolation('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:) ) ! mu(i', j, k)
            mxp(:, j, k) = fox(:)
          end if

          if( k <= d%nc(KK) ) then
            ! gy, x average
            fix( 1 : d%nc(II) ) = f%gy(:, j, k)  ! gy(i, j', k)
            call Get_midp_interpolation('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:) ) ! gy(i', j', k)
            gyxp(:, j, k) = fox(:)
          end if

          if( j <= d%nc(JJ) ) then
            ! gz, x average
            fix( 1 : d%nc(II) ) = f%gz(:, j, k)  ! gz(i, j, k')
            call Get_midp_interpolation('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:) ) ! gz(i', j, k')
            gzxp(:, j, k) = fox(:)
          end if

        end do
      end do
    end if

!_______________________________________________________________________________
! interpolation  operation in y direction
!_______________________________________________________________________________    
    do k = 1, d%np(KK)
      do i = 1, d%np(II)

        if ( k<=d%nc(KK) ) then
          ! qx, y average
          fiy( 1 : d%nc(JJ) ) = f%qx(i, :, k) ! qx(i', j, k)
          call Get_midp_interpolation('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:) ) ! qx(i', j', k)
          qxyp(i, :, k) = foy(:)
        end if
        
        if ( k<=d%nc(KK) .and.  i<=d%nc(II) ) then
          ! qy, y average
          fiy(:) = f%qy(i, :, k)  ! qy(i, j', k)
          call Get_midp_interpolation('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ) ) ! qy(i, j, k)
          qyyc(i, :, k) = foy( 1 : d%nc(JJ) )
        end if

        if ( i<=d%nc(II) ) then
          !qz, y average
          fiy( 1 : d%nc(JJ) ) = f%qz(i, :, k) ! qz(i, j, k')
          call Get_midp_interpolation('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:)) ! qz(i, j', k')
          qzyp(i, :, k) = foy(:)
        end if

      end do
    end do

    if(ithermo == 1) then

      do k = 1, d%np(KK)
        do i = 1, d%np(II)

          if ( k<=d%nc(KK) ) then
            ! gx, y average
            fiy( 1 : d%nc(JJ) ) = f%gx(i, :, k)                                      ! gx(i', j, k)
            call Get_midp_interpolation('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:) ) ! gx(i', j', k)
            gxyp(i, :, k) = foy(:)
          end if
          
          if ( k<=d%nc(KK) .and.  i<=d%nc(II) ) then
            ! gy, y average
            fiy(:) = f%gy(i, :, k)                                                   ! gy(i, j', k)
            call Get_midp_interpolation('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ) ) ! gy(i, j, k)
            gyyc(i, :, k) = foy( 1 : d%nc(JJ) )
            ! mu, y average
            fiy( 1 : d%nc(JJ) ) = f%mVisc(i, :, k)                                   ! mu(i, j, k)
            call Get_midp_interpolation('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:) ) ! mu(i, j', k)
            myp(i, :, k) = foy(:)
          end if
  
          if ( i<=d%nc(II) ) then
            !gz, y average
            fiy( 1 : d%nc(JJ) ) = f%gz(i, :, k)                                     ! gz(i, j, k')
            call Get_midp_interpolation('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:)) ! gz(i, j', k')
            gzyp(i, :, k) = foy(:)
          end if
  
        end do
      end do
    end if
!_______________________________________________________________________________
! interpolation  operation in z direction
!_______________________________________________________________________________
    do j = 1, d%np(JJ)
      do i = 1, d%np(II)

        if( j <= d%nc(JJ) ) then
          ! qx, z average
          fiz( 1 : d%nc(KK) ) = f%qx(i, j, :) ! qx(i', j, k)
          call Get_midp_interpolation('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:) ) ! qx(i', j, k')
          qxzp(i, j, :) = fo(:)
        end if

        if( i <= d%nc(II) ) then
          ! qy, z average
          fiz( 1 : d%nc(KK) ) = f%qy(i, j, :) ! qy(i, j', k)
          call Get_midp_interpolation('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:) ) ! qy(i, j', k')
          qyzp(i, j, :) = foz(:)
        end if

        if( i <= d%nc(II) .and. j <= d%nc(JJ) ) then
          ! qz, z average
          fiz(:) = f%qz(i, j, :) ! qz(i, j, k')
          call Get_midp_interpolation('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) ) ! qz(i, j, k)
          qzzc(i, j, :) = foz( 1 : d%nc(KK) )
        end if

      end do
    end do

    if(ithermo == 1) then
      do j = 1, d%np(JJ)
        do i = 1, d%np(II)
  
          if( j <= d%nc(JJ) ) then
            ! gx, z average
            fiz( 1 : d%nc(KK) ) = f%gx(i, j, :) ! gx(i', j, k)
            call Get_midp_interpolation('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:) ) ! gx(i', j, k')
            gxzp(i, j, :) = foz(:)
          end if
  
          if( i <= d%nc(II) ) then
            ! gy, z average
            fiz( 1 : d%nc(KK) ) = f%gy(i, j, :) ! gy(i, j', k)
            call Get_midp_interpolation('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:) ) ! gy(i, j', k')
            gyzp(i, j, :) = foz(:)
          end if
  
          if( i <= d%nc(II) .and. j <= d%nc(JJ) ) then
            ! gz, z average
            fiz(:) = f%gz(i, j, :) ! gz(i, j, k')
            call Get_midp_interpolation('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) ) ! gz(i, j, k)
            gzzc(i, j, :) = foz( 1 : d%nc(KK) )
            ! mu, z average
            fi( 1 : d%nc(KK) ) = f%mVisc(i, j, :) ! mu(i, j, k)
            call Get_midp_interpolation('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:) ) ! mu(i, j, k')
            mzp(i, j, :) = foz(:)
          end if

        end do
      end do
    end if

!===============================================================================
! dmdx at points (preparation)
!===============================================================================
    if(ithermo == 1) then
!_______________________________________________________________________________
! dmyp/dx & dmzp/dx & d(qx)/dx,  operation in x direction
!_______________________________________________________________________________
      do k = 1, d%np(KK)
        do j = 1, d%np(JJ)
          
          if( k <= d%nc(KK)) then
            ! d(myp) / dx
            fix(1 : d%nc(II)) = myp(:, j, k)   !(i, j', k)
            call Get_1st_derivative( 'x', 'C2C', d, fix(1 : d%nc(II)), fox(1 : d%nc(II)) ) ! (i, j', k)
            dmdx_yp(:, j, k) = fox( 1 : d%nc(II) )
          end if
          
          if( j <= d%nc(JJ) ) then
            ! d(mzp) / dx
            fix(1 : d%nc(II)) = mzp(:, j, k)   !(i, j, k')
            call Get_1st_derivative( 'x', 'C2C', d, fix(1 : d%nc(II)), fox(1 : d%nc(II)) ) ! (i, j', k)
            dmdx_zp(:, j, k) = fox( 1 : d%nc(II) )
          end if

          if( k <= d%nc(KK) .and. j <= d%nc(JJ) ) then
            ! d(qx)/dx
            fix(:) = f%qx(:, j, k)                                              !qx       (i', j, k)
            call Get_1st_derivative('x', 'P2C', d, fix(:), fox(1 : d%nc(II)))   !d(qx)/dx (i, j, k)
            div(:, j, k) = div(:, j, k) + &
                           fox(1 : d%nc(II))
          end if

        end do
      end do
  !_______________________________________________________________________________
  ! dmxp/dy & dmzp/dy & d(qy)/dy, operation in y direction
  !_______________________________________________________________________________
      do k = 1, d%np(KK)
        do i = 1, d%np(II)
          
          if ( k<=d%nc(KK) ) then
            ! d(mxp) / dy
            fiy( 1 : d%nc(JJ) ) = mxp(i, :, k) ! mu(i', j, k)
            call Get_1st_derivative('y', 'C2C', d, fiy( 1 : d%nc(JJ) ), foy( 1 : d%nc(JJ) ) ) ! dm/dy (i', j, k)
            dmdy_xp(i, :, k) = foy(1 : d%nc(JJ))
          end if
          
          if ( i<=d%nc(II) ) then
            ! d(mzp) / dy
            fiy( 1 : d%nc(JJ) ) = mzp(i, :, k) ! mu(i, j, k')
            call Get_1st_derivative('y', 'C2C', d, fiy( 1 : d%nc(JJ) ), foy( 1 : d%nc(JJ) )) ! dm/dy (i, j, k')
            dmdy_zp(i, :, k) = foy(1 : d%nc(JJ))
          end if

          if ( k<=d%nc(KK) .and. i<=d%nc(II)) then
            ! d(qy) / dy
            fiy(:) = f%qy(i, :, k)                                               ! qy        (i, j', k)
            call Get_1st_derivative('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ) ) ! d(qy))/dy (i, j, k)
            div(i, :, k) = div(i, :, k) + &
                           foy(1 : d%nc(JJ))
          end if
  
        end do
      end do
  !_______________________________________________________________________________
  ! dmxp/dz & dmyp/dz operation in z direction
  !_______________________________________________________________________________
      do j = 1, d%np(JJ)
        do i = 1, d%np(II)
          
          if( j <= d%nc(JJ) ) then
            ! d(mxp) / dz
            fiz( 1 : d%nc(KK) ) = mxp(i, j, :)                                                ! mu    (i', j, k)
            call Get_1st_derivative('z', 'C2C', d, fiz( 1 : d%nc(KK) ), foz( 1 : d%nc(KK) ) ) ! dm/dz (i', j, k)
            dmdz_xp(i, j, :) = foz( 1 : d%nc(KK) )
          end if
          
          if( j <= d%nc(JJ) ) then
            ! d(myp) / dz
            fiz( 1 : d%nc(KK) ) = myp(i, j, :)                                                ! mu    (i, j', k)
            call Get_1st_derivative('z', 'C2C', d, fiz( 1 : d%nc(KK) ), foz( 1 : d%nc(KK) ) ) ! dm/dz (i, j', k)
            dmdz_yp(i, j, :) = foz( 1 : d%nc(KK) )
          end if

          if( j <= d%nc(JJ) .and. i <= d%nc(II) ) then
            ! d(qz) / dz
            fiz(:) = f%qz(i, j, :)                                               ! qz       (i, j, k')
            call Get_1st_derivative('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) ) ! d(qz)/dz (i, j, k)
            div(i, j, :) = div(i, j, :) + &
                           foz(1 : d%nc(KK))
          end if
        end do
      end do
    end if

!===============================================================================
! the RHS of momentum equation
!===============================================================================
!-------------------------------------------------------------------------------
! the RHS of momentum equation in x direction
!_______________________________________________________________________________ 
    if(ithermo == 0 ) then
      do k = 1, d%np(KK)
        do j = 1, d%np(JJ)

          if( k <= d%nc(KK) .and. j <= d%nc(JJ) ) then
            ! for x convection term (1/3): d(qx * qx)/dx at (i', j, k)
            fix( 1 : d%nc(II) ) = qxxc(:, j, k) * qxxc(:, j, k)! (qx * qx) (i, j, k)
            call Get_1st_derivative('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:)) ! d(qx*qx)/dx, (i', j, k)
            f%m1_rhs(:, j, k) = f%m1_rhs(:, j, k) + &
                                         fox(:)

            if(iviscous == IVIS_EXPLICIT) then
              ! for x diffusion term (1/1), \mu * LL(u)
              fix(:) = f%qx(:, j, k) ! qx(i', j, k)
              call Get_2nd_derivative( 'x', 'P2P', d, fix(:), fox (:) ) !LL(qx) at (i', j, k)
              f%m1_rhs(:, j, k) =  f%m1_rhs(:, j, k) + &
                                            f%rre * fox (:)
            else if (iviscous == IVIS_SEMIMPLT ) then
              ! for x diffusion term (1/1), \mu * LL(u)
              fix(:) = f%qx(:, j, k) ! qx(i', j, k)
              call Get_2nd_derivative( 'x', 'P2P', d, fix(:), fox (:) ) !LL(qx) at (i', j, k)
              m1_rhs_semimplt(:, j, k) =  m1_rhs_semimplt(:, j, k) + &
                                            f%rre * fox (:)
            end if
          end if

          if( k <= d%nc(KK)) then
            ! for y convection term (1/3), d(qx * qy)/dx at (i, j', k)
            fix(:) = qxyp(:, j, k) * qyxp(:, j, k) !(qx * qy) (i', j', k)
            call Get_1st_derivative('x', 'P2C', d, fix(:), fox( 1 : d%nc(II) )) ! d(qx*qy)/dx, (i, j', k)
            f%m2_rhs(:, j, k) = f%m2_rhs(:, j, k) + &
                                         fox( 1 : d%nc(II) )
          end if

          if( j <= d%nc(JJ)) then
            ! for z convection term (1/3), d(qx * qz)/dx at (i, j, k')
            fix(:) = qxzp(:, j, k) * qzxp(:, j, k) !d(qx * qz)(i', j, k')
            call Get_1st_derivative('x', 'P2C', d, fix(:), fox( 1 : d%nc(II) )) ! d(qx*qz)/dx, (i, j', k)
            f%m3_rhs(:, j, k) =  f%m3_rhs(:, j, k) + &
                                          fox( 1 : d%nc(II) )
          end if

        end do
      end do

    else if(ithermo == 1) then
      do k = 1, d%np(KK)
        do j = 1, d%np(JJ)

          if( k <= d%nc(KK) .and. j <= d%nc(JJ) ) then

            ! for x convection term (1/3), d(gx * qx)/dx at (i', j, k)
            fix( 1 : d%nc(II) ) = gxxc(:, j, k) * qxxc(:, j, k) !(gx * qx) (i, j, k)
            call Get_1st_derivative('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:)) ! d(gx*qx)/dx, (i', j, k)
            f%m1_rhs(:, j, k) = f%m1_rhs(:, j, k) + &
                                         fox(:)

            if(iviscous == IVIS_EXPLICIT) then
            ! for x diffusion term (1/7), \mu * LL(u)
              fix(:) = f%qx(:, j, k) ! qx(i', j, k)
              call Get_2nd_derivative( 'x', 'P2P', d, fix(:), fox (:) ) !LL(qx) at (i', j, k)
              f%m1_rhs(:, j, k) =  f%m1_rhs(:, j, k) + &
                                            f%rre * mxp(:, j, k) * fox (:)
            else if (iviscous == IVIS_SEMIMPLT) then
              ! for x diffusion term (1/7), \mu * LL(u)
              fix(:) = f%qx(:, j, k) ! qx(i', j, k)
              call Get_2nd_derivative( 'x', 'P2P', d, fix(:), fox (:) ) !LL(qx) at (i', j, k)
              m1_rhs_semimplt(:, j, k) =  m1_rhs_semimplt(:, j, k) + &
                                            f%rre * mxp(:, j, k) * fox (:)
            else
            end if

            ! for x diffusion term (2/7), \mu * 1/3 * d (div)/dx at (i', j, k)
            fix( 1 : d%nc(II) ) = div(:, j, k) ! div(i, j, k)
            call Get_1st_derivative( 'x', 'C2P', d, fix( 1 : d%nc(II) ), fox (:) ) !d (div)/dx at (i', j, k)
            f%m1_rhs(:, j, k) =  f%m1_rhs(:, j, k) + &
                                          one_third_rre * mxp(:, j, k) * fox(:)

            ! for x diffusion term (3/7), d\mu/dx * (2*du/dx -2/3 * div(u))
            call Get_midp_interpolation( 'x', 'C2P', d, fix( 1 : d%nc(II) ), foxd(:) ) !(div) at (i', j, k)
            fix(:) = f%qx(:, j, k)  ! qx(i', j, k)            
            call Get_1st_derivative('x', 'P2P', d, fix(:), fox(:))  ! du/dx, (i', j, k)
            fix(1 : d%nc(II)) = f%mVisc(:, j, k)   !mu(i, j, k)              
            call Get_1st_derivative('x', 'C2P', d, fix(:), foxm(:)) ! dmu/dx, (i', j, k)
            f%m1_rhs(:, j, k) = f%m1_rhs(:, j, k) + &
                                         two_rre * foxm(:) * ( fox(:) - ONE_THIRD * foxd(:) )

            ! for x diffusion term (4/7), d(mxp)/dy * d(qyyc)/dx at (i', j, k)
            fix( 1 : d%nc(II) ) = qyyc(:, j, k) !qy(i, j, k)
            call Get_1st_derivative( 'x', 'C2P', d, fix( 1 : d%nc(II) ), fox (:) ) !d(qy)/dx (i', j, k)
            f%m1_rhs(:, j, k) =  f%m1_rhs(:, j, k) + &
                                          f%rre * dmdy_xp(:, j, k) * fox(:)

            ! for x diffusion term (5/7), d(mxp)/dz * d(qzzc)/dx at (i', j, k)
            fix( 1 : d%nc(II) ) = qzzc(:, j, k) ! qz(i, j, k)
            call Get_1st_derivative( 'x', 'C2P', d, fix( 1 : d%nc(II) ), fox (:) ) !d(qz)/dx (i', j, k)
            f%m1_rhs(:, j, k) =  f%m1_rhs(:, j, k) + &
                                          f%rre * dmdz_xp(:, j, k) * fox(:)
          end if

          if( k <= d%nc(KK)) then
            ! for y convection term (1/3), (gx * qy)/dx at (i, j', k)
            fix(:) = gxyp(:, j, k) * qyxp(:, j, k) ! (gx * qy) (i', j', k)
            call Get_1st_derivative('x', 'P2C', d, fix(:), fox( 1 : d%nc(II) )) ! d(gx * qy)/dx (i, j', k)
            f%m2_rhs(:, j, k) = f%m2_rhs(:, j, k) + &
                                         fox( 1 : d%nc(II) )

            ! for y diffusion term (6/7), dmu/dx * dv/dx = d(myp)/dx * d(qyxp) / dx
            fix(:) = qyxp(:, j, k) ! qy(i', j', k)
            call Get_1st_derivative( 'x', 'P2C', d, fix( 1 : d%nc(II) ), fox (1 : d%nc(II)) ) !d(qy)/dx (i, j', k)
            f%m2_rhs(:, j, k) =  f%m2_rhs(:, j, k) + &
                                          f%rre * dmdx_yp(:, j, k) * fox(1 : d%nc(II))
          end if

          if( j <= d%nc(JJ)) then
            ! for z convection term (1/3), dx(gx * qz) at (i, j, k')
            fix(:) = gxzp(:, j, k) * qzxp(:, j, k) !(gx * qz)(i', j, k')
            call Get_1st_derivative('x', 'P2C', d, fix(:), fox( 1 : d%nc(II) )) !d(gx * qz)/dx (i, j, k')
            f%m3_rhs(:, j, k) =  f%m3_rhs(:, j, k) + &
                                          fox( 1 : d%nc(II) )

            ! for z diffusion term (6/7), dmu/dx * dw/dx = d(mzp)/dx * d(qzxp) / dx
            fix(:) = qzxp(:, j, k) !qz(i', j, k')
            call Get_1st_derivative( 'x', 'P2C', d, fix(:), fox (1 : d%nc(II)) ) !d(qz)/dx (i, j, k')
            f%m3_rhs(:, j, k) =  f%m3_rhs(:, j, k) + &
                                          f%rre * dmdx_zp(:, j, k) * fox(1 : d%nc(II))
          end if

        end do
      end do
    else
    end if
    ! pressure gradient in x direction
    do k = 1, d%nc(KK)
      do j = 1, d%nc(JJ)
        ! for dp/dx at (i', j, k)
        fix( 1 : d%nc(II) ) = f%pres(:, j, k)                               ! p (i, j, k)
        call Get_1st_derivative('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:)) ! d(p)/dx, (i', j, k)
        m1_rhs_semimplt(:, j, k) = m1_rhs_semimplt(:, j, k) - &
                                     fox(:)
      end do
    end do

    ! gravity force in x direction
    if( ithermo == 1 .and. ( igravity == 1 .or. igravity == -1) ) then
      do k = 1, d%nc(KK)
        do j = 1, d%nc(JJ)
          ! for dp/dx at (i', j, k)
          fix( 1 : d%nc(II) ) = f%dDens(:, j, k)                                  ! d (i, j, k)
          call Get_midp_interpolation('x', 'C2P', d, fix( 1 : d%nc(II) ), fox(:)) ! d (i', j, k)
          m1_rhs_semimplt(:, j, k) = m1_rhs_semimplt(:, j, k) + &
                                       fox(:) * f%fgravity
        end do
      end do
    end if

!-------------------------------------------------------------------------------
! the RHS of momentum equation in y direction
!_______________________________________________________________________________ 
    if(ithermo == 0 ) then
      do k = 1, d%np(KK)
        do i = 1, d%np(II)
          
          if ( k <= d%nc(KK) ) then
            ! for x convection term (2/3), d(qy * qx)/dy at (i', j, k)
            fiy(:) = qyxp(i, :, k) * qxyp(i, :, k)                               ! (qy * qx)     (i', j', k)
            call Get_1st_derivative('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ) ) ! d(qy * qx)/dy (i', j, k)
            f%m1_rhs(i, :, k) = f%m1_rhs(i, :, k) + &
                                         foy( 1 : d%nc(JJ) )
          end if

          ! for z convection term (2/3), dy(qy * qz) at (i, j, k')
          if ( i <= d%nc(II)) then
            fiy(:) = qyzp(i, :, k) * qzyp(i, :, k)                               ! (qy * qz)    (i, j', k')
            call Get_1st_derivative('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ))  ! d(qy * qz)/dy(i, j, k')
            f%m3_rhs(i, :, k) = f%m3_rhs(i, :, k) + &
                                         foy( 1 : d%nc(JJ) )
          end if

          if ( k <= d%nc(KK) .and. i <= d%nc(II)) then
            ! for y convection term (2/3), d(qy * qy)/dy at (i, j', k)
            fiy( 1 : d%nc(JJ) ) = qyyc(i, :, k) * qyyc(i, :, k)                 ! (qy * qy)     (i, j, k)
            call Get_1st_derivative('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:)) ! d(qy * qy)/dy (i, j', k)
            f%m2_rhs(i, :, k) = f%m2_rhs(i, :, k) + &
                                         foy(:)

            if(iviscous == IVIS_EXPLICIT) then
            ! for y diffusion term (1/1), \mu * LL(v)
              fiy(:) = f%qy(i, :, k)                                    !qy        (i, j', k)
              call Get_2nd_derivative( 'y', 'P2P', d, fiy(:), foy (:) ) !LL(qy) at (i, j', k)
              f%m2_rhs(i, :, k) =  f%m2_rhs(i, :, k) + &
                                            f%rre * foy (:)
            else if (iviscous == IVIS_SEMIMPLT) then
            ! for y diffusion term (1/1), \mu * LL(v)
              fiy(:) = f%qy(i, :, k)                                    !qy        (i, j', k)
              call Get_2nd_derivative( 'y', 'P2P', d, fiy(:), foy (:) ) !LL(qy) at (i, j', k)
              m2_rhs_semimplt(i, :, k) =  m2_rhs_semimplt(i, :, k) + &
                                            f%rre * foy (:)
            else
            end if
          end if

        end do
      end do

    else if (ithermo == 1) then

      do k = 1, d%np(KK)
        do i = 1, d%np(II)
          
          if ( k <= d%nc(KK) ) then
            ! for x convection term(2/3), d(gy * qx)/dy at (i', j, k)
            fiy(:) = gyxp(i, :, k) * qxyp(i, :, k)                                 ! (gy * qx)     (i', j', k)
            call Get_1st_derivative('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ) )   ! d(gy * qx)/dy (i', j, k)
            f%m1_rhs(i, :, k) = f%m1_rhs(i, :, k) + &
                                         foy( 1 : d%nc(JJ) )
            ! for x diffusion term (6/7), dmu/dy * du/dy = d(mxp)/dy * d(qxyp) / dy
            fiy(:) = qxyp(i, :, k)                                                 ! qx      (i', j', k)
            call Get_1st_derivative('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ) )   ! d(qx)/dy(i', j, k)
            f%m1_rhs(i, :, k) = f%m1_rhs(i, :, k) + &
                                         f%rre * dmdy_xp(i, :, k) * foy( 1 : d%nc(JJ) )
          end if

          if ( i <= d%nc(II)) then
            ! for z convection term (2/3), d(gy * qz)/dy at (i, j, k')
            fiy(:) = gyzp(i, :, k) * qzyp(i, :, k)                                  ! (gy * qz)     (i, j', k')
            call Get_1st_derivative('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ))     ! d(gy * qz)/dy (i, j, k')
            f%m3_rhs(i, :, k) = f%m3_rhs(i, :, k) + &
                                         foy( 1 : d%nc(JJ) )
            !for z diffusion term (7/7), dmu/dy * dw/dy = d(mzp)/dy * d(qzyp) / dy
            fiy(:) = qzyp(i, :, k)                                                  ! qz       (i, j', k')
            call Get_1st_derivative('y', 'P2C', d, fiy(:), foy( 1 : d%nc(JJ) ) )    ! d(qz)/dy (i, j, k')
            f%m3_rhs(i, :, k) = f%m3_rhs(i, :, k) + &
                                         f%rre * dmdy_zp(i, :, k) * foy( 1 : d%nc(JJ) )
          end if

          if ( k <= d%nc(KK) .and. i <= d%nc(II)) then
            ! for y convection term (2/3), d(gy * qy)/dy at (i, j', k)
            fiy( 1 : d%nc(JJ) ) = gyyc(i, :, k) * qyyc(i, :, k)                     ! (gy * qy)     (i, j, k)
            call Get_1st_derivative('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:))     ! d(gy * qy)/dy (i, j',k)
            f%m2_rhs(i, :, k) = f%m2_rhs(i, :, k) + &
                                         foy(:)

            if(iviscous == IVIS_EXPLICIT) then
            ! for y diffusion term (1/7), \mu * LL(v)
              fiy(:) = f%qy(i, :, k)                                                ! qy     (i, j', k)
              call Get_2nd_derivative( 'y', 'P2P', d, fiy(:), foy (:) )             ! LL(qv) (i, j', k)
              f%m2_rhs(i, :, k) =  f%m2_rhs(i, :, k) + &
                                            f%rre * myp(i, :, k) * foy (:)
            else if (iviscous == IVIS_SEMIMPLT) then
              ! for y diffusion term (1/7), \mu * LL(v)
              fiy(:) = f%qy(i, :, k)                                                ! qy     (i, j', k)
              call Get_2nd_derivative( 'y', 'P2P', d, fiy(:), foy (:) )             ! LL(qv) (i, j', k)
              m2_rhs_semimplt(i, :, k) =  m2_rhs_semimplt(i, :, k) + &
                                            f%rre * myp(i, :, k) * foy (:)
            else
            end if
            
            ! for y diffusion term (2/7), \mu * 1/3 * d (div)/dy at (i, j', k)
            fiy( 1 : d%nc(JJ) ) = div(i, :, k)                                      ! div           (i, j, k)
            call Get_1st_derivative( 'y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:))    ! d (div)/dy at (i, j', k)
            f%m2_rhs(i, :, k) = f%m2_rhs(i, :, k) + &
                                         one_third_rre * myp(i, :, k) * foy(:)

            ! for y diffusion term (3/7), d \mu /dy * (2*dv/dy -2/3 * div(u))
            call Get_midp_interpolation( 'y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foyd(:) ) ! (div)     (i, j', k)
            fiy(:) = f%qy(i, :, k)                                                     ! qy        (i, j', k)          
            call Get_1st_derivative('y', 'P2P', d, fiy(:), foy(:))                     ! d(qy))/dy (i, j', k)
            fiy(1 : d%nc(JJ)) = f%mVisc(i, :, k)                                       ! mu        (i, j, k)         
            call Get_1st_derivative('y', 'C2P', d, fiy(:), foym(:))                    ! dmu/dy    (i, j', k)
            f%m2_rhs(i, :, k) = f%m2_rhs(i, :, k) + &
                                         two_rre * foym(:) * ( foy(:) - ONE_THIRD * foyd(:) )

            ! for y diffusion term (4/7), dmu/dx * du/dy = d(myp)/dx * d(qxxc)/dy at (i', j, k)
            fiy( 1 : d%nc(JJ) ) = qxxc(i, :, k)                                       ! qx          (i, j, k)
            call Get_1st_derivative( 'y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy (:) )    ! d(qx)/dy    (i, j', k)
            f%m2_rhs(i, :, k) =  f%m2_rhs(i, :, k) + &
                                          f%rre * dmdx_yp(i, :, k) * foy(:)
            
            ! for y diffusion term (5/7), dmu/dz * dw/dy = d(myp)/dz * d(qzzc)/dy at (i', j, k)
            fiy( 1 : d%nc(JJ) ) = qzzc(i, :, k)                                       ! qz          (i, j, k)
            call Get_1st_derivative( 'y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy (:) )    ! d(qz)/dy    (i, j', k)
            f%m2_rhs(i, :, k) =  f%m2_rhs(i, :, k) + &
                                          f%rre * dmdz_yp(i, :, k) * foy(:)
          end if

        end do
      end do
    else
    end if
    ! pressure gradient in y direction
    do k = 1, d%nc(KK)
      do i = 1, d%nc(II)
        ! for d(p)/dy at (i, j', k)
        fiy( 1 : d%nc(JJ) ) = f%pres(i, :, k)                               ! p     (i, j, k)
        call Get_1st_derivative('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:)) ! dp/dy (i, j',k)
        m2_rhs_semimplt(i, :, k) = m2_rhs_semimplt(i, :, k) - &
                                     foy(:)
      end do
    end do

    ! gravity force in y direction
    if( ithermo == 1 .and. ( igravity == 2 .or. igravity == -2) ) then
      do k = 1, d%nc(KK)
        do i = 1, d%nc(II)
          ! for d at (i, j', k)
          fiy( 1 : d%nc(JJ) ) = f%dDens(i, :, k)                              ! d (i, j, k)
          call Get_1st_derivative('y', 'C2P', d, fiy( 1 : d%nc(JJ) ), foy(:)) ! d (i, j', k)
          m2_rhs_semimplt(i, :, k) = m2_rhs_semimplt(i, :, k) + &
                                       foy(:) * f%fgravity
        end do
      end do
    end if
      
!-------------------------------------------------------------------------------
! the RHS of momentum equation in z direction
!_______________________________________________________________________________ 
    if(ithermo == 0) then

      do j = 1, d%np(JJ)
        do i = 1, d%np(II)

          ! for x convection term (3/3), d(qz * qx)/dz at (i', j, k)
          if( j <= d%nc(JJ) ) then
            fiz(:) = qzxp(i, j, :) * qxzp(i, j, :)                                    ! (qz * qx)     (i', j, k')
            call Get_1st_derivative('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) )      ! d(qz * qx)/dz (i', j, k)
            f%m1_rhs(i, j, :) = f%m1_rhs(i, j, :) + &
                                         foz( 1 : d%nc(KK) )
          end if

          ! for y convection term (3/3), d(gz * qy)/dz at (i, j', k)
          if( i <= d%nc(II) ) then
            fiz(:) = qzyp(i, j, :) * qyzp(i, j, :)                                    ! (gz * qy)     (i, j', k')
            call Get_1st_derivative('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) )      ! d(gz * qy)/dz (i, j', k)
            f%m2_rhs(i, j, :) = f%m2_rhs(i, j, :) + &
                                         foz( 1 : d%nc(KK) )
          end if
          
          if( j <= d%nc(JJ) .and. i <= d%nc(II) ) then
            ! for z convection term (3/3), dz(gz * qz) at (i, j, k')
            fiz( 1 : d%nc(KK) ) = qzzc(i, j, :) * qzzc(i, j, :)                       ! (gz * qz)     (i, j, k)
            call Get_1st_derivative('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:))       ! d(gz * qz)/dz (i, j, k')
            f%m3_rhs(i, j, :) = f%m3_rhs(i, j, :) + &
                                         foz(:)
            
            if(iviscous == IVIS_EXPLICIT) then
            ! for z diffusion term (1/1), \mu * LL(w)
              fiz(:) = f%qz(i, j, :)                                                  ! qz     (i, j, k')
              call Get_2nd_derivative( 'z', 'P2P', d, fiz(:), foz (:) )               ! LL(qz) (i, j, k')
              f%m3_rhs(i, j, :) =  f%m3_rhs(i, j, :) + &
                                            f%rre * foz (:)
            else if (iviscous == IVIS_SEMIMPLT) then
              ! for z diffusion term (1/1), \mu * LL(w)
              fiz(:) = f%qz(i, j, :)                                                  ! qz     (i, j, k')
              call Get_2nd_derivative( 'z', 'P2P', d, fiz(:), foz (:) )               ! LL(qz) (i, j, k')
              m3_rhs_semimplt(i, j, :) =  m3_rhs_semimplt(i, j, :) + &
                                            f%rre * foz (:)
            else
            end if
          end if

        end do
      end do

    else if (ithermo == 1) then

      do j = 1, d%np(JJ)
        do i = 1, d%np(II)

          if( j <= d%nc(JJ) ) then
            ! for x convection term (3/3), dz(gz * qx) at (i', j, k)
            fiz(:) = gzxp(i, j, :) * qxzp(i, j, :)                                   ! (gz * qx)     (i', j, k')
            call Get_1st_derivative('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) )     ! d(gz * qx)/dz (i', j, k)
            f%m1_rhs(i, j, :) = f%m1_rhs(i, j, :) + &
                                         foz( 1 : d%nc(KK) )
            ! for x diffusion term (7/7), dmu/dz * du/dz = d(mxp)/dz * d(qxzp) / dz
            fiz(:) = qxzp(i, j, :)                                                   ! qx       (i', j, k')
            call Get_1st_derivative('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) )     ! d(qx)/dz (i', j, k)
            f%m1_rhs(i, j, :) = f%m1_rhs(i, j, :) + &
                                         f%rre * dmdz_xp(i, j, :) * foz( 1 : d%nc(KK) )
          end if

          if( i <= d%nc(II) ) then
            ! for y convection term (3/3), dz(gz * qy) at (i, j', k)
            fiz(:) = gzyp(i, j, :) * qyzp(i, j, :)                                   ! (gz * qy)     (i, j', k')
            call Get_1st_derivative('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) )     ! d(gz * qy)/dz (i, j', k)
            f%m2_rhs(i, j, :) = f%m2_rhs(i, j, :) + &
                                         foz( 1 : d%nc(KK) )
            ! for y diffusion term (7/7), dmu/dz * dv/dz = d(myp)/dz * d(qyzp) / dz
            fiz(:) = qyzp(i, j, :)                                                   ! qy       (i, j', k')
            call Get_1st_derivative('z', 'P2C', d, fiz(:), foz( 1 : d%nc(KK) ) )     ! d(qy)/dz (i, j', k)
            f%m2_rhs(i, j, :) = f%m2_rhs(i, j, :) + &
                                         f%rre * dmdz_yp(i, j, :) * foz( 1 : d%nc(KK) )
          end if

          
          if( j <= d%nc(JJ) .and. i <= d%nc(II) ) then
            ! for z convection term (3/3), dz(gz * qz) at (i, j, k')
            fiz( 1 : d%nc(KK) ) = gzzc(i, j, :) * qzzc(i, j, :)                      ! (gz * qz)     (i, j, k)
            call Get_1st_derivative('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:))      ! d(gz * qz)/dz (i, j, k')
            f%m3_rhs(i, j, :) = f%m3_rhs(i, j, :) + &
                                         foz(:)
            
            if(iviscous == IVIS_EXPLICIT) then
            ! for z diffusion term (1/7), \mu * LL(w)
              fiz(:) = f%qz(i, j, :)                                                 ! qz             (i, j, k')
              call Get_2nd_derivative( 'z', 'P2P', d, fiz(:), foz (:) )              ! LL(qz)         (i, j, k')
              f%m3_rhs(i, j, :) =  f%m3_rhs(i, j, :) + &
                                            f%rre * mzp(i, j, :) * foz (:)
            else if (iviscous == IVIS_SEMIMPLT)
              ! for z diffusion term (1/7), \mu * LL(w)
              fiz(:) = f%qz(i, j, :)                                                 ! qz             (i, j, k')
              call Get_2nd_derivative( 'z', 'P2P', d, fiz(:), foz (:) )              ! LL(qz)         (i, j, k')
              m3_rhs_semimplt(i, j, :) =  m3_rhs_semimplt(i, j, :) + &
                                            f%rre * mzp(i, j, :) * foz (:)
            else
            end if
    
            ! for z diffusion term (2/7), \mu * 1/3 * d (div)/dz at (i, j, k')
            fiz( 1 : d%nc(KK) ) = div(i, j, :)                                       ! div            (i, j, k)             
            call Get_1st_derivative('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:) )     ! d (div)/dz     (i, j, k')
            f%m3_rhs(i, j, :) = f%m3_rhs(i, j, :) + &
                                         one_third_rre * mzp(i, j, :) * foz(:)

            ! for z diffusion term (3/7), d \mu /dz * (2*dw/dz -2/3 * div(u))
            call Get_midp_interpolation( 'z', 'C2P', d, fiz( 1 : d%nc(KK) ), fozd(:) ) !(div)         (i, j, k')
            fiz(:) = f%qz(i, j, :)                                                     ! qz           (i, j, k')      
            call Get_1st_derivative('z', 'P2P', d, fiz(:), foz(:))                     ! d(qz)/dz,    (i, j, k')
            fiz(1 : d%nc(KK)) = f%mVisc(i, j, :)                                       ! mu           (i, j, k)
            call Get_1st_derivative('z', 'C2P', d, fiz(:), fozm(:))                    ! dmu/dz,      (i, j, k')
            f%m3_rhs(i, j, :) = f%m3_rhs(i, j, :) + &
                                         two_rre * fozm(:) * ( foz(:) - ONE_THIRD * fozd(:) )

            ! for z diffusion term (4/7), dmu/dx * du/dz = d(mzp)/dx * d(qxxc)/dz at (i, j, k')
            fiz( 1 : d%nc(JJ) ) = qxxc(i, :, k)                                        !qx            (i, j, k)
            call Get_1st_derivative( 'z', 'C2P', d, fiz( 1 : d%nc(JJ) ), foz (:) )     !d(qx)/dz      (i, j, k')
            f%m3_rhs(i, j, :) = f%m3_rhs(i, j, :) + &
                                         f%rre * dmdx_zp(i, j, :) * foz(:)

            ! for z diffusion term (5/7), dmu/dy * dv/dz = d(mzp)/dy * d(qyyc)/dz at (i, j, k')
            fiz( 1 : d%nc(JJ) ) = qyyc(i, :, k)                                        !qy            (i, j, k)
            call Get_1st_derivative( 'z', 'C2P', d, fiz( 1 : d%nc(JJ) ), foz (:) )     !d(qy)/dz      (i, j, k')
            f%m3_rhs(i, j, :) = f%m3_rhs(i, j, :) + &
                                         f%rre * dmdy_zp(i, j, :) * foz(:)
          end if  
        end do
      end do

    else 
    end if

    ! pressure gradient in z direction
    do j = 1, d%nc(JJ)
      do i = 1, d%nc(II)
        ! for d(p)/dz at (i, j, k')
        fiz( 1 : d%nc(KK) ) = f%pres(i, j, :)                               ! p     (i, j, k)
        call Get_1st_derivative('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:)) ! dp/dz (i, j, k')
        m3_rhs_semimplt(i, j, :) = m3_rhs_semimplt(i, j, :) - &
                                     foz(:)
      end do
    end do

    ! gravity force in z direction
    if( ithermo == 1 .and. ( igravity == 3 .or. igravity == -3) ) then
      do j = 1, d%nc(JJ)
        do i = 1, d%nc(II)
          ! for d at (i, j, k')
          fiz( 1 : d%nc(KK) ) = f%dDens(i, j, :)                              ! d (i, j, k)
          call Get_1st_derivative('z', 'C2P', d, fiz( 1 : d%nc(KK) ), foz(:)) ! d (i, j, k')
          m3_rhs_semimplt(i, j, :) = m3_rhs_semimplt(i, j, :) + &
                                       foz(:) * f%fgravity
        end do
      end do
    end if

!-------------------------------------------------------------------------------
! to build up rhs in total
!_______________________________________________________________________________ 
    ! x-momentum
    call Calculate_momentum_fractional_step(f%m1_rhs0, f%m1_rhs1, m1_rhs_semimplt, isub)
    if(idriven /= IDRVF_NO) call Calculate_momentum_driven_source(f%m1_rhs, d, isub) 

    ! y-momentum
    call Calculate_momentum_fractional_step(f%m2_rhs0, f%m2_rhs1, m2_rhs_semimplt, isub)

    ! z-momentum
    call Calculate_momentum_fractional_step(f%m3_rhs0, f%m3_rhs1, m3_rhs_semimplt, isub)
 
    return
  end subroutine Compute_momentum_rhs


  subroutine Calculate_provisional_mvar(rhs, u)
    use precision_mod
    implicit none
    real(WP), dimension(:, :, :), intent(in) :: rhs, u

    u(:, :, :) = u(:, :, :) + rhs(:, :, :)

    return
  end subroutine 

  subroutine Solve_momentum_eq(f, d, isub)
    use input_general_mod, only: iviscous, IVIS_SEMIMPLT, IVIS_EXPLICIT, ithermo
    use udf_type_mod, only: t_flow, t_domain
    implicit none
    type(t_flow),   intent(inout) :: f
    type(t_domain), intent(in   ) :: d
    integer(4),     intent(in   ) :: isub

!-------------------------------------------------------------------------------
! to calculate the rhs of the momenturn equation
!_______________________________________________________________________________ 
    call Compute_momentum_rhs(f, d, isub)
!-------------------------------------------------------------------------------
! to calculate provisional (q) or (g)
!_______________________________________________________________________________ 
    if(iviscous == IVIS_EXPLICIT) then

      if(ithermo == 0) then 
        call Calculate_provisional_mvar(f%m1_rhs, f%qx)
        call Calculate_provisional_mvar(f%m2_rhs, f%qy)
        call Calculate_provisional_mvar(f%m3_rhs, f%qz)
      else
        call Calculate_provisional_mvar(f%m1_rhs, f%gx)
        call Calculate_provisional_mvar(f%m2_rhs, f%gy)
        call Calculate_provisional_mvar(f%m3_rhs, f%gz)
      end if

    else if(iviscous == IVIS_SEMIMPLT) then
      !in order for a high order spacial accuracy
      ! to use Alternating direction implicit method
      ! ref: Cui2013: Convergence analysis of high-order compact 
      ! alternating direction implicit schemes for the two-dimensional 
      ! time fractional equation
    else 

    end if
!-------------------------------------------------------------------------------
! to calculate the provisional divergence constrains
!_______________________________________________________________________________
     

    return
  end subroutine

end module eq_momentum_mod