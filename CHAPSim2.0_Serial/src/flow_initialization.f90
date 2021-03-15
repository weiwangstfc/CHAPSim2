module flow_variables_mod
  use precision_mod
  implicit none

  real(WP), save, allocatable, dimension(:, :, :) :: qx, qy, qz
  real(WP), save, allocatable, dimension(:, :, :) :: gx, gy, gz
  real(WP), save, allocatable, dimension(:, :, :) :: pres
  real(WP), save, allocatable, dimension(:, :, :) :: pcor

  real(WP), save, allocatable, dimension(:, :, :) :: dDens
  real(WP), save, allocatable, dimension(:, :, :) :: mVisc

  real(WP), save, allocatable, dimension(:, :, :) :: dh
  real(WP), save, allocatable, dimension(:, :, :) :: hEnth
  real(WP), save, allocatable, dimension(:, :, :) :: kCond
  real(WP), save, allocatable, dimension(:, :, :) :: tTemp
  
  private :: Allocate_variables
  private :: Generate_poiseuille_flow_profile
  private :: Initialize_poiseuille_flow
  private :: Initialize_vortexgreen_flow
  private :: Initialize_thermal_variables

  public  :: Initialize_flow_variables
  
contains

  subroutine Allocate_variables
    use input_general_mod, only : ithermo
    use geometry_mod
    use parameters_constant_mod, only : ZERO, ONE
    implicit none

    allocate ( qx ( 1 : domain%np(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; qx = ZERO
    allocate ( qy ( 1 : domain%nc(1), 1 : domain%np(2), 1 : domain%nc(3) )  ) ; qy = ZERO
    allocate ( qz ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%np(3) )  ) ; qz = ZERO

    allocate ( gx ( 1 : domain%np(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; gx = ZERO
    allocate ( gy ( 1 : domain%nc(1), 1 : domain%np(2), 1 : domain%nc(3) )  ) ; gy = ZERO
    allocate ( gz ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%np(3) )  ) ; gz = ZERO

    allocate ( pres ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; pres = ZERO
    allocate ( pcor ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; pcor = ZERO

    allocate ( dDens ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; dDens = ONE
    allocate ( mVisc ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; mVisc = ONE


    if(ithermo == 1) then
      allocate ( dh    ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; dh = ZERO
      allocate ( hEnth ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; hEnth = ZERO
      allocate ( kCond ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; kCond = ONE
      allocate ( tTemp ( 1 : domain%nc(1), 1 : domain%nc(2), 1 : domain%nc(3) )  ) ; tTemp = ONE
    end if

  end subroutine Allocate_variables

  subroutine Initialize_thermal_variables (d_dummy, m_dummy, dh_dummy, h_dummy, k_dummy, t_dummy)
    use input_general_mod, only : tiRef, t0Ref
    use input_thermo_mod, only : tpIni
    implicit none
    real(WP), intent(inout) :: d_dummy(:, :, :)
    real(WP), intent(inout) :: m_dummy(:, :, :)
    real(WP), intent(inout) :: dh_dummy(:, :, :)
    real(WP), intent(inout) :: h_dummy(:, :, :)
    real(WP), intent(inout) :: k_dummy(:, :, :)
    real(WP), intent(inout) :: t_dummy(:, :, :)
  
      tpIni%t = tiRef / t0Ref
      call tpIni%Refresh_thermal_properties_from_T()

      d_dummy(:, :, :)  = tpIni%d
      m_dummy(:, :, :)  = tpIni%m

      dh_dummy(:, :, :) = tpIni%dh
      h_dummy(:, :, :)  = tpIni%h
      k_dummy(:, :, :)  = tpIni%k
      t_dummy(:, :, :)  = tpIni%t
      return
  end subroutine Initialize_thermal_variables

  subroutine Generate_poiseuille_flow_profile(u_laminar, d)
    use parameters_constant_mod, only : ZERO, ONE, ONEPFIVE, TWO, MAXP, TRUNCERR
    use input_general_mod
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), intent(out) :: u_laminar(:)
    
    real(WP) :: a, b, c, yy, ymax, ymin, umean
    integer(4) :: j


    u_laminar (:) = ZERO

    ymax = d%yp( d%np_geo(2) )
    ymin = d%yp( 1 )
    if (d%case == ICASE_CHANNEL) then
      a = (ymax - ymin) / TWO
      b = ZERO
      c = ONEPFIVE
    else if (d%case == ICASE_PIPE) then
      a = (ymax - ymin)
      b = ZERO
      c = TWO
    else if (d%case == ICASE_ANNUAL) then
      a = (ymax - ymin) / TWO
      b = (ymax + ymin) / TWO
      c = TWO
    else 
      a = MAXP
      b = ZERO
      c = ONE
    end if

    do j = 1, d%nc(2)
      yy = d%yc(j)
      u_laminar(j) = ( ONE - ( (yy - b)**2 ) / a / a ) * c
    end do

    ! scale the bulk velocity to be one
    umean = ZERO
    do j = 1, d%nc(2)
      umean = umean + u_laminar(j) * (d%yp(j + 1) - d%yp(j) )
    end do
    umean = umean / (ymax - ymin)

    u_laminar(:) = u_laminar(:) / umean

    ! check the bulk velocity is one
    umean = ZERO
    do j = 1, d%nc(2)
      umean = umean + u_laminar(j) * (d%yp(j + 1) - d%yp(j) )
    end do
    umean = umean / (ymax - ymin)
    if ( abs_wp(umean - ONE) > TRUNCERR) then
      write(*, *) umean
      call Print_error_msg("Error in poiseuille_flow_profile in Subroutine" &
            // "Generate_poiseuille_flow_profile")
    end if

    return
  end subroutine Generate_poiseuille_flow_profile

  subroutine Initialize_poiseuille_flow(ux_dummy, uy_dummy, uz_dummy, pres_dummy, d)
    use random_number_generation_mod
    use parameters_constant_mod, only : ZERO, ONE
    use input_general_mod
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), intent(out) :: ux_dummy(:, :, :)
    real(WP), intent(out) :: uy_dummy(:, :, :)
    real(WP), intent(out) :: uz_dummy(:, :, :)
    real(WP), intent(out) :: pres_dummy(:, :, :)
    
    real(WP), allocatable, dimension(:) :: u_laminar
    integer :: i, j, k
    integer :: seed
    real(WP) :: rd(3)

    ! to get the profile
    allocate ( u_laminar ( d%nc(2) ) ); u_laminar(:) = ZERO
    call Generate_poiseuille_flow_profile ( u_laminar, d )

    pres_dummy(:, :, :) =  ZERO
    seed = 0 ! real random
    do k = 1, d%nc(3)
      do j = 1, d%nc(2)
        do i = 1, d%nc(1)
          seed = seed + k + j + i ! repeatable random
          call Initialize_random_number ( seed )
          call Generate_rvec_random( -ONE, ONE, 3, rd)
          ux_dummy(i, j, k) = initNoise * rd(1) + u_laminar (j)
          uy_dummy(i, j, k) = initNoise * rd(2)
          uz_dummy(i, j, k) = initNoise * rd(3)
        end do
      end do
    end do

    uy_dummy(:, 1,       :) = d%ubc(1, 2)
    uy_dummy(:, d%np(2), :) = d%ubc(2, 2)
    
    deallocate (u_laminar)
    return
  end subroutine  Initialize_poiseuille_flow

  subroutine  Initialize_vortexgreen_flow(ux_dummy, uy_dummy, uz_dummy, pres_dummy, d)
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), intent(out) :: ux_dummy(:, :, :)
    real(WP), intent(out) :: uy_dummy(:, :, :)
    real(WP), intent(out) :: uz_dummy(:, :, :)
    real(WP), intent(out) :: pres_dummy(:, :, :)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k

    do k = 1, d%nc(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          ux_dummy(i, j, k) =  sin_wp ( xp ) * cos_wp ( yc ) * cos_wp ( zc )
        end do
      end do
    end do

    do k = 1, d%nc(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%np(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%nc(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          uy_dummy(i, j, k) = -cos_wp ( xc ) * sin_wp ( yp ) * cos_wp ( zc )
        end do
      end do
    end do

    do k = 1, d%np(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%nc(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          uz_dummy(i, j, k) =  ZERO
        end do
      end do
    end do

    do k = 1, d%nc(3)
      zp = d%h(3) * real(k - 1, WP)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yp = d%yp(j)
        yc = d%yc(j)
        do i = 1, d%nc(1)
          xp = d%h(1) * real(i - 1, WP)
          xc = d%h(1) * (real(i - 1, WP) + HALF)
          pres_dummy(i, j, k)= ( cos_wp( TWO * xc       ) + cos_wp( TWO * yc       ) ) * &
                               ( cos_wp( TWO * zc + TWO ) ) / SIXTEEN
        end do
      end do
    end do
    
    return
  end subroutine Initialize_vortexgreen_flow

  subroutine  Initialize_sinetest_flow(ux_dummy, uy_dummy, uz_dummy, pres_dummy, d)
    use parameters_constant_mod, only : HALF, ZERO, SIXTEEN, TWO
    use udf_type_mod
    use math_mod
    implicit none
    type(t_domain), intent(in) :: d
    real(WP), intent(out) :: ux_dummy(:, :, :)
    real(WP), intent(out) :: uy_dummy(:, :, :)
    real(WP), intent(out) :: uz_dummy(:, :, :)
    real(WP), intent(out) :: pres_dummy(:, :, :)
    real(WP) :: xc, yc, zc
    real(WP) :: xp, yp, zp
    integer(4) :: i, j, k

    do k = 1, d%nc(3)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do j = 1, d%nc(2)
        yc = d%yc(j)
        do i = 1, d%np(1)
          xp = d%h(1) * real(i - 1, WP)
          ux_dummy(i, j, k) =  sin_wp ( xp ) + sin_wp(yc) + sin_wp(zc)
        end do 
      end do
    end do

    do k = 1, d%nc(3)
      zc = d%h(3) * (real(k - 1, WP) + HALF)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do j = 1, d%np(2)
          yp = d%yp(j)
          uy_dummy(i, j, k) = sin_wp ( xc ) + sin_wp(yp) + sin_wp(zc)
        end do
      end do
    end do

    
    do j = 1, d%nc(2)
      yc = d%yc(j)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do k = 1, d%np(3)
          zp = d%h(3) * real(k - 1, WP)
          uz_dummy(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zp)
        end do
      end do
    end do

    do j = 1, d%nc(2)
      yc = d%yc(j)
      do i = 1, d%nc(1)
        xc = d%h(1) * (real(i - 1, WP) + HALF)
        do k = 1, d%nc(3)
          zc = d%h(3) * (real(k - 1, WP) + HALF)
          pres_dummy(i, j, k) = sin_wp ( xc ) + sin_wp(yc) + sin_wp(zc)
        end do
      end do
    end do
    
    return
  end subroutine Initialize_sinetest_flow

  subroutine Initialize_flow_variables( )
    use geometry_mod
    use input_general_mod
    use parameters_constant_mod
    implicit none

    interface 
       subroutine Display_vtk_slice(d, str, varnm, vartp, var0)
        use udf_type_mod
        type(t_domain), intent( in ) :: d
        integer(4) :: vartp
        character( len = *), intent( in ) :: str
        character( len = *), intent( in ) :: varnm
        real(WP), intent( in ) :: var0(:, :, :)
       end subroutine Display_vtk_slice
    end interface

    ! allocate variables
    call Allocate_variables

    ! to initialize thermal variables 
    if (ithermo == 1) then
      call Initialize_thermal_variables (dDens, mVisc, dh, hEnth, kCond, tTemp)
    else
      dDens(:, :, :) = ONE
      mVisc(:, :, :) = ONE
    end if

    ! to initialize flow velocity and 
    if ( (icase == ICASE_CHANNEL) .or. &
         (icase == ICASE_PIPE) .or. &
         (icase == ICASE_ANNUAL) ) then

      call Initialize_poiseuille_flow (qx, qy, qz, pres, domain)

    else if (icase == ICASE_TGV) then
      
      call Initialize_vortexgreen_flow (qx, qy, qz, pres, domain)
    else if (icase == ICASE_SINETEST) then
      call Initialize_sinetest_flow (qx, qy, qz, pres, domain)
    else 
      call Print_error_msg("No such case defined in Subroutine: "//"Initialize_flow_variables" )
    end if
    ! to initialize pressure correction term
    pcor(:, :, :) = ZERO

    !call Display_vtk_slice(domain, 'xy', 'u', 1, qx)
    !call Display_vtk_slice(domain, 'xy', 'v', 2, qy)
    call Display_vtk_slice(domain, 'xy', 'p', 0, pres)

    !call Display_vtk_slice(domain, 'yz', 'v', 2, qy)
    !call Display_vtk_slice(domain, 'yz', 'w', 3, qz)
    call Display_vtk_slice(domain, 'yz', 'p', 0, pres)

    !call Display_vtk_slice(domain, 'zx', 'u', 1, qx)
    !call Display_vtk_slice(domain, 'zx', 'w', 3, qz)
    call Display_vtk_slice(domain, 'zx', 'p', 0, pres)

    ! to update mass flux terms 
    !call Refresh_massflux (flow_dummy, thermo_dummy, d)

    return
  end subroutine

end module flow_variables_mod