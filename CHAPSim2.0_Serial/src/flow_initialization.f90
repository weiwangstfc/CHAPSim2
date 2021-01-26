module flow_variables_mod
  use precision_mod
  use geometry_mod
  use input_general_mod
  implicit none

  real(WP), save, allocatable, dimension(:, :, :, :) :: u_xpencil
  real(WP), save, allocatable, dimension(:, :, :, :) :: g_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: pre_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: phi_xpencil

  real(WP), save, allocatable, dimension(:, :, :) :: d_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: m_xpencil

  real(WP), save, allocatable, dimension(:, :, :) :: dh_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: h_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: k_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: t_xpencil
  
  private
  private :: Allocate_variables_xpencil
  private :: Generate_poiseuille_flow_profile
  private :: Initialize_poiseuille_flow
  private :: Initialize_thermal_variables
  private :: Initialize_thermo_variables

  public  :: Initialize_flow_variables
  
contains

  subroutine Allocate_variables
    use input_general_mod, only : ithermo, NDIM
    use domain_decomposition_mod
    use parameters_constant_mod
    implicit none

    integer(4) :: i0, i1
    integer(4) :: j0, j1
    integer(4) :: k0, k1

    
    i0 = local_xpencil%irange(1)
    i1 = local_xpencil%irange(2)

    j0 = local_xpencil%jrange(1)
    j1 = local_xpencil%jrange(2)

    k0 = local_xpencil%krange(1)
    k1 = local_xpencil%krange(2)

    ! allocate xpencil based data
    allocate ( u_xpencil ( i0 : i1, j0 : j1, k0 : k1, NDIM )  ) ; u_xpencil = ZERO
    allocate ( g_xpencil ( i0 : i1, j0 : j1, k0 : k1, NDIM )  ) ; g_xpencil = ZERO

    allocate ( pre_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; pre_xpencil = ZERO
    allocate ( phi_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; phi_xpencil = ZERO
    
    allocate ( d_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; d_xpencil = ONE
    allocate ( m_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; m_xpencil = ONE


    if(ithermo == 1) then
      allocate ( dh_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; dh_xpencil = ZERO
      allocate ( h_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; h_xpencil = ZERO
      allocate ( k_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; k_xpencil = ONE
      allocate ( t_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; t_xpencil = ONE
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

  subroutine Generate_poiseuille_flow_profile(u_laminar)
    use input_general_mod, only: lyt, lyb, icase
    use geometry_mod, only: domain
    implicit none
    real(WP), intent(out) :: u_laminar(:)
    real(WP), intent(in) :: u_dummy(:)
    
    real(WP) :: a, b, c, yy

    u_laminar (:) = ZERO

    if (icase == ICASE_CHANNEL) then
      a = (lyt - lyb) / TWO
      b = ZERO
      c = ONEPFIVE
    else if (icase == ICASE_PIPE) then
      a = lyt - lyb
      b = ZERO
      c = TWO
    else if (icase == ICASE_ANNUAL) then
      a = (lyt - lyb) / TWO
      b = (lyt + lyb) / TWO
      c = TWO
    else 
      a = MAXP
      b = ZERO
      c = ONE
    end if

    do j = 1, domain%nc(2)
      yy = domain%yc(j)
      u_laminar(j) = ( ONE - ( (yy - b)**2 ) / a / a ) * c
    end do

    ! check the mean velocity
    umean = ZERO
    do j = 1, domain%nc(2)
      umean = umean + u_laminar(j) * (domain%yp(j + 1) - domain%yp(j) )
    end do
    umean = umean / (lyt - lyb)
    if ( abs_wp(umean - ONE) > TRUNCERR) &
    call Print_error_msg("Error in poiseuille_flow_profile.")

    return
  end subroutine Generate_poiseuille_flow_profile

  subroutine Initialize_poiseuille_flow(pencil_dummy, u_dummy, pre_dummy)
    use random_number_generation_mod
    use geometry_mod, only: domain
    use domain_decomposition_mod
    implicit none
    type(pencil_t), intent(in) :: pencil_dummy
    real(WP), intent(out) :: u_dummy(:, :, :, :)
    real(WP), intent(out) :: pre_dummy(:, :, :)
    
    real(WP), allocatable, dimension(:) :: u_laminar
    integer :: i, j, k
    integer :: seed
    real(WP) :: rd(3)

    allocate ( u_laminar (domain%nc(2) ) ); u_laminar(:) = ZERO
    call Generate_poiseuille_flow_profile ( u_laminar)

    pre_dummy(:, :, :) =  ZERO
    seed = 0 ! real random
    do k = pencil_dummy%krange(1), pencil_dummy%krange(2)
      do j = pencil_dummy%jrange(1), pencil_dummy%jrange(2)
        do i = pencil_dummy%irange(1), pencil_dummy%irange(2)
          seed = seed + k + j + i ! repeatable random
          call Initialize_random_number ( seed )
          call Generate_rvec_random ( -ONE, ONE, 3, rd)
          u_dummy(:, :, :, 1) = initNoise * rd(1) + u_laminar (j)
          u_dummy(:, :, :, 2) = initNoise * rd(2)
          u_dummy(:, :, :, 3) = initNoise * rd(3)
        end do
      end do
    end do

    deallocate (u_laminar)
    return
  end subroutine  Initialize_poiseuille_flow

  subroutine  Initialize_vortexgreen_flow(pencil_dummy, u_dummy, pre_dummy)
    use domain_decomposition_mod
    use geometry_mod
    implicit none
    type(pencil_t), intent(in) :: pencil_dummy
    real(WP), intent(out) :: u_dummy(:, :, :, :)
    real(WP), intent(out) :: pre_dummy(:, :, :)


    do k = pencil_dummy%krange(1), pencil_dummy%krange(2)
      do j = pencil_dummy%jrange(1), pencil_dummy%jrange(2)
        do i = pencil_dummy%irange(1), pencil_dummy%irange(2)
          u_dummy(i, j, k, 1) = sin_wp ( node_dummy(i, j, k)%x ) * &
                                cos_wp ( cell_dummy(i, j, k)%y ) * &
                                cos_wp ( cell_dummy(i, j, k)%z )
          u_dummy(i, j, k, 2) = -cos_wp ( cell_dummy(i, j, k)%x ) * &
                                sin_wp ( node_dummy(i, j, k)%y ) * &
                                cos_wp ( cell_dummy(i, j, k)%z )
          u_dummy(i, j, k, 3) = ZERO
          pre_dummy(i, j, k)  = ( cos_wp( TWO * cell_dummy(i, j, k)%x       ) + &
                                  cos_wp( TWO * cell_dummy(i, j, k)%y       ) ) * &
                                ( cos_wp( TWO * cell_dummy(i, j, k)%z + TWO ) ) / SIXTEEN
        end do
      end do
    end do
    return
  end subroutine Initialize_vortexgreen_flow

  subroutine Initialize_flow_variables( )
    use geometry_mod
    use input_general_mod
    implicit none
    ! allocate variables
    call Allocate_variables

    ! to initialize thermal variables 
    if (ithermo == 1) then
      call Initialize_thermal_variables (d_xpencil, m_xpencil, dh_xpencil, &
            h_xpencil, k_xpencil, t_xpencil)
    else
      d_xpencil(:) = ONE
      m_xpencil(:) = ONE
    end if

    ! to initialize flow variables
    if ( (icase == ICASE_CHANNEL) .or. &
         (icase == ICASE_PIPE) .or. &
         (icase == ICASE_ANNUAL) ) then

      call Initialize_poiseuille_flow (local_xpencil, u_xpencil, pre_xpencil)

    else if (icase == ICASE_TGV) then
      
      call vortexgreen_flow_initialization (local_xpencil, u_xpencil, pre_xpencil)
    
    else 
      call Print_error_msg("No such case defined.")
    end if

    ! to initialize 

    ! to update mass flux terms 
    !call Refresh_massflux (flow_dummy, thermo_dummy, domain_dummy)

    return
  end subroutine


  


  

  

  

end module flow_variables_mod