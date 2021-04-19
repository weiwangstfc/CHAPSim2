module variable_initialization_mod
  implicit none
  private
  public :: Initialize_flow_variables
  public :: Initialize_thermo_variables

contains
  subroutine Initialize_thermo_hkt_variables(dh_dummy, h_dummy, k_dummy, t_dummy)
    use precision_mod
    use input_general_mod, only: ithermo, t0Ref, tiRef
    use input_thermo_mod
    real(WP), allocatable, dimension(:, :, :), intent(out) :: dh_dummy
    real(WP), allocatable, dimension(:, :, :), intent(out) :: h_dummy
    real(WP), allocatable, dimension(:, :, :), intent(out) :: k_dummy
    real(WP), allocatable, dimension(:, :, :), intent(out) :: t_dummy

    if (ithermo /= 1) return

    tpIni%t = tiRef / t0Ref
    call tpIni%Refresh_thermal_properties_from_T()

    dh_dummy(:, :, :) = tpIni%dh
    h_dummy(:, :, :)  = tpIni%h
    k_dummy(:, :, :)  = tpIni%k
    t_dummy(:, :, :)  = tpIni%t

  end subroutine Initialize_thermo_variables


  subroutine Initialize_thermo_dm_variables(d_dummy, m_dummy, dp_dummy, mp_dummy)
    use precision_mod
    use input_general_mod, only: ithermo, t0Ref, tiRef
    use input_thermo_mod
    real(WP), allocatable, dimension(:, :, :), intent(out) :: d_dummy
    real(WP), allocatable, dimension(:, :, :), intent(out) :: m_dummy
    real(WP), allocatable, dimension(:, :, :), intent(out) :: dp_dummy
    real(WP), allocatable, dimension(:, :, :), intent(out) :: mp_dummy

    if (ithermo == 1) then

      tpIni%t = tiRef / t0Ref
      call tpIni%Refresh_thermal_properties_from_T()
      d_dummy(:, :, :)  = tpIni%d
      m_dummy(:, :, :)  = tpIni%m

      call 

    else 

      d_dummy(:, :, :) = ONE
      m_dummy(:, :, :) = ONE
      dp_dummy(:, :, :) = ONE
      mp_dummy(:, :, :) = ONE

    end if


  end subroutine Initialize_dm_variables

  subroutine Initialize_variables(, domain_dummy, node_dummy, cell_dummy)
    use input_general_mod
    use input_thermo_mod
    type(thermoProperty_t), intent(out) :: thermo_dummy(:, :, :)
    type(flow_t), intent(out) :: flow_dummy(:, :, :)
    type(domain_t), intent(in) :: domain_dummy
    type(cell_t), intent(in) :: cell_dummy(:, :, :)
    type(node_t), intent(in) :: node_dummy(:, :, :)

    ! to initialize thermal variables 
    if (ithermo == 1) then

      tpIni%t = tpIni0 / tpRef0
      call tpIni%Refresh_thermal_properties_from_T()
      d_dummy(:, :, :)  = tpIni%d
      dh_dummy(:, :, :) = tpIni%dh
      h_dummy(:, :, :)  = tpIni%h
      k_dummy(:, :, :)  = tpIni%k
      m_dummy(:, :, :)  = tpIni%m
      t_dummy(:, :, :)  = tpIni%t

    else
      call thermo_dummy(:, :, :)%Get_initialized_thermal_properties()
    end if

    ! to initialize flow variables
    if ( (icase == ICASE_CHANNEL) .or. &
         (icase == ICASE_PIPE) .or. &
         (icase == ICASE_ANNUAL) ) then

      call Initialize_poiseuille_flow (flow_dummy, domain_dummy, node_dummy, cell_dummy)

    else if (icase == ICASE_TGV) then
      
      call vortexgreen_flow_initialization (flow_dummy, domain_dummy, node_dummy, cell_dummy)
    
    else 
      
    end if

    ! to update mass flux terms 
    call Refresh_massflux (flow_dummy, thermo_dummy, domain_dummy)


  end subroutine


end module