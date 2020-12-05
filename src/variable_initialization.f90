module variable_initialization_mod
  implicit none
  private
  public :: Initialize_variables

contains
  subroutine Initialize_variables(flow_dummy, thermo_dummy, domain_dummy, node_dummy, cell_dummy)
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
      thermo_dummy(:, :, :) = tpIni
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