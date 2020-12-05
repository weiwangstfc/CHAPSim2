module variables_refresh
  implicit none
  private
  public :: Refresh_massflux

contains
  subroutine Refresh_massflux (flow_dummy, thermo_dummy, domain_dummy)
    type(thermoProperty_t), intent(in) :: thermo_dummy(:, :, :)
    type(flow_t), intent(inout) :: flow_dummy(:, :, :)
    type(domain_t), intent(in) :: domain_dummy


  end subroutine Refresh_massflux

end module 