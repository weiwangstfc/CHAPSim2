module flow_variables_mod
  use precision_mod
  use input_thermo_mod, only : thermoProperty_t
  implicit none

  type flow_t
    real(WP) :: u(3)
    real(WP) :: g(3)
    real(WP) :: pre
    real(WP) :: phi
  end type flow_t

  type(flow_t),           save, allocatable, dimension(:, :, :) :: flow_xpencil
  type(thermoProperty_t), save, allocatable, dimension(:, :, :) :: thermo_xpencil

  private
  public :: Allocate_variables

contains

  subroutine Allocate_variables ()
    use input_general_mod, only : ithermo
    use domain_decomposition_mod
    implicit none
    
    ! allocate xpencil based data
    allocate ( flow_xpencil ( domain_xpencil%irange(1) : domain_xpencil%irange(2), &
                              domain_xpencil%jrange(1) : domain_xpencil%jrange(2), &
                              domain_xpencil%krange(1) : domain_xpencil%krange(2) ) )
    

    allocate ( thermo_xpencil ( domain_xpencil%irange(1) : domain_xpencil%irange(2), &
                                domain_xpencil%jrange(1) : domain_xpencil%jrange(2), &
                                domain_xpencil%krange(1) : domain_xpencil%krange(2) ) )

  end subroutine Allocate_variables

end module flow_variables_mod