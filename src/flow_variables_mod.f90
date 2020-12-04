module flow_variables_mod
  use precision_mod
  use input_thermo_mod, only : thermoProperty_t
  implicit none

  type flow_t
    real(WP) :: ux
    real(WP) :: uy
    real(WP) :: uz

    real(WP) :: gx
    real(WP) :: gy
    real(WP) :: gz

    real(WP) :: pre
    real(WP) :: phi
  end type flow_t

  type(flow_t), save, allocatable, dimension(:, :, :) :: flow
  type(thermoProperty_t), save, allocatable, dimension(:, :, :) :: thermo


  public :: Allocate_flow_variables

contains
  subroutine Allocate_flow_variables ()
    use input_general_mod, only : ithermo
    use domain_decomposition_mod

    !allocate (flow (iStart_xpencil : iEnd_xpencil, &
    !                jStart_xpencil : jEnd_xpencil, &
    !                kStart_xpencil : kEnd_xpencil) )
    

  end subroutine Allocate_flow_variables

end module flow_variables_mod