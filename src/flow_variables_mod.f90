module flow_variables_mod
  use precision_mod
  use input_general_mod, only: ithermo
  implicit none

  real(WP), save, allocatable, dimension(:, :, :) :: ux, uy, uz
  real(WP), save, allocatable, dimension(:, :, :) :: gx, gy, gz
  real(WP), save, allocatable, dimension(:, :, :) :: pre, phi

  real(WP), save, allocatable, dimension(:, :, :) :: massEnthalpy
  real(WP), save, allocatable, dimension(:, :, :) :: enthalpy
  real(WP), save, allocatable, dimension(:, :, :) :: density
  real(WP), save, allocatable, dimension(:, :, :) :: temperature
  real(WP), save, allocatable, dimension(:, :, :) :: thermalConductivity
  real(WP), save, allocatable, dimension(:, :, :) :: dynamicViscosity

  public :: Allocate_flow_variables

contains
  subroutine Allocate_flow_variables ()
    use decomp_2d

    call alloc_x (ux,  opt_global=.true.)
    call alloc_x (uy,  opt_global=.true.)
    call alloc_x (uz,  opt_global=.true.)

    call alloc_x (pre, opt_global=.true.)
    call alloc_x (phi, opt_global=.true.)

    call alloc_x (gx,  opt_global=.true.)
    call alloc_x (gy,  opt_global=.true.)
    call alloc_x (gz,  opt_global=.true.)

    if(ithermo == 1) then
      call alloc_x (massEnthalpy,        opt_global=.true.)
      call alloc_x (enthalpy,            opt_global=.true.)
      call alloc_x (density,             opt_global=.true.)
      call alloc_x (temperature,         opt_global=.true.)
      call alloc_x (thermalConductivity, opt_global=.true.)
      call alloc_x (dynamicViscosity,    opt_global=.true.)
    end if
    !if(nrank == 0) print *, shape(ux) !test

  end subroutine Allocate_flow_variables

end module flow_variables_mod