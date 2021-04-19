module flow_variables_mod
  use precision_mod
  use input_general_mod, only : ithermo
  implicit none

  real(WP), save, allocatable, dimension(:, :, :) :: u1_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: u2_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: u3_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: g1_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: g2_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: g3_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: pre_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: phi_xpencil

  real(WP), save, allocatable, dimension(:, :, :) :: d_xpencil
  real(WP), save, allocatable, dimension(:, :, :) :: m_xpencil
  real(WP), save, allocatable, dimension(:, :, :, :) :: dp_xpencil
  real(WP), save, allocatable, dimension(:, :, :, :) :: mp_xpencil

  if(ithermo == 1) then
    real(WP), save, allocatable, dimension(:, :, :) :: dh_xpencil
    real(WP), save, allocatable, dimension(:, :, :) :: h_xpencil
    real(WP), save, allocatable, dimension(:, :, :) :: k_xpencil
    real(WP), save, allocatable, dimension(:, :, :) :: t_xpencil
  end if
  

  private
  public :: Allocate_variables_xpencil

contains

  subroutine Allocate_variables_xpencil
    use input_general_mod, only : ithermo, NDIM
    use domain_decomposition_mod
    implicit none

    integer :: i0, i1
    integer :: j0, j1
    integer :: k0, k1

    
    i0 = domain_xpencil%irange(1)
    i1 = domain_xpencil%irange(2)

    j0 = domain_xpencil%jrange(1)
    j1 = domain_xpencil%jrange(2)

    k0 = domain_xpencil%krange(1)
    k1 = domain_xpencil%krange(2)

    ! allocate xpencil based data
    allocate ( u1_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; u1_xpencil = ZERO
    allocate ( u2_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; u2_xpencil = ZERO
    allocate ( u3_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; u3_xpencil = ZERO

    allocate ( g1_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; g1_xpencil = ZERO
    allocate ( g2_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; g2_xpencil = ZERO
    allocate ( g3_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; g3_xpencil = ZERO

    allocate ( pre_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; pre_xpencil = ZERO
    allocate ( phi_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; phi_xpencil = ZERO
    
    allocate ( d_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; d_xpencil = ONE
    allocate ( m_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; m_xpencil = ONE

    allocate ( dp_xpencil ( i0 : i1, j0 : j1, k0 : k1, NDIM )  ) ; dp_xpencil = ONE
    allocate ( mp_xpencil ( i0 : i1, j0 : j1, k0 : k1, NDIM )  ) ; mp_xpencil = ONE

    if(ithermo == 1) then
      allocate ( dh_xpencil ( i0 : i1, j0 : j1, k0 : k1 )  ) ; dh_xpencil = ZERO
      allocate ( h_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; h_xpencil = ZERO
      allocate ( k_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; k_xpencil = ONE
      allocate ( t_xpencil  ( i0 : i1, j0 : j1, k0 : k1 )  ) ; t_xpencil = ONE
    end if


  end subroutine Allocate_variables_xpencil

end module flow_variables_mod