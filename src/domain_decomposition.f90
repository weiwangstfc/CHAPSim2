!##############################################################################
module domain_decompistion_mod
  implicit none

  private
  public :: Initialize_domain_decompsition

contains

  subroutine Initialize_domain_decompsition ()
    use input_general_mod, only : npx, npy, npz, p_row, p_col
    use decomp_2d

    call decomp_2d_init( npx, npy, npz, p_row, p_col )

  end subroutine Initialize_domain_decompsition

end module domain_decompistion_mod