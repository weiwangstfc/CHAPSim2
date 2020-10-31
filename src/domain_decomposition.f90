!##############################################################################
subroutine Initialize_domain_decompsition ()
  use input_mod, only : npx, npy, npz, p_row, p_col
  use decomp_2d
  implicit none

  call decomp_2d_init( npx, npy, npz, p_row, p_col )

end subroutine Initialize_domain_decompsition
