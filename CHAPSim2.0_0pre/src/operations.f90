module operations
  implicit none

  private
  public :: Interpolate_stencil

contains

  subroutine Interpolate_1stencil_c2p(array0, array1, dir)
    real(WP), dimension(:, :, :), intent(in) :: array0
    real(WP), dimension(:, :, :), intent(in) :: array1
    integer :: dir
    

    ! x interpolation

    



  end subroutine Interpolate_stencil


end module