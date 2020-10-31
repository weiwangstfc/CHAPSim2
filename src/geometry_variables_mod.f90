module geometry_mod
  use input_mod
  use parameters_constant_mod
  use math_mod
  use VTK_mod, only : Generate_vtk_mesh_slice
  implicit none

  real(wp) :: dx, dz
  real(wp) :: dx2, dz2
  real(wp) :: dxi, dzi

  real(WP), allocatable, dimension(:) :: xc, xp
  real(WP), allocatable, dimension(:) :: yc, yp
  real(WP), allocatable, dimension(:) :: zc, zp

  public :: Initialize_geometry_variables

contains
  subroutine Initialize_geometry_variables ()
    use mpi_mod

    integer :: i, j, k
    real(WP) :: s, yy, c1, c2, c3, c4

    dx = lxx / real(ncx, WP)
    dz = lzz / real(ncz, WP)

    dx2 = dx * dx
    dz2 = dz * dz

    dxi = ONE / dx
    dzi = ONE / dz

    allocate ( xc(ncx) )
    allocate ( yc(ncy) )
    allocate ( zc(ncz) )
    xc(:) = ZERO
    yc(:) = ZERO
    zc(:) = ZERO

    allocate ( xp(npx) )
    allocate ( yp(npy) )
    allocate ( zp(npz) )
    xp(:) = ZERO
    yp(:) = ZERO
    zp(:) = ZERO

    block_xcoordinate: do i = 1, npx
      xp(i) = real( (i - 1), WP ) * dx
      if (i < npx) xc(i) = (real( i, WP ) - HALF) * dx
    end do block_xcoordinate

    block_zcoordinate: do k = 1, npz
      zp(k) = real( (k - 1), WP ) * dz
      if (k < npz) zc(k) = (real( k, WP ) - HALF) * dz
    end do block_zcoordinate

    block_ycnst: if (istret == ISTRET_SIDES) then
      c1 = rstret * HALF
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = HALF
    else if (istret == ISTRET_BOTTOM) then
      c1 = rstret * ONE
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = ONE
    else if (istret == ISTRET_TOP) then
      c1 = rstret * ZERO
      c2 = tanh_wp (rstret)
      c3 = ZERO
      c4 = ONE
    else
      c1 = rstret * HALF
      c2 = tanh_wp (c1)
      c3 = ONE
      c4 = HALF
    end if block_ycnst

    block_ynd: do j = 1, npy
      yy = real ((j - 1), WP) / real ( (npy - 1), WP)
      if(istret == ISTRET_NO) then
        s = yy
      else 
        s = (tanh_wp( (rstret * yy) - c1 ) / c2 + c3) * c4
      end if
      yp(j) = s * (lyt - lyb) + lyb
    end do block_ynd 

    block_ycl: do j = 1, ncy
      yc(j) = ( yp(j) + yp(j + 1) ) * HALF
    end do block_ycl

    if(nrank == 0) then
      call Generate_vtk_mesh_slice ( npx, npy, xp(:), yp(:), 'xy' )
      call Generate_vtk_mesh_slice ( npx, npz, xp(:), zp(:), 'xz' )
      call Generate_vtk_mesh_slice ( npy, npz, yp(:), zp(:), 'yz' )
    end if

  end subroutine  Initialize_geometry_variables
  
end module geometry_mod


