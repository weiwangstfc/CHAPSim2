module geometry_variables_mod
  use precision_mod
  use parameters_input_mod
  use parameters_constant_mod
  use math_mod

  real(wp) :: dx, dz
  real(wp) :: dx2, dz2
  real(wp) :: dxi, dzi

  real(WP), allocatable, dimension(:) :: xc, xp
  real(WP), allocatable, dimension(:) :: yc, yp
  real(WP), allocatable, dimension(:) :: zc, zp


contains
  subroutine Initialize_geometry_variables ()
    implicit none
    integer :: i, j, k
    real(WP) :: s, yy

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


    block_xcoordinate: do i = 1, npx - 1
      xp(i) = real( (i - 1), WP ) * dx
      xc(i) = (real( i, WP ) - HALF) * dx
    end do block_xcoordinate
    if(is_x_periodic) then
      xp(npx) = xp(1)
      xc(ncx) = xp(npx) + HALF * dx
    else 
      xp(npx) = lxx
    end if 

    block_zcoordinate: do k = 1, npz - 1
      zp(k) = real( (k - 1), WP ) * dz
      zc(k) = zp(i) + HALF * dz
    end do block_zcoordinate
    if(is_z_periodic) then
      zp(npz) = zp(1)
    else 
      zp(npz) = lzz
    end if 

    block_ycoordinate: if( istret == ISTRET_NO ) then

      block_uniform_y: do j = 1, npy
        s = real ((j - 1), WP) / real ( (npy - 1), WP)
        yp(j) = s * (lyt - lyb)  + lyb
      end do block_uniform_y

    else if (istret == ISTRET_CENTRE) then

      yp(1) = lyb
      yp(npy) = lyt

    else if (istret == ISTRET_SIDES ) then
      
      yp(1) = lyb
      yp(npy) = lyt

      do j = 2, npy - 1
        yy = real ((j - 1), WP) / real ( (npy - 1), WP)
        s = ONE + (tanh_wp (rstret * (yy - ONE ) ) ) / (tanh_wp ( rstret ) )
        yp(j) = s * (yp(npy) - yp(1)) + yp(1)
      end do

    else if (istret == ISTRET_BOTTOM ) then 
    else if (istret == ISTRET_TOP ) then 
    else 
    end if block_ycoordinate


  end subroutine  Initialize_geometry_variables
  
end module geometry_variables_mod


