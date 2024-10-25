module convert_primary_conservative_mod
  public  :: calculate_mflux_from_velo_bc
  public  :: calculate_mflux_from_velo_domain
  public  :: calcuate_velo_from_mflux_domain

 contains 

 subroutine calculate_mflux_from_velo_bc(fl, dm)
  use udf_type_mod
  use operations
  use decomp_2d
  use parameters_constant_mod
  implicit none
  type(t_domain), intent(inout) :: dm
  type(t_flow  ), intent(inout) :: fl

  real(WP), dimension( dm%dpcc%xsz(1), &
                       dm%dpcc%xsz(2), &
                       dm%dpcc%xsz(3)) :: d_pcc
  real(WP), dimension( dm%dcpc%ysz(1), &
                       dm%dcpc%ysz(2), &
                       dm%dcpc%ysz(3)) :: d_cpc_ypencil
  real(WP), dimension( dm%dccp%zsz(1), &
                       dm%dccp%zsz(2), &
                       dm%dccp%zsz(3)) :: d_ccp_zpencil
  real(WP), dimension( dm%dcpc%ysz(1), &
                       dm%dcpc%ysz(2), &
                       dm%dcpc%ysz(3)) :: qy_ypencil
  real(WP), dimension( dm%dcpc%ysz(1), &
                       dm%dcpc%ysz(2), &
                       dm%dcpc%ysz(3)) :: gy_ypencil
  real(WP), dimension( dm%dccp%ysz(1), &
                       dm%dccp%ysz(2), &
                       dm%dccp%ysz(3)) :: qz_ypencil
  real(WP), dimension( dm%dccp%ysz(1), &
                       dm%dccp%ysz(2), &
                       dm%dccp%ysz(3)) :: gz_ypencil
  real(WP), dimension( dm%dccc%ysz(1), &
                       dm%dccc%ysz(2), &
                       dm%dccc%ysz(3)) ::  d_ypencil
  real(WP), dimension( dm%dccp%zsz(1), &
                       dm%dccp%zsz(2), &
                       dm%dccp%zsz(3)) :: qz_zpencil
  real(WP), dimension( dm%dccp%zsz(1), &
                       dm%dccp%zsz(2), &
                       dm%dccp%zsz(3)) :: gz_zpencil
  real(WP), dimension( dm%dccc%zsz(1), &
                       dm%dccc%zsz(2), &
                       dm%dccc%zsz(3)) ::  d_zpencil

  real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_4cc 
  real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_c4c
  real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_cc4
  
  if(.not. dm%is_thermo) return
!----------------------------------------------------------------------------------------------------------
! x-pencil : u1 -> g1 = u1_pcc * d_pcc
!----------------------------------------------------------------------------------------------------------
  dm%fbcx_gx(:, :, :) = dm%fbcx_qx(:, :, :) * dm%fbcx_ftp(:, :, :)%d
  dm%fbcx_gy(:, :, :) = dm%fbcx_qy(:, :, :) * dm%fbcx_ftp(:, :, :)%d
!----------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2 = u2_cpc * d_cpc
!----------------------------------------------------------------------------------------------------------
  call transpose_x_to_y(fl%qy,    qy_ypencil, dm%dcpc)
  call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)
  fbcy_c4c(:, :, :) = dm%fbcy_ftp(:, :, :)%d
  call Get_y_midp_C2P_3D (d_ypencil, d_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)
  gy_ypencil = qy_ypencil * d_cpc_ypencil
  call transpose_y_to_x(gy_ypencil, fl%gy, dm%dcpc)
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3 = u3_ccp * d_ccp
!----------------------------------------------------------------------------------------------------------
  call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
  call transpose_x_to_y(fl%qz,      qz_ypencil, dm%dccp)
  call transpose_y_to_z(qz_ypencil, qz_zpencil, dm%dccp)
  fbcz_cc4(:, :, :) = dm%fbcz_ftp(:, :, :)%d
  call Get_z_midp_C2P_3D (d_zpencil, d_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp, fbcz_cc4)
  gz_zpencil = qz_zpencil * d_ccp_zpencil

  call transpose_z_to_y(gz_zpencil, gz_ypencil, dm%dccp)
  call transpose_y_to_x(gz_ypencil, fl%gz,      dm%dccp)

  return
end subroutine calculate_mflux_from_velo_bc

  subroutine calculate_mflux_from_velo_domain(fl, dm)
    use udf_type_mod
    use operations
    use decomp_2d
    use parameters_constant_mod
    implicit none
    type(t_domain), intent(inout)   :: dm
    type(t_flow  ), intent(inout) :: fl

    real(WP), dimension( dm%dpcc%xsz(1), &
                         dm%dpcc%xsz(2), &
                         dm%dpcc%xsz(3)) :: d_pcc
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: d_cpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: d_ccp_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: qy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: gy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: qz_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: gz_ypencil
    real(WP), dimension( dm%dccc%ysz(1), &
                         dm%dccc%ysz(2), &
                         dm%dccc%ysz(3)) ::  d_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: qz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: gz_zpencil
    real(WP), dimension( dm%dccc%zsz(1), &
                         dm%dccc%zsz(2), &
                         dm%dccc%zsz(3)) ::  d_zpencil

    real(WP), dimension( dm%dcpc%zsz(1), &
                         dm%dcpc%zsz(2), &
                         dm%dcpc%zsz(3)) ::  d_cpc_zpencil
    real(WP), dimension( dm%dcpc%xsz(1), &
                         dm%dcpc%xsz(2), &
                         dm%dcpc%xsz(3)) ::  d_cpc_xpencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) ::  d_ccp_ypencil
    real(WP), dimension( dm%dccp%xsz(1), &
                         dm%dccp%xsz(2), &
                         dm%dccp%xsz(3)) ::  d_ccp_xpencil

    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_4cc 
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_c4c
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_cc4
    
    if(.not. dm%is_thermo) return
!----------------------------------------------------------------------------------------------------------
! x-pencil : u1 -> g1 = u1_pcc * d_pcc
!----------------------------------------------------------------------------------------------------------
    fbcx_4cc(:, :, :) = dm%fbcx_ftp(:, :, :)%d
    call Get_x_midp_C2P_3D (fl%dDens, d_pcc, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4cc)
    fl%gx = fl%qx * d_pcc
!----------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2 = u2_cpc * d_cpc
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy,    qy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens,  d_ypencil, dm%dccc)
    fbcy_c4c(:, :, :) = dm%fbcy_ftp(:, :, :)%d
    call Get_y_midp_C2P_3D (d_ypencil, d_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)
    gy_ypencil = qy_ypencil * d_cpc_ypencil
    call transpose_y_to_x(gy_ypencil, fl%gy, dm%dcpc)
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3 = u3_ccp * d_ccp
!----------------------------------------------------------------------------------------------------------
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_x_to_y(fl%qz,      qz_ypencil, dm%dccp)
    call transpose_y_to_z(qz_ypencil, qz_zpencil, dm%dccp)
    fbcz_cc4(:, :, :) = dm%fbcz_ftp(:, :, :)%d
    call Get_z_midp_C2P_3D (d_zpencil, d_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp, fbcz_cc4)
    gz_zpencil = qz_zpencil * d_ccp_zpencil

    call transpose_z_to_y(gz_zpencil, gz_ypencil, dm%dccp)
    call transpose_y_to_x(gz_ypencil, fl%gz,      dm%dccp)

!----------------------------------------------------------------------------------------------------------
! BC: - x pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_qx(1) == IBC_DIRICHLET .or. dm%ibcx_qx(2) == IBC_DIRICHLET) then 
      dm%fbcx_gx(:, :, :) = dm%fbcx_qx(:, :, :) * dm%fbcx_ftp(:, :, :)%d
    end if

    if(dm%ibcx_qy(1) == IBC_DIRICHLET .or. dm%ibcx_qy(2) == IBC_DIRICHLET) then 
      call transpose_y_to_x(d_cpc_ypencil, d_cpc_xpencil, dm%dcpc)
      dm%fbcx_gy(1, :, :) = dm%fbcx_qy(1, :, :) * d_cpc_xpencil(1,              :, :)
      dm%fbcx_gy(2, :, :) = dm%fbcx_qy(2, :, :) * d_cpc_xpencil(dm%dcpc%xsz(1), :, :)
      dm%fbcx_gy(3, :, :) = dm%fbcx_gy(1, :, :)
      dm%fbcx_gy(4, :, :) = dm%fbcx_gy(2, :, :)
    end if

    if(dm%ibcx_qz(1) == IBC_DIRICHLET .or. dm%ibcx_qz(2) == IBC_DIRICHLET) then 
      call transpose_z_to_y(d_ccp_zpencil, d_ccp_ypencil, dm%dccp)
      call transpose_y_to_x(d_ccp_ypencil, d_ccp_xpencil, dm%dccp)
      dm%fbcx_gz(1, :, :) = dm%fbcx_qz(1, :, :) * d_ccp_xpencil(1,              :, :)
      dm%fbcx_gz(2, :, :) = dm%fbcx_qz(2, :, :) * d_ccp_xpencil(dm%dcpc%xsz(1), :, :)
      dm%fbcx_gz(3, :, :) = dm%fbcx_gz(1, :, :)
      dm%fbcx_gz(4, :, :) = dm%fbcx_gz(2, :, :)
    end if
!----------------------------------------------------------------------------------------------------------
! BC: - y pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qx(1) == IBC_DIRICHLET .or. dm%ibcy_qx(2) == IBC_DIRICHLET) then 
      dm%fbcy_gx(:, :, :) = dm%fbcy_qx(:, :, :) * dm%fbcy_ftp(:, :, :)%d
    end if

    if(dm%ibcy_qy(1) == IBC_DIRICHLET .or. dm%ibcy_qy(1) == IBC_DIRICHLET) then 
      dm%fbcy_gy(:, 1, :) = dm%fbcy_qy(:, 1, :) * d_cpc_ypencil(:, 1,              :)
      dm%fbcy_gy(:, 2, :) = dm%fbcy_qy(:, 2, :) * d_cpc_ypencil(:, dm%dcpc%xsz(1), :)
      dm%fbcy_gy(:, 3, :) = dm%fbcy_gy(:, 1, :)
      dm%fbcy_gy(:, 4, :) = dm%fbcy_gy(:, 2, :)
    end if

    if(dm%ibcy_qz(1) == IBC_DIRICHLET .or. dm%ibcy_qz(1) == IBC_DIRICHLET) then 
      call transpose_z_to_y(d_ccp_zpencil, d_ccp_ypencil, dm%dccp)
      dm%fbcy_gz(:, 1, :) = dm%fbcy_qz(:, 1, :) * d_ccp_ypencil(:, 1,              :)
      dm%fbcy_gz(:, 2, :) = dm%fbcy_qz(:, 2, :) * d_ccp_ypencil(:, dm%dccp%xsz(1), :)
      dm%fbcy_gz(:, 3, :) = dm%fbcy_gz(:, 1, :)
      dm%fbcy_gz(:, 4, :) = dm%fbcy_gz(:, 2, :)
    end if
!----------------------------------------------------------------------------------------------------------
! BC: - z pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcz_qx(1) == IBC_DIRICHLET .or. dm%ibcz_qx(1) == IBC_DIRICHLET) then 
      dm%fbcz_gx(:, :, :) = dm%fbcz_qx(:, :, :) * dm%fbcz_ftp(:, :, :)%d
    end if

    if(dm%ibcz_qy(1) == IBC_DIRICHLET .or. dm%ibcz_qy(1) == IBC_DIRICHLET) then 
      call transpose_y_to_x(d_cpc_ypencil, d_cpc_zpencil, dm%dcpc)
      dm%fbcz_gy(:, :, 1) = dm%fbcz_qy(:, :, 1) * d_cpc_zpencil(:, :, 1             )
      dm%fbcz_gy(:, :, 2) = dm%fbcz_qy(:, :, 2) * d_cpc_zpencil(:, :, dm%dcpc%xsz(1))
      dm%fbcz_gy(:, :, 3) = dm%fbcz_gy(:, :, 1)
      dm%fbcz_gy(:, :, 4) = dm%fbcz_gy(:, :, 2)
    end if

    if(dm%ibcz_qz(1) == IBC_DIRICHLET .or. dm%ibcz_qz(2) == IBC_DIRICHLET) then 
      dm%fbcz_gz(:, :, 1) = dm%fbcz_qz(:, :, 1) * d_ccp_zpencil(:, :, 1             )
      dm%fbcz_gz(:, :, 2) = dm%fbcz_qz(:, :, 2) * d_ccp_zpencil(:, :, dm%dccp%xsz(1))
      dm%fbcz_gz(:, :, 3) = dm%fbcz_gz(:, :, 1)
      dm%fbcz_gz(:, :, 4) = dm%fbcz_gz(:, :, 2)
    end if


    return
  end subroutine calculate_mflux_from_velo_domain
  !==========================================================================================================
  subroutine calcuate_velo_from_mflux_domain(fl, dm)
    use udf_type_mod
    use operations
    use decomp_2d
    use parameters_constant_mod
    implicit none
    type(t_domain), intent(inout)   :: dm
    type(t_flow  ), intent(inout) :: fl

    real(WP), dimension( dm%dpcc%xsz(1), &
                         dm%dpcc%xsz(2), &
                         dm%dpcc%xsz(3)) :: d_pcc
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: d_cpc_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: d_ccp_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: qy_ypencil
    real(WP), dimension( dm%dcpc%ysz(1), &
                         dm%dcpc%ysz(2), &
                         dm%dcpc%ysz(3)) :: gy_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: qz_ypencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) :: gz_ypencil
    real(WP), dimension( dm%dccc%ysz(1), &
                         dm%dccc%ysz(2), &
                         dm%dccc%ysz(3)) ::  d_ypencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: qz_zpencil
    real(WP), dimension( dm%dccp%zsz(1), &
                         dm%dccp%zsz(2), &
                         dm%dccp%zsz(3)) :: gz_zpencil
    real(WP), dimension( dm%dccc%zsz(1), &
                         dm%dccc%zsz(2), &
                         dm%dccc%zsz(3)) ::  d_zpencil
    real(WP), dimension( dm%dcpc%zsz(1), &
                         dm%dcpc%zsz(2), &
                         dm%dcpc%zsz(3)) ::  d_cpc_zpencil
    real(WP), dimension( dm%dcpc%xsz(1), &
                         dm%dcpc%xsz(2), &
                         dm%dcpc%xsz(3)) ::  d_cpc_xpencil
    real(WP), dimension( dm%dccp%ysz(1), &
                         dm%dccp%ysz(2), &
                         dm%dccp%ysz(3)) ::  d_ccp_ypencil
    real(WP), dimension( dm%dccp%xsz(1), &
                         dm%dccp%xsz(2), &
                         dm%dccp%xsz(3)) ::  d_ccp_xpencil

    real(WP), dimension( 4, dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: fbcx_4cc 
    real(WP), dimension( dm%dcpc%ysz(1), 4, dm%dcpc%ysz(3) ) :: fbcy_c4c
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), 4 ) :: fbcz_cc4
    
    if(.not. dm%is_thermo) return
!----------------------------------------------------------------------------------------------------------
! x-pencil : g1 -> u1
!----------------------------------------------------------------------------------------------------------
    fbcx_4cc(:, :, :) = dm%fbcx_ftp(:, :, :)%d
    call Get_x_midp_C2P_3D (fl%dDens, d_pcc, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4cc) ! check, d_pcc is not changed in mom eqs, could be save seperatedly to save calc.
    fl%qx = fl%gx / d_pcc
!----------------------------------------------------------------------------------------------------------
! y-pencil : u2 -> g2
!----------------------------------------------------------------------------------------------------------
    fbcy_c4c(:, :, :) = dm%fbcy_ftp(:, :, :)%d
    call transpose_x_to_y(fl%gy,   gy_ypencil, dm%dcpc)
    call transpose_x_to_y(fl%dDens, d_ypencil, dm%dccc)
    call Get_y_midp_C2P_3D (d_ypencil, d_cpc_ypencil, dm, dm%iAccuracy, dm%ibcy_ftp, fbcy_c4c)
    qy_ypencil = gy_ypencil / d_cpc_ypencil
    call transpose_y_to_x(qy_ypencil, fl%qy, dm%dcpc)
!----------------------------------------------------------------------------------------------------------
! Z-pencil : u3 -> g3
!----------------------------------------------------------------------------------------------------------
    fbcz_cc4(:, :, :) = dm%fbcz_ftp(:, :, :)%d
    call transpose_y_to_z( d_ypencil,  d_zpencil, dm%dccc)
    call transpose_x_to_y(fl%gz,      gz_ypencil, dm%dccp)
    call transpose_y_to_z(gz_ypencil, gz_zpencil, dm%dccp)
    call Get_z_midp_C2P_3D (d_zpencil, d_ccp_zpencil, dm, dm%iAccuracy, dm%ibcz_ftp, fbcz_cc4)
    qz_zpencil = gz_zpencil / d_ccp_zpencil
    call transpose_z_to_y(qz_zpencil, qz_ypencil, dm%dccp)
    call transpose_y_to_x(qz_ypencil, fl%qz, dm%dccp)

    fl%gx0 = fl%gx
    fl%gy0 = fl%gy
    fl%gz0 = fl%gz


!----------------------------------------------------------------------------------------------------------
! BC: - x pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcx_qx(1) == IBC_DIRICHLET .and. dm%ibcx_qx(2) == IBC_DIRICHLET) then 
      dm%fbcx_qx(:, :, :) = dm%fbcx_gx(:, :, :) / dm%fbcx_ftp(:, :, :)%d
    end if

    if(dm%ibcx_qy(1) == IBC_DIRICHLET .or. dm%ibcx_qy(2) == IBC_DIRICHLET) then 
      call transpose_y_to_x(d_cpc_ypencil, d_cpc_xpencil, dm%dcpc)
      dm%fbcx_qy(1, :, :) = dm%fbcx_gy(1, :, :) / d_cpc_xpencil(1,              :, :)
      dm%fbcx_qy(2, :, :) = dm%fbcx_gy(2, :, :) / d_cpc_xpencil(dm%dcpc%xsz(1), :, :)
      dm%fbcx_qy(3, :, :) = dm%fbcx_qy(1, :, :)
      dm%fbcx_qy(4, :, :) = dm%fbcx_qy(2, :, :)
    end if

    if(dm%ibcx_qz(1) == IBC_DIRICHLET .or. dm%ibcx_qz(2) == IBC_DIRICHLET) then 
      call transpose_z_to_y(d_ccp_zpencil, d_ccp_ypencil, dm%dccp)
      call transpose_y_to_x(d_ccp_ypencil, d_ccp_xpencil, dm%dccp)
      dm%fbcx_qz(1, :, :) = dm%fbcx_gz(1, :, :) / d_ccp_xpencil(1,              :, :)
      dm%fbcx_qz(2, :, :) = dm%fbcx_gz(2, :, :) / d_ccp_xpencil(dm%dcpc%xsz(1), :, :)
      dm%fbcx_qz(3, :, :) = dm%fbcx_qz(1, :, :)
      dm%fbcx_qz(4, :, :) = dm%fbcx_qz(2, :, :)
    end if
!----------------------------------------------------------------------------------------------------------
! BC: - y pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcy_qx(1) == IBC_DIRICHLET .or. dm%ibcy_qx(2) == IBC_DIRICHLET) then 
      dm%fbcy_qx(:, :, :) = dm%fbcy_gx(:, :, :) / dm%fbcy_ftp(:, :, :)%d
    end if

    if(dm%ibcy_qy(1) == IBC_DIRICHLET .or. dm%ibcy_qy(2) == IBC_DIRICHLET) then 
      dm%fbcy_qy(:, 1, :) = dm%fbcy_gy(:, 1, :) / d_cpc_ypencil(:, 1,              :)
      dm%fbcy_qy(:, 2, :) = dm%fbcy_gy(:, 2, :) / d_cpc_ypencil(:, dm%dcpc%xsz(1), :)
      dm%fbcy_qy(:, 3, :) = dm%fbcy_qy(:, 1, :)
      dm%fbcy_qy(:, 4, :) = dm%fbcy_qy(:, 2, :)
    end if

    if(dm%ibcy_qz(1) == IBC_DIRICHLET .or. dm%ibcy_qz(2) == IBC_DIRICHLET) then 
      call transpose_z_to_y(d_ccp_zpencil, d_ccp_ypencil, dm%dccp)
      dm%fbcy_qz(:, 1, :) = dm%fbcy_gz(:, 1, :) / d_ccp_ypencil(:, 1,              :)
      dm%fbcy_qz(:, 2, :) = dm%fbcy_gz(:, 2, :) / d_ccp_ypencil(:, dm%dccp%xsz(1), :)
      dm%fbcy_qz(:, 3, :) = dm%fbcy_qz(:, 1, :)
      dm%fbcy_qz(:, 4, :) = dm%fbcy_qz(:, 2, :)
    end if
!----------------------------------------------------------------------------------------------------------
! BC: - z pencil
!----------------------------------------------------------------------------------------------------------
    if(dm%ibcz_qx(1) == IBC_DIRICHLET .or. dm%ibcz_qx(2) == IBC_DIRICHLET) then 
      dm%fbcz_qx(:, :, :) = dm%fbcz_gx(:, :, :) / dm%fbcz_ftp(:, :, :)%d
    end if

    if(dm%ibcz_qy(1) == IBC_DIRICHLET .or. dm%ibcz_qy(2) == IBC_DIRICHLET) then 
      call transpose_y_to_x(d_cpc_ypencil, d_cpc_zpencil, dm%dcpc)
      dm%fbcz_qy(:, :, 1) = dm%fbcz_gy(:, :, 1) / d_cpc_zpencil(:, :, 1             )
      dm%fbcz_qy(:, :, 2) = dm%fbcz_gy(:, :, 2) / d_cpc_zpencil(:, :, dm%dcpc%xsz(1))
      dm%fbcz_qy(:, :, 3) = dm%fbcz_qy(:, :, 1)
      dm%fbcz_qy(:, :, 4) = dm%fbcz_qy(:, :, 2)
    end if

    if(dm%ibcz_qz(1) == IBC_DIRICHLET .or. dm%ibcz_qz(2) == IBC_DIRICHLET) then 
      dm%fbcz_qz(:, :, 1) = dm%fbcz_gz(:, :, 1) / d_ccp_zpencil(:, :, 1             )
      dm%fbcz_qz(:, :, 2) = dm%fbcz_gz(:, :, 2) / d_ccp_zpencil(:, :, dm%dccp%xsz(1))
      dm%fbcz_qz(:, :, 3) = dm%fbcz_qz(:, :, 1)
      dm%fbcz_qz(:, :, 4) = dm%fbcz_qz(:, :, 2)
    end if

    return
  end subroutine calcuate_velo_from_mflux_domain

end module