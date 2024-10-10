module convert_primary_conservative_mod

  public  :: calculate_mflux_from_velo_domain
  public  :: calcuate_velo_from_mflux_domain

 contains 

  subroutine calculate_mflux_from_velo_domain(fl, dm)
    use udf_type_mod
    use operations
    use decomp_2d
    use parameters_constant_mod
    implicit none
    type(t_domain), intent(in )   :: dm
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

    return
  end subroutine calculate_mflux_from_velo_domain
  !==========================================================================================================
  subroutine calcuate_velo_from_mflux_domain(fl, dm)
    use udf_type_mod
    use operations
    use decomp_2d
    use parameters_constant_mod
    implicit none
    type(t_domain), intent(in )   :: dm
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
! x-pencil : g1 -> u1
!----------------------------------------------------------------------------------------------------------
    fbcx_4cc(:, :, :) = dm%fbcx_ftp(:, :, :)%d
    call Get_x_midp_C2P_3D (fl%dDens, d_pcc, dm, dm%iAccuracy, dm%ibcx_ftp, fbcx_4cc)
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

    return
  end subroutine calcuate_velo_from_mflux_domain

end module