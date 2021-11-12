





!-------------------------------------------------------------------------------
!    qx --> qx_ypencil --> qx_yppc_ypencil --> qx_yppc
!                     |--> qx_zpencil --> qx_zpcp_zpencil --> qx_zpcp_ypencil --> qx_zpcp
!------------------------------------------------------------------------------- 
    call transpose_x_to_y(f%qx, qx_ypencil, dm%dpcc)                 ! used in x-mom, w   thermal
    call Get_y_midp_C2P_3dArray ( qx_ypencil, d, qx_yppc_ypencil )  ! used in x-mom, w/o thermal
    call transpose_y_to_x(qx_yppc_ypencil, qx_yppc, dm%dppc)         ! used in x-mom, wo  thermal

    call transpose_y_to_z(qx_ypencil, qx_zpencil, dm%dpcc)           ! used in x-mom, w thermal
    call Get_z_midp_C2P_3dArray ( qx_zpencil, d, qx_zpcp_zpencil )  ! used in x-mon, w/o thermal
    call transpose_z_to_y(qx_zpcp_zpencil, qx_zpcp_ypencil, dm%dpcp) ! intermediate, not used.
    call transpose_y_to_x(qx_zpcp_ypencil, qx_zpcp,         dm%dpcp) ! used in z-mom, wo  thermal

!-------------------------------------------------------------------------------
!    qy --> qy_xppc --> qy_xppc_ypencil
!      |--> qy_ypencil --> qy_zpencil --> qy_zcpp_zpencil --> qy_zcpp_ypencil
!------------------------------------------------------------------------------- 

    call Get_x_midp_C2P_3dArray ( f%qy, d, qy_xppc )                ! used in mom-x, w/o thermal
    call transpose_x_to_y(f%qy, qy_ypencil, dm%dcpc)                 ! used in y-mom, w/o thermal 
    call transpose_y_to_z(qy_ypencil, qy_zpencil, dm%dcpc)           ! used in y-mom, w   thermal
    call Get_z_midp_C2P_3dArray ( qy_zpencil, d, qy_zcpp_zpencil )  ! used in y-mon, w   thermal
    call transpose_z_to_y(qy_zcpp_zpencil, qy_zcpp_ypencil, dm%dcpp) ! used in z-mom, wo  thermal




    
    if(ithermo == 1) then
      call Get_x_midp_P2C_3dArray ( f%qx,    d, qx_xccc )
      call Get_x_midp_P2C_3dArray ( f%gx,    d, gx_xccc )
      call Get_x_midp_C2P_3dArray ( f%gy,    d, gy_xppc )
      call Get_x_midp_C2P_3dArray ( f%gz,    d, gz_xpcp )
      call Get_x_midp_C2P_3dArray ( f%mVisc, d,  m_xpcc )
    end if

    

!-------------------------------------------------------------------------------
!    qz --> qz_xpcp --> qz_xpcp_ypencil --> qz_xpcp_zpencil
!      |--> qz_ypencil --> qz_zpencil --> qz_zccc_zpencil --> qz_zccc_ypencil --> qz_zccc
!                     |--> qz_ycpp_ypencil
      
!-------------------------------------------------------------------------------
    call Get_x_midp_C2P_3dArray ( f%qz, d, qz_xpcp )                ! used in z-mom, wo  thermal
    call transpose_x_to_y(qz_xpcp,         qz_xpcp_ypencil, dm%dpcp) ! intermediate, not used.
    call transpose_y_to_z(qz_xpcp_ypencil, qz_xpcp_zpencil, dm%dpcp) ! used in z-mom, wo  thermal

    call transpose_x_to_y(f%qz, qz_ypencil, dm%dccp)                   ! used in z-mom, w   thermal
    call transpose_y_to_z(qz_ypencil, qz_zpencil, dm%dccp)             ! used in z-mom, w   thermal
    call Get_y_midp_C2P_3dArray ( qz_ypencil, d, qz_ycpp_ypencil )    ! used in mom-z, w/o thermal
    if(ithermal == 1) then
      call Get_z_midp_P2C_3dArray ( qz_zpencil, d, qz_zccc_zpencil )  ! intermediate, not used.
      call transpose_z_to_y(qz_zccc_zpencil, qz_zccc_ypencil, dm%dccc) ! used in y-mom, w   thermal
      call transpose_y_to_x(qz_zccc_ypencil, qz_zccc,         dm%dccc) ! used in x-mom, w   thermal
    end if
    
    if(ithermo == 1) then
      call transpose_x_to_y(f%gx,   gx_ypencil, dm%dpcc)
      call Get_y_midp_C2P_3dArray ( gx_ypencil, d, gx_yppc_ypencil )

      call transpose_x_to_y(f%gy,   gy_ypencil, dm%dcpc)
      call Get_y_midp_P2C_3dArray ( gy_ypencil, d, gy_yccc_ypencil )

      call transpose_x_to_y(f%gz,   gz_ypencil, dm%dccp)
      call Get_y_midp_C2P_3dArray ( gz_ypencil, d, gz_ycpp_ypencil )
      
      call transpose_x_to_y(f%mVisc, m_ypencil, dm%dccc)
      call Get_y_midp_C2P_3dArray (  m_ypencil, d,  m_ycpc_ypencil )

      call transpose_x_to_y(gy_xppc, gy_xppc_ypencil, dm%dppc)
      call transpose_x_to_y(qx_xccc, qx_xccc_ypencil, dm%dccc)

      call Get_y_midp_P2C_3dArray ( qy_ypencil, d, qy_yccc_ypencil )! used in z-mom, w thermal
      call transpose_y_to_x(qy_yccc_ypencil, qy_yccc, dm%dccc) ! used in x-mom, w thermal

      call transpose_y_to_x(gx_yppc_ypencil, gx_yppc, dm%dppc)
    end if

    
    
    

    
    
    
    if(ithermo == 1) then
      call transpose_y_to_z(gx_ypencil, gx_zpencil, dm%dpcc)
      call transpose_y_to_z(gy_ypencil, gy_zpencil, dm%dcpc)
      call transpose_y_to_z(gz_ypencil, gz_zpencil, dm%dccp)
      call transpose_y_to_z( m_ypencil,  m_zpencil, dm%dccc)

      call Get_z_midp_C2P_3dArray ( gx_zpencil, d, gx_zpcp_zpencil )
      call Get_z_midp_C2P_3dArray ( gy_zpencil, d, gy_zcpp_zpencil )
      call Get_z_midp_P2C_3dArray ( gz_zpencil, d, gz_zccc_zpencil )
      c

      call transpose_z_to_y(gx_zpcp_zpencil, gx_zpcp_ypencil, dm%dpcp)
      call transpose_z_to_y(gx_zpcp_ypencil, gx_zpcp,         dm%dpcp)
      call transpose_z_to_y(gy_zcpp_zpencil, gy_zcpp_ypencil, dm%dcpp)

      call Get_z_midp_P2C_3dArray ( qz_zpencil, d, qz_zccc_zpencil )  ! intermediate, not used.
      call transpose_z_to_y(qz_zccc_zpencil, qz_zccc_ypencil, dm%dccc) ! used in y-mom, w  thermal
      call transpose_y_to_x(qz_zccc_ypencil, qz_zccc,         dm%dccc) ! used in x-mom, wo thermal

    end if
!_______________________________________________________________________________
! 1st Deriviation : dmdx at points (preparation)
!_______________________________________________________________________________
    if(ithermo == 1) then
!-------------------------------------------------------------------------------
! In x-pencil, dm_ycpc/dx & dm_zccp/dx & d(qx)/dx,  operation in x direction
!-------------------------------------------------------------------------------
      call transpose_z_to_y(m_zccp_zpencil, m_zccp_ypencil, dm%dccp)
      call transpose_y_to_x(m_zccp_ypencil, m_zccp,         dm%dccp)
      call transpose_y_to_x(m_ycpc_ypencil, m_ycpc,         dm%dcpc)

      call Get_x_1st_derivative_C2P_3dArray( f%mVisc, d, dmdx_xpcc )
      call Get_x_1st_derivative_C2C_3dArray(  m_ycpc, d, dmdx_ycpc )
      call Get_x_1st_derivative_C2C_3dArray(  m_zccp, d, dmdx_zccp )
      call Get_x_1st_derivative_P2C_3dArray(    f%qx, d, div0      )
      div = div + div0
!-------------------------------------------------------------------------------
! In Y-pencil, dm_xpcc/dy & dm_zccp/dy & d(qy)/dy, operation in y direction
!-------------------------------------------------------------------------------
      call transpose_x_to_y(m_xpcc,       m_xpcc_ypencil, dm%dpcc)
      call transpose_x_to_y(dmdx_ycpc, dmdx_ycpc_ypencil, dm%dcpc)
      call Get_y_1st_derivative_C2C_3dArray( m_xpcc_ypencil, d, dmdy_xpcc_ypencil )
      call Get_y_1st_derivative_C2P_3dArray( m_ypencil,      d, dmdy_ycpc_ypencil )
      call Get_y_1st_derivative_C2C_3dArray( m_zccp_ypencil, d, dmdy_zccp_ypencil )
      call Get_y_1st_derivative_P2C_3dArray( qy_ypencil,     d,       div_ypencil )
      call transpose_y_to_x(div_ypencil,       div0,      dm%dccc)
      call transpose_y_to_x(dmdy_xpcc_ypencil, dmdy_xpcc, dm%dpcc)
      div = div + div0
!-------------------------------------------------------------------------------
! In Z-pencil, dm_xpcc/dz & dm_ycpc/dz operation in z direction
!-------------------------------------------------------------------------------
      call transpose_y_to_z(m_xpcc_ypencil, m_xpcc_zpencil, dm%dpcc)
      call transpose_y_to_z(m_ycpc_ypencil, m_ycpc_zpencil, dm%dcpc)
      call Get_z_1st_derivative_C2C_3dArray( m_xpcc_zpencil, d, dmdz_xpcc_zpencil )
      call Get_z_1st_derivative_C2C_3dArray( m_ycpc_zpencil, d, dmdz_ycpc_zpencil )
      call Get_z_1st_derivative_C2P_3dArray(      m_zpencil, d, dmdz_zccp_zpencil )
      call Get_z_1st_derivative_P2C_3dArray(     qz_zpencil, d,       div_zpencil )

      call transpose_z_to_y(dmdz_ycpc_zpencil, dmdz_ycpc_ypencil, dm%dcpc)

      call transpose_z_to_y(div_zpencil, div_ypencil, dm%dccc)
      call transpose_y_to_x(div_ypencil, div0,        dm%dccc)
      div = div + div0
      call transpose_x_to_y(div,         div_ypencil, dm%dccc)
      call transpose_y_to_z(div_ypencil, div_zpencil, dm%dccc)
    end if