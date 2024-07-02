module eq_vof_mod
  use precision_mod
  use operations
  use decomp_2d
  implicit none

  real(WP), parameter :: eps = 1.e-16_WP

!  public :: cal_colour_function
  private :: pxyz
  private :: itg_h_hat
  private :: cal_flux
  private :: clip_vof
  public  :: Solve_vof_eq
  public :: cal_colour_function

  contains

!===============================================================================
!===============================================================================
  !calculate the colour function
  subroutine cal_colour_function(dm, vf)
    use parameters_constant_mod
    use operations
    use udf_type_mod

    implicit none

    type(t_domain), intent(IN)    :: dm
    type(t_vof),    intent(INOUT) :: vf

    real(WP),dimension(dm%dccc%xsz(1),dm%dccc%xsz(2),dm%dccc%xsz(3)) ::        &
                                                                   mx_cccx,    &
                                                                   my_cccx,    &
                                                                   mz_cccx,    &
                                                                   modcccx

    real(WP),dimension(dm%dppp%xsz(1),dm%dppp%xsz(2),dm%dppp%xsz(3)) ::        &
                                                                   mx_pppx,    &
                                                                   my_pppx,    &
                                                                   mz_pppx,    &
                                                                   nx_pppx,    &
                                                                   ny_pppx,    &
                                                                   nz_pppx,    &
                                                                   modpppx

    real(WP),dimension(dm%dccc%xsz(1),dm%dccc%xsz(2),dm%dccc%xsz(3)) :: cccx
    real(WP),dimension(dm%dccc%ysz(1),dm%dccc%ysz(2),dm%dccc%ysz(3)) :: cccy
    real(WP),dimension(dm%dccc%zsz(1),dm%dccc%zsz(2),dm%dccc%zsz(3)) :: cccz
    real(WP),dimension(dm%dpcc%xsz(1),dm%dpcc%xsz(2),dm%dpcc%xsz(3)) :: pccx
    real(WP),dimension(dm%dpcc%ysz(1),dm%dpcc%ysz(2),dm%dpcc%ysz(3)) :: pccy
    real(WP),dimension(dm%dpcc%zsz(1),dm%dpcc%zsz(2),dm%dpcc%zsz(3)) :: pccz
    real(WP),dimension(dm%dcpc%ysz(1),dm%dcpc%ysz(2),dm%dcpc%ysz(3)) :: cpcy
    real(WP),dimension(dm%dcpc%zsz(1),dm%dcpc%zsz(2),dm%dcpc%zsz(3)) :: cpcz
    real(WP),dimension(dm%dccp%ysz(1),dm%dccp%ysz(2),dm%dccp%ysz(3)) :: ccpy
    real(WP),dimension(dm%dccp%zsz(1),dm%dccp%zsz(2),dm%dccp%zsz(3)) :: ccpz
    real(WP),dimension(dm%dcpp%xsz(1),dm%dcpp%xsz(2),dm%dcpp%xsz(3)) :: cppx
    real(WP),dimension(dm%dcpp%ysz(1),dm%dcpp%ysz(2),dm%dcpp%ysz(3)) :: cppy
    real(WP),dimension(dm%dcpp%zsz(1),dm%dcpp%zsz(2),dm%dcpp%zsz(3)) :: cppz
    real(WP),dimension(dm%dpcp%ysz(1),dm%dpcp%ysz(2),dm%dpcp%ysz(3)) :: pcpy
    real(WP),dimension(dm%dpcp%zsz(1),dm%dpcp%zsz(2),dm%dpcp%zsz(3)) :: pcpz
    real(WP),dimension(dm%dppc%ysz(1),dm%dppc%ysz(2),dm%dppc%ysz(3)) :: ppcy
    real(WP),dimension(dm%dppc%zsz(1),dm%dppc%zsz(2),dm%dppc%zsz(3)) :: ppcz
    real(WP),dimension(dm%dppp%ysz(1),dm%dppp%ysz(2),dm%dppp%ysz(3)) :: pppy
    real(WP),dimension(dm%dppp%zsz(1),dm%dppp%zsz(2),dm%dppp%zsz(3)) :: pppz

    real(WP) :: fbcx(4, dm%np(2), dm%np(3))
    real(WP) :: fbcy(dm%np(1), 4, dm%np(3))
    real(WP) :: fbcz(dm%np(1), dm%np(2), 4)
    integer  :: ibcx(2)
    integer  :: ibcy(2)
    integer  :: ibcz(2)

    integer  :: n, i, j, k

    real(WP) :: hh, beta
    real(WP) :: lnx, lny, lnz, llxx, llyy, llzz, llxy, llxz, llyz, phi
    real(WP) :: y1, y2
    real(WP) :: a200, a020, a002, a110, a011, a101, a100, a010, a001
    real(WP) :: cx, cy, cz, A, B1, B2, C1, C2, Q
    real(WP) :: aa0, aa1, aa2, aa3, aa4
    real(WP) :: D

    real(WP) :: A2d, B12d, B22d, Q2d, aa2d, bb2d, cc2d, solu1, solu2, D2d

    real(WP) :: limit

    hh = dm%h(1)
    beta = real(vf%ibeta,WP)

    !calculate mx on x-pencil
    if(vf%igrad==1) then
      call Get_x_midp_C2P_3D(vf%phi, pccx, dm, dm%ibcx(:,6),                   &
                             dm%fbcx_vof(:,:,:))
      call Get_x_1st_derivative_P2C_3D(pccx, mx_cccx, dm, dm%ibcx(:,6),        &
                                       dm%fbcx_vof(:,:,:))
    else if (vf%igrad==2) then
      call transpose_x_to_y(vf%phi, cccy, dm%dccc)
      call Get_y_midp_C2P_3D(cccy, cpcy, dm, dm%ibcy(:,6),                     &
                             dm%fbcy_vof(:,:,:))
      call transpose_y_to_z(cpcy, cpcz, dm%dcpc)
      call Get_z_midp_C2P_3D(cpcz, cppz, dm, dm%ibcz(:,6),                     &
                             dm%fbcz_vof(:,:,:))
      call transpose_z_to_y(cppz, cppy, dm%dcpp)
      call transpose_y_to_x(cppy, cppx, dm%dcpp)
      call Get_x_1st_derivative_C2P_3D(cppx, mx_pppx, dm, dm%ibcx(:,6),        &
                                       dm%fbcx_vof(:,:,:))

      if(dm%ibcx(1,6)==IBC_PERIODIC)then
        ibcx = IBC_PERIODIC
      else
        ibcx = IBC_DIRICHLET
      end if
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcy = IBC_PERIODIC
      else
        ibcy = IBC_DIRICHLET
      end if
      if(dm%ibcz(1,6)==IBC_PERIODIC)then
        ibcz = IBC_PERIODIC
      else
        ibcz = IBC_DIRICHLET
      end if
      call Get_x_midp_P2C_3D(mx_pppx, cppx, dm, ibcx)
      call transpose_x_to_y(cppx, cppy, dm%dcpp)
      call Get_y_midp_P2C_3D(cppy, ccpy, dm, ibcy)
      call transpose_y_to_z(ccpy, ccpz, dm%dccp)
      call Get_z_midp_P2C_3D(ccpz, cccz, dm, ibcz)
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, mx_cccx, dm%dccc)
    end if

    !calculate my on y-pencil
    if(vf%igrad==1) then
      call transpose_x_to_y(vf%phi, cccy, dm%dccc)
      call Get_y_midp_C2P_3D(cccy, cpcy, dm, dm%ibcy(:,6),                     &
                             dm%fbcy_vof(:,:,:))
      call Get_y_1st_derivative_P2C_3D(cpcy, cccy, dm, dm%ibcy(:,6),           &
                                       dm%fbcy_vof(:,:,:))
      call transpose_y_to_x(cccy, my_cccx, dm%dccc)
    else if (vf%igrad==2) then
      call Get_x_midp_C2P_3D(vf%phi, pccx, dm, dm%ibcx(:,6),                   &
                             dm%fbcx_vof(:,:,:))
      call transpose_x_to_y(pccx, pccy, dm%dpcc)
      call transpose_y_to_z(pccy, pccz, dm%dpcc)
      call Get_z_midp_C2P_3D(pccz, pcpz, dm, dm%ibcz(:,6),                     &
                             dm%fbcz_vof(:,:,:))
      call transpose_z_to_y(pcpz, pcpy, dm%dpcp)
      call Get_y_1st_derivative_C2P_3D(pcpy, pppy, dm, dm%ibcy(:,6),           &
                                       dm%fbcy_vof(:,:,:))
      call transpose_y_to_x(pppy, my_pppx, dm%dppp)

      call Get_y_midp_P2C_3D(pppy, pcpy, dm, ibcy)
      call transpose_y_to_z(pcpy, pcpz, dm%dpcp)
      call Get_z_midp_P2C_3D(pcpz, pccz, dm, ibcz)
      call transpose_z_to_y(pccz, pccy, dm%dpcc)
      call transpose_y_to_x(pccy, pccx, dm%dpcc)
      call Get_x_midp_P2C_3D(pccx, my_cccx, dm, ibcx)
    end if

    !calculate mz on z-pencil
    if(vf%igrad==1) then
      call transpose_x_to_y(vf%phi, cccy, dm%dccc)
      call transpose_y_to_z(cccy, cccz, dm%dccc)
      call Get_z_midp_C2P_3D(cccz, ccpz, dm, dm%ibcz(:,6),                     &
                             dm%fbcz_vof(:,:,:))
      call Get_z_1st_derivative_P2C_3D(ccpz, cccz, dm, dm%ibcz(:,6),           &
                                       dm%fbcz_vof(:,:,:))
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, mz_cccx, dm%dccc)
    else if (vf%igrad==2) then
      call Get_x_midp_C2P_3D(vf%phi, pccx, dm, dm%ibcx(:,6),                   &
                             dm%fbcx_vof(:,:,:))
      call transpose_x_to_y(pccx, pccy, dm%dpcc)
      call Get_y_midp_C2P_3D(pccy, ppcy, dm, dm%ibcy(:,6),                     &
                             dm%fbcy_vof(:,:,:))
      call transpose_y_to_z(ppcy, ppcz, dm%dppc)
      call Get_z_1st_derivative_C2P_3D(ppcz, pppz, dm, dm%ibcz(:,6),           &
                                       dm%fbcz_vof(:,:,:))
      call transpose_z_to_y(pppz, pppy, dm%dppp)
      call transpose_y_to_x(pppy, mz_pppx, dm%dppp)

      call Get_z_midp_P2C_3D(pppz, ppcz, dm, ibcz)
      call transpose_z_to_y(ppcz, ppcy, dm%dppc)
      call Get_y_midp_P2C_3D(ppcy, pccy, dm, ibcy)
      call transpose_y_to_x(pccy, pccx, dm%dpcc)
      call Get_x_midp_P2C_3D(pccx, mz_cccx, dm, ibcx)
    end if

    !calculate normal vector (nx, ny, nz) at cell centres, based on the global 
    !coordinate
    modcccx = sqrt(mx_cccx**2+my_cccx**2+mz_cccx**2)
    vf%lnx = mx_cccx/(modcccx+eps)
    vf%lny = my_cccx/(modcccx+eps)
    vf%lnz = mz_cccx/(modcccx+eps)

    !calculate normal vector (nx, ny, nz) at cell nodes, based on the global 
    !coordinate
    if(vf%igrad==2) then
      modpppx = sqrt(mx_pppx**2+my_pppx**2+mz_pppx**2)
      nx_pppx = mx_pppx/(modpppx+eps)
      ny_pppx = my_pppx/(modpppx+eps)
      nz_pppx = mz_pppx/(modpppx+eps)
    end if

    !calculate the Cartesian curvature tensor (lXX, lYY, lZZ, lXY, lYZ, lXZ)
    !component lXX
    if(vf%igrad==1) then
      if(dm%ibcx(1,6)==IBC_PERIODIC)then
        ibcx = dm%ibcx(:,6)
        fbcx = dm%fbcx_vof(:,:,:)
      else
        ibcx = IBC_INTRPL
        fbcx = dm%fbcx_vof(:,:,:)
      end if
      call Get_x_midp_C2P_3D(vf%lnx, pccx, dm, ibcx, fbcx)
      if(dm%ibcx(1,6)==IBC_PERIODIC)then
        ibcx = dm%ibcx(:,6)
        fbcx = dm%fbcx_vof(:,:,:)
      else
        ibcx = IBC_DIRICHLET
        fbcx(1,1:dm%nc(2),1:dm%nc(3)) = pccx(1,:,:)
        fbcx(2,1:dm%nc(2),1:dm%nc(3)) = pccx(dm%np(1),:,:)
        fbcx(3,:,:) = fbcx(1,:,:)
        fbcx(4,:,:) = fbcx(2,:,:)
      end if
      call Get_x_1st_derivative_P2C_3D(pccx, cccx, dm, ibcx, fbcx)
    else if(vf%igrad==2) then
      call transpose_x_to_y(nx_pppx, pppy, dm%dppp)
      call Get_y_midp_P2C_3D(pppy, pcpy, dm, ibcy)
      call transpose_y_to_z(pcpy, pcpz, dm%dpcp)
      call Get_z_midp_P2C_3D(pcpz, pccz, dm, ibcz)
      call transpose_z_to_y(pccz, pccy, dm%dpcc)
      call transpose_y_to_x(pccy, pccx, dm%dpcc)
      call Get_x_1st_derivative_P2C_3D(pccx, cccx, dm, ibcx)
    end if

    vf%llxx = cccx*(dm%h(1)**2)/hh

    !component lYY
    if(vf%igrad==1) then
      call transpose_x_to_y(vf%lny, cccy, dm%dccc)
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcy = dm%ibcy(:,6)
        fbcy = dm%fbcy_vof(:,:,:)
      else
        ibcy = IBC_INTRPL
        fbcy = dm%fbcy_vof(:,:,:)
      end if
      call Get_y_midp_C2P_3D(cccy, cpcy, dm, ibcy, fbcy)
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcy = dm%ibcy(:,6)
        fbcy = dm%fbcy_vof(:,:,:)
      else
        ibcy = IBC_DIRICHLET
        fbcy(1:dm%nc(1),1,1:dm%nc(3)) = cpcy(:,1,:)
        fbcy(1:dm%nc(1),2,1:dm%nc(3)) = cpcy(:,dm%np(2),:)
        fbcy(:,3,:) = fbcy(:,1,:)
        fbcy(:,4,:) = fbcy(:,2,:)
      end if
      call Get_y_1st_derivative_P2C_3D(cpcy, cccy, dm, ibcy, fbcy)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
    else if(vf%igrad==2) then
      call Get_x_midp_P2C_3D(ny_pppx, cppx, dm, ibcx)
      call transpose_x_to_y(cppx, cppy, dm%dcpp)
      call transpose_y_to_z(cppy, cppz, dm%dcpp)
      call Get_z_midp_P2C_3D(cppz, cpcz, dm, ibcz)
      call transpose_z_to_y(cpcz, cpcy, dm%dcpc)
      call Get_y_1st_derivative_P2C_3D(cpcy, cccy, dm, ibcy)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
    end if

    vf%llyy = cccx*(dm%h(2)**2)/hh

    !component lZZ
    if(vf%igrad==1) then
      call transpose_x_to_y(vf%lnz, cccy, dm%dccc)
      call transpose_y_to_z(cccy, cccz, dm%dccc)
      if(dm%ibcz(1,6)==IBC_PERIODIC)then
        ibcz = dm%ibcz(:,6)
        fbcz = dm%fbcz_vof(:,:,:)
      else
        ibcz = IBC_INTRPL
        fbcz = dm%fbcz_vof(:,:,:)
      end if
      call Get_z_midp_C2P_3D(cccz, ccpz, dm, ibcz, fbcz)
      if(dm%ibcz(1,6)==IBC_PERIODIC)then
        ibcz = dm%ibcz(:,6)
        fbcz = dm%fbcz_vof(:,:,:)
      else
        ibcz = IBC_DIRICHLET
        fbcz(1:dm%nc(1),1:dm%nc(2),1) = ccpz(:,:,1)
        fbcz(1:dm%nc(1),1:dm%nc(2),2) = ccpz(:,:,dm%np(3))
        fbcz(:,:,3) = fbcz(:,:,1)
        fbcz(:,:,4) = fbcz(:,:,2)
      end if
      call Get_z_1st_derivative_P2C_3D(ccpz, cccz, dm, ibcz, fbcz)
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
    else if(vf%igrad==2) then
      call Get_x_midp_P2C_3D(nz_pppx, cppx, dm, ibcx)
      call transpose_x_to_y(cppx, cppy, dm%dcpp)
      call Get_y_midp_P2C_3D(cppy, ccpy, dm, ibcy)
      call transpose_y_to_z(ccpy, ccpz, dm%dccp)
      call Get_z_1st_derivative_P2C_3D(ccpz, cccz, dm, ibcz)
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
    end if

    vf%llzz = cccx*(dm%h(3)**2)/hh

    !component lXY
    if(vf%igrad==1) then
      if(dm%ibcx(1,6)==IBC_PERIODIC)then
        ibcx = dm%ibcx(:,6)
        fbcx = dm%fbcx_vof(:,:,:)
      else
        ibcx = IBC_INTRPL
        fbcx = dm%fbcx_vof(:,:,:)
      end if
      call Get_x_midp_C2P_3D(vf%lny, pccx, dm, ibcx, fbcx)
      if(dm%ibcx(1,6)==IBC_PERIODIC)then
        ibcx = dm%ibcx(:,6)
        fbcx = dm%fbcx_vof(:,:,:)
      else
        ibcx = IBC_DIRICHLET
        fbcx(1,1:dm%nc(1),1:dm%nc(3)) = pccx(1,:,:)
        fbcx(2,1:dm%nc(1),1:dm%nc(3)) = pccx(dm%np(1),:,:)
        fbcx(3,:,:) = fbcx(1,:,:)
        fbcx(4,:,:) = fbcx(2,:,:)
      end if
      call Get_x_1st_derivative_P2C_3D(pccx, cccx, dm, ibcx, fbcx)
      vf%llxy = cccx*dm%h(1)*dm%h(2)/hh*HALF
      call transpose_x_to_y(vf%lnx, cccy, dm%dccc)
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcy = dm%ibcy(:,6)
        fbcy = dm%fbcy_vof(:,:,:)
      else
        ibcy = IBC_INTRPL
        fbcy = dm%fbcy_vof(:,:,:)
      end if
      call Get_y_midp_C2P_3D(cccy, cpcy, dm, ibcy, fbcy)
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcy = dm%ibcy(:,6)
        fbcy = dm%fbcy_vof(:,:,:)
      else
        ibcy = IBC_DIRICHLET
        fbcy(1:dm%nc(1),1,1:dm%nc(3)) = cpcy(:,1,:)
        fbcy(1:dm%nc(1),2,1:dm%nc(3)) = cpcy(:,dm%np(2),:)
        fbcy(:,3,:) = fbcy(:,1,:)
        fbcy(:,4,:) = fbcy(:,2,:)
      end if
      call Get_y_1st_derivative_P2C_3D(cpcy, cccy, dm, ibcy, fbcy)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llxy = vf%llxy+cccx*dm%h(1)*dm%h(2)/hh*HALF
    else if(vf%igrad==2) then
      call Get_x_midp_P2C_3D(nx_pppx, cppx, dm, ibcx)
      call transpose_x_to_y(cppx, cppy, dm%dcpp)
      call transpose_y_to_z(cppy, cppz, dm%dcpp)
      call Get_z_midp_P2C_3D(cppz, cpcz, dm, ibcz)
      call transpose_z_to_y(cpcz, cpcy, dm%dcpc)
      call Get_y_1st_derivative_P2C_3D(cpcy, cccy, dm, ibcy)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llxy = cccx*dm%h(1)*dm%h(2)/hh*HALF
      call transpose_x_to_y(ny_pppx, pppy, dm%dppp)
      call Get_y_midp_P2C_3D(pppy, pcpy, dm, ibcy)
      call transpose_y_to_z(pcpy, pcpz, dm%dpcp)
      call Get_z_midp_P2C_3D(pcpz, pccz, dm, ibcz)
      call transpose_z_to_y(pccz, pccy, dm%dpcc)
      call transpose_y_to_x(pccy, pccx, dm%dpcc)
      call Get_x_1st_derivative_P2C_3D(pccx, cccx, dm, ibcx)
      vf%llxy = vf%llxy+cccx*dm%h(1)*dm%h(2)/hh*HALF
    end if

    !component lXZ
    if(vf%igrad==1) then
      if(dm%ibcx(1,6)==IBC_PERIODIC)then
        ibcx = dm%ibcx(:,6)
        fbcx = dm%fbcx_vof(:,:,:)
      else
        ibcx = IBC_INTRPL
        fbcx = dm%fbcx_vof(:,:,:)
      end if
      call Get_x_midp_C2P_3D(vf%lnz, pccx, dm, ibcx, fbcx)
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcx = dm%ibcx(:,6)
        fbcx = dm%fbcx_vof(:,:,:)
      else
        ibcx = IBC_DIRICHLET
        fbcx(1,1:dm%nc(2),1:dm%nc(3)) = pccx(1,:,:)
        fbcx(2,1:dm%nc(2),1:dm%nc(3)) = pccx(dm%np(1),:,:)
        fbcx(3,:,:) = fbcx(1,:,:)
        fbcx(4,:,:) = fbcx(2,:,:)
      end if
      call Get_x_1st_derivative_P2C_3D(pccx, cccx, dm, ibcx, fbcx)
      vf%llxz = cccx*dm%h(1)*dm%h(3)/hh*HALF
      call transpose_x_to_y(vf%lnx, cccy, dm%dccc)
      call transpose_y_to_z(cccy, cccz, dm%dccc)
      if(dm%ibcz(1,6)==IBC_PERIODIC)then
        ibcz = dm%ibcz(:,6)
        fbcz = dm%fbcz_vof(:,:,:)
      else
        ibcz = IBC_INTRPL
        fbcz = dm%fbcz_vof(:,:,:)
      end if
      call Get_z_midp_C2P_3D(cccz, ccpz, dm, ibcz, fbcz)
      if(dm%ibcz(1,6)==IBC_PERIODIC)then
        ibcz = dm%ibcz(:,6)
        fbcz = dm%fbcz_vof(:,:,:)
      else
        ibcz = IBC_DIRICHLET
        fbcz(1:dm%nc(1),1:dm%nc(2),1) = ccpz(:,:,1)
        fbcz(1:dm%nc(1),1:dm%nc(2),2) = ccpz(:,:,dm%np(3))
        fbcz(:,:,3) = fbcz(:,:,1)
        fbcz(:,:,4) = fbcz(:,:,2)
      end if
      call Get_z_1st_derivative_P2C_3D(ccpz, cccz, dm, ibcz, fbcz)
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llxz = vf%llxz+cccx*dm%h(1)*dm%h(3)/hh*HALF
    else if(vf%igrad==2) then
      call Get_x_midp_P2C_3D(nx_pppx, cppx, dm, ibcx)
      call transpose_x_to_y(cppx, cppy, dm%dcpp)
      call Get_y_midp_P2C_3D(cppy, ccpy, dm, ibcy)
      call transpose_y_to_z(ccpy, ccpz, dm%dccp)
      call Get_z_1st_derivative_P2C_3D(ccpz, cccz, dm, ibcz)
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llxz = cccx*dm%h(1)*dm%h(3)/hh*HALF
      call transpose_x_to_y(nz_pppx, pppy, dm%dppp)
      call Get_y_midp_P2C_3D(pppy, pcpy, dm, ibcy)
      call transpose_y_to_z(pcpy, pcpz, dm%dpcp)
      call Get_z_midp_P2C_3D(pcpz, pccz, dm, ibcz)
      call transpose_z_to_y(pccz, pccy, dm%dpcc)
      call transpose_y_to_x(pccy, pccx, dm%dpcc)
      call Get_x_1st_derivative_P2C_3D(pccx, cccx, dm, ibcx)
      vf%llxz = vf%llxz+cccx*dm%h(1)*dm%h(3)/hh*HALF
    end if

    !component lYZ
    if(vf%igrad==1) then
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcy = dm%ibcy(:,6)
        fbcy = dm%fbcy_vof(:,:,:)
      else
        ibcy = IBC_INTRPL
        fbcy = dm%fbcy_vof(:,:,:)
      end if
      call transpose_x_to_y(vf%lnz, cccy, dm%dccc)
      call Get_y_midp_C2P_3D(cccy, cpcy, dm, ibcy, fbcy)
      if(dm%ibcy(1,6)==IBC_PERIODIC)then
        ibcy = dm%ibcy(:,6)
        fbcy = dm%fbcy_vof(:,:,:)
      else
        ibcy = IBC_DIRICHLET
        fbcy(1:dm%nc(1),1,1:dm%nc(3)) = cpcy(:,1,:)
        fbcy(1:dm%nc(1),2,1:dm%nc(3)) = cpcy(:,dm%np(2),:)
        fbcy(:,3,:) = fbcy(:,1,:)
        fbcy(:,4,:) = fbcy(:,2,:)
      end if
      call Get_y_1st_derivative_P2C_3D(cpcy, cccy, dm, ibcy, fbcy)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llyz = cccx*dm%h(2)*dm%h(3)/hh*HALF
      call transpose_x_to_y(vf%lny, cccy, dm%dccc)
      call transpose_y_to_z(cccy, cccz, dm%dccc)
      if(dm%ibcz(1,6)==IBC_PERIODIC)then
        ibcz = dm%ibcz(:,6)
        fbcz = dm%fbcz_vof(:,:,:)
      else
        ibcz = IBC_INTRPL
        fbcz = dm%fbcz_vof(:,:,:)
      end if
      call Get_z_midp_C2P_3D(cccz, ccpz, dm, ibcz, fbcz)
      if(dm%ibcz(1,6)==IBC_PERIODIC)then
        ibcz = dm%ibcz(:,6)
        fbcz = dm%fbcz_vof(:,:,:)
      else
        ibcz = IBC_DIRICHLET
        fbcz(1:dm%nc(1),1:dm%nc(2),1) = ccpz(:,:,1)
        fbcz(1:dm%nc(1),1:dm%nc(2),2) = ccpz(:,:,dm%np(3))
        fbcz(:,:,3) = fbcz(:,:,1)
        fbcz(:,:,4) = fbcz(:,:,2)
      end if
      call Get_z_1st_derivative_P2C_3D(ccpz, cccz, dm, ibcz, fbcz)
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llyz = vf%llyz+cccx*dm%h(2)*dm%h(3)/hh*HALF
    else if(vf%igrad==2) then
      call Get_x_midp_P2C_3D(ny_pppx, cppx, dm, ibcx)
      call transpose_x_to_y(cppx, cppy, dm%dcpp)
      call Get_y_midp_P2C_3D(cppy, ccpy, dm, ibcy)
      call transpose_y_to_z(ccpy, ccpz, dm%dccp)
      call Get_z_1st_derivative_P2C_3D(ccpz, cccz, dm, ibcz)
      call transpose_z_to_y(cccz, cccy, dm%dccc)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llyz = cccx*dm%h(2)*dm%h(3)/hh*HALF
      call Get_x_midp_P2C_3D(nz_pppx, cppx, dm, ibcx)
      call transpose_x_to_y(cppx, cppy, dm%dcpp)
      call transpose_y_to_z(cppy, cppz, dm%dcpp)
      call Get_z_midp_P2C_3D(cppz, cpcz, dm, ibcz)
      call transpose_z_to_y(cpcz, cpcy, dm%dcpc)
      call Get_y_1st_derivative_P2C_3D(cpcy, cccy, dm, ibcy)
      call transpose_y_to_x(cccy, cccx, dm%dccc)
      vf%llyz = vf%llyz+cccx*dm%h(2)*dm%h(3)/hh*HALF
    end if

    !convert normal vector (nx, ny, nz) to local-coordinate (nX, nY, nZ)
    vf%lnx = vf%lnx*dm%h(1)/hh
    vf%lny = vf%lny*dm%h(2)/hh
    vf%lnz = vf%lnz*dm%h(3)/hh

    !skip calculaton of P(x) when vof is too close to 0 or 1
    limit = vf%voflim
    where (vf%phi<limit .or. vf%phi>(ONE-limit))
      vf%lnx = ZERO
      vf%lny = ZERO
      vf%lnz = ZERO
      vf%llxx = ZERO
      vf%llyy = ZERO
      vf%llzz = ZERO
      vf%llxy = ZERO
      vf%llxz = ZERO
      vf%llyz = ZERO
    end where

    ! in case of linear reconstrcution
    if(vf%ireconstruct==1) then
      vf%llxx = ZERO
      vf%llyy = ZERO
      vf%llzz = ZERO
      vf%llxy = ZERO
      vf%llxz = ZERO
      vf%llyz = ZERO
    end if

    !calculate face curvature kappa
    vf%kappa = -hh*(vf%llxx/dm%h(1)**2+vf%llyy/dm%h(2)**2+vf%llzz/dm%h(3)**2)

    !calculate d
    y1 = (ONE-ONE/sqrt(THREE))/TWO
    y2 = (ONE+ONE/sqrt(THREE))/TWO
    do k=1, dm%dccc%xsz(3)
      do j=1, dm%dccc%xsz(2)
        do i=1, dm%dccc%xsz(1)

          lnx = vf%lnx(i,j,k)
          lny = vf%lny(i,j,k)
          lnz = vf%lnz(i,j,k)
          phi = vf%phi(i,j,k)
          llxx = vf%llxx(i,j,k)
          llyy = vf%llyy(i,j,k)
          llzz = vf%llzz(i,j,k)
          llxy = vf%llxy(i,j,k)
          llxz = vf%llxz(i,j,k)
          llyz = vf%llyz(i,j,k)

          if(phi<limit .or. phi>(ONE-limit)) then
            vf%dd(i,j,k) = -1000._WP
            cycle
          end if
          if(abs(lnx)<eps .and. abs(lny)<eps .and. abs(lnz)<eps) then
            vf%dd(i,j,k) = -1000._WP
            cycle
          end if

          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz))) then
            cx = ZERO
            cy = ONE
            cz = ONE
          else if (abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
            cx = ONE
            cy = ZERO
            cz = ONE
          else
            cx = ONE
            cy = ONE
            cz = ZERO
          end if

          vf%a100(i,j,k) = lnx-HALF*(cx*llxx+cx*cy*llxy+cx*cz*llxz)
          vf%a010(i,j,k) = lny-HALF*(cy*llyy+cx*cy*llxy+cy*cz*llyz)
          vf%a001(i,j,k) = lnz-HALF*(cz*llzz+cx*cz*llxz+cy*cz*llyz)
          vf%a200(i,j,k) = HALF*cx*llxx
          vf%a020(i,j,k) = HALF*cy*llyy
          vf%a002(i,j,k) = HALF*cz*llzz
          vf%a110(i,j,k) = cx*cy*llxy
          vf%a101(i,j,k) = cx*cz*llxz
          vf%a011(i,j,k) = cy*cz*llyz

          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz))) then
            a100 = vf%a100(i,j,k)
            a010 = vf%a010(i,j,k)
            a001 = vf%a001(i,j,k)
            a020 = vf%a020(i,j,k)
            a002 = vf%a002(i,j,k)
            a011 = vf%a011(i,j,k)
            A  = exp(TWO*beta*a100)
            B1 = exp(TWO*beta*(a020*(y1**2)+a002*(y1**2)+a011*(y1**2)+         &
                               a010*y1+a001*y1))
            B2 = exp(TWO*beta*(a020*(y1**2)+a002*(y2**2)+a011*(y1*y2)+         &
                               a010*y1+a001*y2))
            C1 = exp(TWO*beta*(a020*(y2**2)+a002*(y1**2)+a011*(y1*y2)+         &
                               a010*y2+a001*y1))
            C2 = exp(TWO*beta*(a020*(y2**2)+a002*(y2**2)+a011*(y2**2)+         &
                               a010*y2+a001*y2))
            Q  = exp(FOUR*beta*a100*(TWO*phi-ONE))

            !for 2D (x-y)
            A2d  = exp(TWO*beta*a100)
            B12d = exp(TWO*beta*(a020*(y1**2)+a010*y1))
            B22d = exp(TWO*beta*(a020*(y2**2)+a010*y2))
            Q2d  = exp(TWO*beta*a100*(TWO*phi-ONE))

          else if (abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
            a100 = vf%a100(i,j,k)
            a010 = vf%a010(i,j,k)
            a001 = vf%a001(i,j,k)
            a200 = vf%a200(i,j,k)
            a002 = vf%a002(i,j,k)
            a101 = vf%a101(i,j,k)
            A  = exp(TWO*beta*a010)
            B1 = exp(TWO*beta*(a200*(y1**2)+a002*(y1**2)+a101*(y1**2)+         &
                               a100*y1+a001*y1))
            B2 = exp(TWO*beta*(a200*(y1**2)+a002*(y2**2)+a101*(y1*y2)+         &
                               a100*y1+a001*y2))
            C1 = exp(TWO*beta*(a200*(y2**2)+a002*(y1**2)+a101*(y1*y2)+         &
                               a100*y2+a001*y1))
            C2 = exp(TWO*beta*(a200*(y2**2)+a002*(y2**2)+a101*(y2**2)+         &
                               a100*y2+a001*y2))
            Q  = exp(FOUR*beta*a010*(TWO*phi-ONE))

            !for 2D (x-y)
            A2d  = exp(TWO*beta*a010)
            B12d = exp(TWO*beta*(a200*(y1**2)+a100*y1))
            B22d = exp(TWO*beta*(a200*(y2**2)+a100*y2))
            Q2d  = exp(TWO*beta*a010*(TWO*phi-ONE))

          else
            a100 = vf%a100(i,j,k)
            a010 = vf%a010(i,j,k)
            a001 = vf%a001(i,j,k)
            a200 = vf%a200(i,j,k)
            a020 = vf%a020(i,j,k)
            a110 = vf%a110(i,j,k)
            A  = exp(TWO*beta*a001)
            B1 = exp(TWO*beta*(a200*(y1**2)+a020*(y1**2)+a110*(y1**2)+         &
                               a100*y1+a010*y1))
            B2 = exp(TWO*beta*(a200*(y1**2)+a020*(y2**2)+a110*(y1*y2)+         &
                               a100*y1+a010*y2))
            C1 = exp(TWO*beta*(a200*(y2**2)+a020*(y1**2)+a110*(y1*y2)+         &
                               a100*y2+a010*y1))
            C2 = exp(TWO*beta*(a200*(y2**2)+a020*(y2**2)+a110*(y2**2)+         &
                               a100*y2+a010*y2))
            Q  = exp(FOUR*beta*a001*(TWO*phi-ONE))
          end if

          !for 2D (x-y)
          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz)) .or.                   &
              abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
            aa2d = A2d*B12d*B22d*(A2d-Q2d)
            bb2d = A2d*(B12d+B22d)*(ONE-Q2d)
            cc2d = ONE-A2d*Q2d

            solu1 = (-bb2d+sqrt(bb2d**2-FOUR*aa2d*cc2d))/(TWO*aa2d+eps)
            solu2 = (-bb2d-sqrt(bb2d**2-FOUR*aa2d*cc2d))/(TWO*aa2d+eps)

            if (solu1>ZERO)then
              D2d = solu1
            else if (solu2>ZERO)then
              D2d = solu2
            else
              write(*,*)'Error in calculating d at cell: ',i,j,k
              D2d = ONE
            end if
          end if

          aa0 = ONE-(A**2)*Q
          aa1 = A*(B1+B2+C1+C2)*(ONE-A*Q)
          aa2 = (A**2)*((B1+B2)*(C1+C2)+B1*B2+C1*C2)*(ONE-Q)
          aa3 = (A**2)*(B1*B2*(C1+C2)+(B1+B2)*C1*C2)*(A-Q)
          aa4 = (A**2)*B1*B2*C1*C2*(A**2-Q)

          D = quartic(aa0, aa1, aa2, aa3, aa4)
          vf%dd(i,j,k) = log(D)/(TWO*beta)
        enddo
      enddo
    enddo

    contains

    real(WP) function cbrt(var)

      implicit none

      real(WP), intent(IN) :: var

      if(var>=ZERO) then
        cbrt = var**(ONE/THREE)
      else
        cbrt = -(-var)**(ONE/THREE)
      end if

    end function cbrt

    real(WP) function quartic(aa0, aa1, aa2, aa3, aa4)

      implicit none

      real(WP), intent(IN) :: aa0, aa1, aa2, aa3, aa4
      real(WP) :: bb0, bb1, bb2, bb3, cc0, cc1, cc2
      real(WP) :: aa, bb, z1, delta
      real(WP) :: ae, be, ce, de2, de

      bb0 = aa0/aa4
      bb1 = aa1/aa4
      bb2 = aa2/aa4
      bb3 = aa3/aa4

      cc0 = bb0*(FOUR*bb2-bb3**2)-bb1**2
      cc1 = bb1*bb3-FOUR*bb0
      cc2 = -bb2

      aa = -cc2**2/NINE+cc1/THREE
      bb = TWO*(cc2**3)/(THREE*NINE)-cc1*cc2/THREE+cc0
      delta = bb**2+FOUR*(aa**3)

      if (delta >= ZERO)then
        z1 = cbrt(HALF*(-bb+sqrt(delta)))-cbrt(HALF*(bb+sqrt(delta)))-cc2/THREE
      else
        z1 = TWO*sqrt(-aa)*cos(atan(-sqrt(-delta)/bb)/THREE)-cc2/THREE
      end if

      ae = HALF*bb3
      be = HALF*z1
      de2 = be**2-bb0
      if(de2<=ZERO) then
        de = limit
      else
        de = sqrt(de2)
      end if
      ce = (-HALF*bb1+ae*be)/de

      quartic = HALF*(-(ae-ce)+sqrt((ae-ce)**2-FOUR*(be-de)))

    end function

  end subroutine

!===============================================================================
!===============================================================================
  !calculate the numerical fluxes
  real(WP) function pxyz(x, y, z,                                              &
                         a100, a010, a001,                                     &
                         a200, a020, a002,                                     &
                         a110, a011, a101)

    implicit none

    real(WP), intent(IN) :: x, y, z
    real(WP), intent(IN) :: a100, a010, a001, a200, a020, a002, a110, a011, a101

    pxyz = a200*(x**2)+a020*(y**2)+a002*(z**2)+                                &
           a110*x*y+a011*y*z+a101*x*z+                                         &
           a100*x+a010*y+a001*z

    end function pxyz

!===============================================================================
  real(WP) function itg_h_hat(var1, var2,                                      &
                              a100, a010, a001,                                &
                              a200, a020, a002,                                &
                              a110, a011, a101,                                &
                              d, beta,                                         &
                              xitg_a, xitg_b,                                  &
                              yitg_a, yitg_b,                                  &
                              zitg_a, zitg_b,                                  &
                              dirct)

    use parameters_constant_mod

    implicit none

    real(WP), intent(IN) :: var1, var2

    real(WP), intent(IN) :: a100, a010, a001, a200, a020, a002, a110, a011, a101
    real(WP), intent(IN) :: d, beta
    real(WP), intent(IN) :: xitg_a, xitg_b, yitg_a, yitg_b, zitg_a, zitg_b
    integer , intent(IN) :: dirct

    real(WP) :: term1, dpdvar, term2coef, term2a, term2b
    real(WP), parameter :: eps = 1.e-16_WP

    if (dirct==1) then
      term1 = xitg_b-xitg_a
      dpdvar = a100
      term2coef = ONE/(beta*dpdvar+eps)
      term2a = cosh(beta*(pxyz(xitg_a, var1, var2,                             &
                               a100, a010, a001,                               &
                               a200, a020, a002,                               &
                               a110, a011, a101)+d))
      term2b = cosh(beta*(pxyz(xitg_b, var1, var2,                             &
                               a100, a010, a001,                               &
                               a200, a020, a002,                               &
                               a110, a011, a101)+d))
    else if (dirct==2) then
      term1 = yitg_b-yitg_a
      dpdvar = a010
      term2coef = ONE/(beta*dpdvar+eps)
      term2a = cosh(beta*(pxyz(var1, yitg_a, var2,                             &
                               a100, a010, a001,                               &
                               a200, a020, a002,                               &
                               a110, a011, a101)+d))
      term2b = cosh(beta*(pxyz(var1, yitg_b, var2,                             &
                               a100, a010, a001,                               &
                               a200, a020, a002,                               &
                               a110, a011, a101)+d))
    else
      term1 = zitg_b-zitg_a
      dpdvar = a001
      term2coef = ONE/(beta*dpdvar+eps)
      term2a = cosh(beta*(pxyz(var1, var2, zitg_a,                             &
                               a100, a010, a001,                               &
                               a200, a020, a002,                               &
                               a110, a011, a101)+d))
      term2b = cosh(beta*(pxyz(var1, var2, zitg_b,                             &
                               a100, a010, a001,                               &
                               a200, a020, a002,                               &
                               a110, a011, a101)+d))
    end if

    itg_h_hat = HALF*(term1+term2coef*log(term2b/term2a))

  end function itg_h_hat

!===============================================================================
  real(WP) function cal_flux(a100, a010, a001,                                 &
                             a200, a020, a002,                                 &
                             a110, a011, a101,                                 &
                             d, beta,                                          &
                             xitg_a, xitg_b,                                   &
                             yitg_a, yitg_b,                                   &
                             zitg_a, zitg_b,                                   &
                             dirct)

    use parameters_constant_mod

    implicit none

    real(WP), intent(IN) :: a100, a010, a001, a200, a020, a002, a110, a011, a101
    real(WP), intent(IN) :: d, beta
    real(WP), intent(IN) :: xitg_a, xitg_b, yitg_a, yitg_b, zitg_a, zitg_b
    integer , intent(IN) :: dirct

    real(WP) :: var11, var21, var12, var22
    real(WP) :: delta1, delta2
    real(WP) :: itg11, itg21, itg12, itg22

    if (dirct==1) then
      var11 = HALF*(yitg_b-yitg_a)/sqrt(THREE)+HALF*(yitg_a+yitg_b)
      var21 = HALF*(yitg_a-yitg_b)/sqrt(THREE)+HALF*(yitg_a+yitg_b)
      delta1 = yitg_b-yitg_a
      var12 = HALF*(zitg_b-zitg_a)/sqrt(THREE)+HALF*(zitg_a+zitg_b)
      var22 = HALF*(zitg_a-zitg_b)/sqrt(THREE)+HALF*(zitg_a+zitg_b)
      delta2 = zitg_b-zitg_a
    else if (dirct==2) then
      var11 = HALF*(xitg_b-xitg_a)/sqrt(THREE)+HALF*(xitg_a+xitg_b)
      var21 = HALF*(xitg_a-xitg_b)/sqrt(THREE)+HALF*(xitg_a+xitg_b)
      delta1 = xitg_b-xitg_a
      var12 = HALF*(zitg_b-zitg_a)/sqrt(THREE)+HALF*(zitg_a+zitg_b)
      var22 = HALF*(zitg_a-zitg_b)/sqrt(THREE)+HALF*(zitg_a+zitg_b)
      delta2 = zitg_b-zitg_a
    else
      var11 = HALF*(xitg_b-xitg_a)/sqrt(THREE)+HALF*(xitg_a+xitg_b)
      var21 = HALF*(xitg_a-xitg_b)/sqrt(THREE)+HALF*(xitg_a+xitg_b)
      delta1 = xitg_b-xitg_a
      var12 = HALF*(yitg_b-yitg_a)/sqrt(THREE)+HALF*(yitg_a+yitg_b)
      var22 = HALF*(yitg_a-yitg_b)/sqrt(THREE)+HALF*(yitg_a+yitg_b)
      delta2 = yitg_b-yitg_a
    end if

    itg11 = itg_h_hat(var11, var12,                                            &
                      a100, a010, a001,                                        &
                      a200, a020, a002,                                        &
                      a110, a011, a101,                                        &
                      d, beta,                                                 &
                      xitg_a, xitg_b,                                          &
                      yitg_a, yitg_b,                                          &
                      zitg_a, zitg_b,                                          &
                      dirct)

    itg12 = itg_h_hat(var11, var22,                                            &
                      a100, a010, a001,                                        &
                      a200, a020, a002,                                        &
                      a110, a011, a101,                                        &
                      d, beta,                                                 &
                      xitg_a, xitg_b,                                          &
                      yitg_a, yitg_b,                                          &
                      zitg_a, zitg_b,                                          &
                      dirct)

    itg21 = itg_h_hat(var21, var12,                                            &
                      a100, a010, a001,                                        &
                      a200, a020, a002,                                        &
                      a110, a011, a101,                                        &
                      d, beta,                                                 &
                      xitg_a, xitg_b,                                          &
                      yitg_a, yitg_b,                                          &
                      zitg_a, zitg_b,                                          &
                      dirct)

    itg22 = itg_h_hat(var21, var22,                                            &
                      a100, a010, a001,                                        &
                      a200, a020, a002,                                        &
                      a110, a011, a101,                                        &
                      d, beta,                                                 &
                      xitg_a, xitg_b,                                          &
                      yitg_a, yitg_b,                                          &
                      zitg_a, zitg_b,                                          &
                      dirct)

    cal_flux = QUARTER*delta1*delta2*(itg11+itg12+itg21+itg22)

  end function cal_flux

!===============================================================================
!===============================================================================
  !clip non-physical vof
  subroutine clip_vof(dm, vf)
    use parameters_constant_mod
    use operations
    use udf_type_mod

    implicit none

    type(t_domain), intent(IN)    :: dm
    type(t_vof),    intent(INOUT) :: vf

    integer :: i, j, k

    do k=1, dm%dccc%xsz(3)
      do j=1, dm%dccc%xsz(2)
        do i=1, dm%dccc%xsz(1)
          vf%phi(i,j,k) = min(max(ZERO, vf%phi(i,j,k)), ONE)
        end do
      end do
    end do

  end subroutine clip_vof

!===============================================================================
!===============================================================================
  !update the vof field
  subroutine Solve_vof_eq(dm, fl, vf)
    use udf_type_mod
    use solver_tools_mod
    use parameters_constant_mod
    use thermo_info_mod
    use vof_info_mod

    implicit none

    type(t_domain), intent(IN)    :: dm
    type(t_flow),   intent(INOUT) :: fl
    type(t_vof),    intent(INOUT) :: vf

    real(WP),dimension(dm%dccc%xsz(1),dm%dccc%xsz(2),dm%dccc%xsz(3)) ::        &
                                                                   rhsx,       &
                                                                   rhsy,       &
                                                                   rhsz
                                                                   
    real(WP),dimension(dm%dccc%ysz(1),dm%dccc%ysz(2),dm%dccc%ysz(3)) ::        &
                                                                   tmp_cccy,   &
                                                                   phi_cccy,   &
                                                                   phi2_cccy,  &
                                                                   lnx_cccy,   &
                                                                   lny_cccy,   &
                                                                   lnz_cccy,   &
                                                                   a100_cccy,  &
                                                                   a010_cccy,  &
                                                                   a001_cccy,  &
                                                                   a200_cccy,  &
                                                                   a020_cccy,  &
                                                                   a002_cccy,  &
                                                                   a110_cccy,  &
                                                                   a011_cccy,  &
                                                                   a101_cccy,  &
                                                                   dd_cccy

    real(WP),dimension(dm%dcpc%ysz(1),dm%dcpc%ysz(2),dm%dcpc%ysz(3)) ::        &
                                                                   qy_cpcy,    &
                                                                   gly_cpcy

    real(WP),dimension(dm%dccp%ysz(1),dm%dccp%ysz(2),dm%dccp%ysz(3)) ::        &
                                                                   qz_ccpy

    real(WP),dimension(dm%dccc%zsz(1),dm%dccc%zsz(2),dm%dccc%zsz(3)) ::        &
                                                                   tmp_cccz,   &
                                                                   phi_cccz,   &
                                                                   phi3_cccz,  &
                                                                   lnx_cccz,   &
                                                                   lny_cccz,   &
                                                                   lnz_cccz,   &
                                                                   a100_cccz,  &
                                                                   a010_cccz,  &
                                                                   a001_cccz,  &
                                                                   a200_cccz,  &
                                                                   a020_cccz,  &
                                                                   a002_cccz,  &
                                                                   a110_cccz,  &
                                                                   a011_cccz,  &
                                                                   a101_cccz,  &
                                                                   dd_cccz

    real(WP),dimension(dm%dccp%zsz(1),dm%dccp%zsz(2),dm%dccp%zsz(3)) ::        &
                                                                   qz_ccpz,    &
                                                                   hlz_ccpz


    integer :: i, j, k
    integer :: dirct
    real(WP) :: u, v, w
    real(WP) :: lnx, lny, lnz
    real(WP) :: xitg_a, xitg_b, yitg_a, yitg_b, zitg_a, zitg_b
    real(WP) :: phi, coef, d, beta
    real(WP) :: a200, a020, a002, a110, a011, a101, a100, a010, a001

    type(DECOMP_INFO) :: dtmp

    beta = real(vf%ibeta,WP)

    !calculate the x-direction numerical flux over the x-pencil (serial/parallel)
!    call cal_colour_function(dm, vf)

    dtmp = dm%dpcc
    do k=1, dtmp%xsz(3)
      do j=1, dtmp%xsz(2)
        do i=1, dtmp%xsz(1)
          u = fl%qx(i,j,k)
          lnx =  ZERO
          lny =  ZERO
          lnz =  ZERO
          a100 = ZERO
          a010 = ZERO
          a001 = ZERO
          a200 = ZERO
          a020 = ZERO
          a002 = ZERO
          a110 = ZERO
          a011 = ZERO
          a101 = ZEro
          d    = ZERO
          phi  = ZERO
          if (u>=ZERO) then
            coef = ONE
            if (i==1 .and. dm%ibcx(1,6)/=IBC_PERIODIC) then
              coef = ZERO
            else if (i==1 .and. dm%ibcx(1,6)==IBC_PERIODIC) then
              lnx  = vf%lnx(dtmp%xsz(1),j,k)
              lny  = vf%lny(dtmp%xsz(1),j,k)
              lnz  = vf%lnz(dtmp%xsz(1),j,k)
              a100 = vf%a100(dtmp%xsz(1),j,k)
              a010 = vf%a010(dtmp%xsz(1),j,k)
              a001 = vf%a001(dtmp%xsz(1),j,k)
              a200 = vf%a200(dtmp%xsz(1),j,k)
              a020 = vf%a020(dtmp%xsz(1),j,k)
              a002 = vf%a002(dtmp%xsz(1),j,k)
              a110 = vf%a110(dtmp%xsz(1),j,k)
              a011 = vf%a011(dtmp%xsz(1),j,k)
              a101 = vf%a101(dtmp%xsz(1),j,k)
              d    = vf%dd  (dtmp%xsz(1),j,k)
              phi  = vf%phi (dtmp%xsz(1),j,k)
            else
              lnx  = vf%lnx(i-1,j,k)
              lny  = vf%lny(i-1,j,k)
              lnz  = vf%lnz(i-1,j,k)
              a100 = vf%a100(i-1,j,k)
              a010 = vf%a010(i-1,j,k)
              a001 = vf%a001(i-1,j,k)
              a200 = vf%a200(i-1,j,k)
              a020 = vf%a020(i-1,j,k)
              a002 = vf%a002(i-1,j,k)
              a110 = vf%a110(i-1,j,k)
              a011 = vf%a011(i-1,j,k)
              a101 = vf%a101(i-1,j,k)
              d    = vf%dd  (i-1,j,k)
              phi  = vf%phi (i-1,j,k)
            end if

            xitg_a = ONE-dm%dt/dm%h(1)*u
            xitg_b = ONE
            yitg_a = ZERO
            yitg_b = ONE
            zitg_a = ZERO
            zitg_b = ONE

          else
            coef = -ONE
            if (i==dtmp%xsz(1) .and. dm%ibcx(2,6)/=IBC_PERIODIC) then
              coef = ZERO
            else
              lnx  = vf%lnx(i,j,k)
              lny  = vf%lny(i,j,k)
              lnz  = vf%lnz(i,j,k)
              a100 = vf%a100(i,j,k)
              a010 = vf%a010(i,j,k)
              a001 = vf%a001(i,j,k)
              a200 = vf%a200(i,j,k)
              a020 = vf%a020(i,j,k)
              a002 = vf%a002(i,j,k)
              a110 = vf%a110(i,j,k)
              a011 = vf%a011(i,j,k)
              a101 = vf%a101(i,j,k)
              d    = vf%dd  (i,j,k)
              phi  = vf%phi (i,j,k)
            end if

            xitg_a = ZERO
            xitg_b = -dm%dt/dm%h(1)*u
            yitg_a = ZERO
            yitg_b = ONE
            zitg_a = ZERO
            zitg_b = ONE

          end if

          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz))) then
            dirct = 1
          else if (abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
            dirct = 2
          else
            dirct = 3
          end if

          if (phi<vf%voflim .or. phi>(ONE-vf%voflim)) then
            vf%flx(i,j,k) = (xitg_b-xitg_a)*phi
          else if (abs(lnx)<eps .and. abs(lny)<eps .and. abs(lnz)<eps) then
            vf%flx(i,j,k) = (xitg_b-xitg_a)*phi
          else
            vf%flx(i,j,k) = cal_flux(a100, a010, a001,                         &
                                     a200, a020, a002,                         &
                                     a110, a011, a101,                         &
                                     d, beta,                                  &
                                     xitg_a, xitg_b,                           &
                                     yitg_a, yitg_b,                           &
                                     zitg_a, zitg_b,                           &
                                     dirct)
          end if

          vf%flx(i,j,k) = vf%flx(i,j,k)*dm%h(1)*coef

        enddo
      enddo
    enddo

    dtmp = dm%dccc
    do k=1, dtmp%xsz(3)
      do j=1, dtmp%xsz(2)
        do i=1, dtmp%xsz(1)
          if(i==dtmp%xsz(1) .and. dm%ibcx(2,6)==IBC_PERIODIC)then
            vf%phi(i,j,k) = vf%phi(i,j,k)-ONE/dm%h(1)*                         &
                            (vf%flx(1,j,k)-vf%flx(i,j,k))
            vf%phi(i,j,k) = vf%phi(i,j,k)/(ONE-dm%dt/dm%h(1)*                  &
                            (fl%qx(1,j,k)-fl%qx(i,j,k)))
            rhsx(i,j,k) = vf%phi(i,j,k)*(fl%qx(1,j,k)-fl%qx(i,j,k))
          else
            vf%phi(i,j,k) = vf%phi(i,j,k)-ONE/dm%h(1)*                         &
                            (vf%flx(i+1,j,k)-vf%flx(i,j,k))
            vf%phi(i,j,k) = vf%phi(i,j,k)/(ONE-dm%dt/dm%h(1)*                  &
                            (fl%qx(i+1,j,k)-fl%qx(i,j,k)))
            rhsx(i,j,k) = vf%phi(i,j,k)*(fl%qx(i+1,j,k)-fl%qx(i,j,k))
          end if
          vf%phi1(i,j,k) = vf%phi(i,j,k)
        enddo
      enddo
    enddo

    call clip_vof(dm, vf)

    !calculate the y-direction numerical flux over the y-pencil (parallel)
    call cal_colour_function(dm, vf)

    call transpose_x_to_y(fl%qy,   qy_cpcy,   dm%dcpc)
    call transpose_x_to_y(vf%phi,  phi_cccy,  dm%dccc)
    call transpose_x_to_y(vf%lnx,  lnx_cccy,  dm%dccc)
    call transpose_x_to_y(vf%lny,  lny_cccy,  dm%dccc)
    call transpose_x_to_y(vf%lnz,  lnz_cccy,  dm%dccc)
    call transpose_x_to_y(vf%a100, a100_cccy, dm%dccc)
    call transpose_x_to_y(vf%a010, a010_cccy, dm%dccc)
    call transpose_x_to_y(vf%a001, a001_cccy, dm%dccc)
    call transpose_x_to_y(vf%a200, a200_cccy, dm%dccc)
    call transpose_x_to_y(vf%a020, a020_cccy, dm%dccc)
    call transpose_x_to_y(vf%a002, a002_cccy, dm%dccc)
    call transpose_x_to_y(vf%a110, a110_cccy, dm%dccc)
    call transpose_x_to_y(vf%a011, a011_cccy, dm%dccc)
    call transpose_x_to_y(vf%a101, a101_cccy, dm%dccc)
    call transpose_x_to_y(vf%dd,   dd_cccy,   dm%dccc)

    dtmp = dm%dcpc
    do k=1, dtmp%ysz(3)
      do j=1, dtmp%ysz(2)
        do i=1, dtmp%ysz(1)
          v = qy_cpcy(i,j,k)
          lnx =  ZERO
          lny =  ZERO
          lnz =  ZERO
          a100 = ZERO
          a010 = ZERO
          a001 = ZERO
          a200 = ZERO
          a020 = ZERO
          a002 = ZERO
          a110 = ZERO
          a011 = ZERO
          a101 = ZERO
          d    = ZERO
          phi  = ZERO
          if (v>=ZERO) then
            coef = ONE
            if (j==1 .and. dm%ibcy(1,6)/=IBC_PERIODIC) then
              coef = ZERO
            else if (j==1 .and. dm%ibcy(1,6)==IBC_PERIODIC) then
              lnx  = lnx_cccy(i,dtmp%ysz(2),k)
              lny  = lny_cccy(i,dtmp%ysz(2),k)
              lnz  = lnz_cccy(i,dtmp%ysz(2),k)
              a100 = a100_cccy(i,dtmp%ysz(2),k)
              a010 = a010_cccy(i,dtmp%ysz(2),k)
              a001 = a001_cccy(i,dtmp%ysz(2),k)
              a200 = a200_cccy(i,dtmp%ysz(2),k)
              a020 = a020_cccy(i,dtmp%ysz(2),k)
              a002 = a002_cccy(i,dtmp%ysz(2),k)
              a110 = a110_cccy(i,dtmp%ysz(2),k)
              a011 = a011_cccy(i,dtmp%ysz(2),k)
              a101 = a101_cccy(i,dtmp%ysz(2),k)
              d    = dd_cccy  (i,dtmp%ysz(2),k)
              phi  = phi_cccy (i,dtmp%ysz(2),k)
            else
              lnx  = lnx_cccy(i,j-1,k)
              lny  = lny_cccy(i,j-1,k)
              lnz  = lnz_cccy(i,j-1,k)
              a100 = a100_cccy(i,j-1,k)
              a010 = a010_cccy(i,j-1,k)
              a001 = a001_cccy(i,j-1,k)
              a200 = a200_cccy(i,j-1,k)
              a020 = a020_cccy(i,j-1,k)
              a002 = a002_cccy(i,j-1,k)
              a110 = a110_cccy(i,j-1,k)
              a011 = a011_cccy(i,j-1,k)
              a101 = a101_cccy(i,j-1,k)
              d    = dd_cccy  (i,j-1,k)
              phi  = phi_cccy (i,j-1,k)
            end if

            xitg_a = ZERO
            xitg_b = ONE
            yitg_a = ONE-dm%dt/dm%h(2)*v
            yitg_b = ONE
            zitg_a = ZERO
            zitg_b = ONE

          else
            coef = -ONE
            if (j==dtmp%ysz(2) .and. dm%ibcy(2,6)/=IBC_PERIODIC) then
              coef = ZERO
            else
              lnx  =d- lnx_cccy(i,j,k)
              lny  = lny_cccy(i,j,k)
              lnz  = lnz_cccy(i,j,k)
              a100 = a100_cccy(i,j,k)
              a010 = a010_cccy(i,j,k)
              a001 = a001_cccy(i,j,k)
              a200 = a200_cccy(i,j,k)
              a020 = a020_cccy(i,j,k)
              a002 = a002_cccy(i,j,k)
              a110 = a110_cccy(i,j,k)
              a011 = a011_cccy(i,j,k)
              a101 = a101_cccy(i,j,k)
              d    = dd_cccy  (i,j,k)
              phi  = phi_cccy (i,j,k)
            end if

            xitg_a = ZERO
            xitg_b = ONE
            yitg_a = ZERO
            yitg_b = -dm%dt/dm%h(2)*v
            zitg_a = ZERO
            zitg_b = ONE

          end if

          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz))) then
            dirct = 1
          else if (abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
            dirct = 2
          else
            dirct = 3
          end if

          if (phi<vf%voflim .or. phi>(ONE-vf%voflim)) then
            gly_cpcy(i,j,k) = (yitg_b-yitg_a)*phi
          else if (abs(lnx)<eps .and. abs(lny)<eps .and. abs(lnz)<eps) then
            gly_cpcy(i,j,k) = (yitg_b-yitg_a)*phi
          else
            gly_cpcy(i,j,k) = cal_flux(a100, a010, a001,                       &
                                       a200, a020, a002,                       &
                                       a110, a011, a101,                       &
                                       d, beta,                                &
                                       xitg_a, xitg_b,                         &
                                       yitg_a, yitg_b,                         &
                                       zitg_a, zitg_b,                         &
                                       dirct)
          end if

          gly_cpcy(i,j,k) = gly_cpcy(i,j,k)*dm%h(2)*coef

        enddo
      enddo
    enddo

    dtmp = dm%dccc
    do k=1, dtmp%ysz(3)
      do j=1, dtmp%ysz(2)
        do i=1, dtmp%ysz(1)
          if(j==dtmp%ysz(2) .and. dm%ibcy(2,6)==IBC_PERIODIC) then
            phi_cccy(i,j,k) = phi_cccy(i,j,k)-ONE/dm%h(2)*                     &
                              (gly_cpcy(i,1,k)-gly_cpcy(i,j,k))
            phi_cccy(i,j,k) = phi_cccy(i,j,k)/(ONE-dm%dt/dm%h(2)*              &
                              (qy_cpcy(i,1,k)-qy_cpcy(i,j,k)))
            tmp_cccy(i,j,k) = phi_cccy(i,j,k)*(qy_cpcy(i,1,k)-qy_cpcy(i,j,k))
          else
            phi_cccy(i,j,k) = phi_cccy(i,j,k)-ONE/dm%h(2)*                     &
                              (gly_cpcy(i,j+1,k)-gly_cpcy(i,j,k))
            phi_cccy(i,j,k) = phi_cccy(i,j,k)/(ONE-dm%dt/dm%h(2)*              &
                              (qy_cpcy(i,j+1,k)-qy_cpcy(i,j,k)))
            tmp_cccy(i,j,k) = phi_cccy(i,j,k)*(qy_cpcy(i,j+1,k)-qy_cpcy(i,j,k))
          end if
          phi2_cccy(i,j,k) = phi_cccy(i,j,k)
        enddo
      enddo
    enddo

    call transpose_y_to_x(phi_cccy,  vf%phi,   dm%dccc)
    call transpose_y_to_x(phi2_cccy, vf%phi2,  dm%dccc)
    call transpose_y_to_x(tmp_cccy,  rhsy,     dm%dccc)

    call clip_vof(dm, vf)

    !calculate the z-direction numerical flux over the z-pencil (parallel)
    call cal_colour_function(dm, vf)

    call transpose_x_to_y(fl%qz,     qz_ccpy,   dm%dccp)
    call transpose_y_to_z(qz_ccpy,   qz_ccpz,   dm%dccp)
    call transpose_y_to_z(phi_cccy,  phi_cccz,  dm%dccc)
    call transpose_x_to_y(vf%lnx,    lnx_cccy,  dm%dccc)
    call transpose_y_to_z(lnx_cccy,  lnx_cccz,  dm%dccc)
    call transpose_x_to_y(vf%lny,    lny_cccy,  dm%dccc)
    call transpose_y_to_z(lny_cccy,  lny_cccz,  dm%dccc)
    call transpose_x_to_y(vf%lnz,    lnz_cccy,  dm%dccc)
    call transpose_y_to_z(lnz_cccy,  lnz_cccz,  dm%dccc)
    call transpose_x_to_y(vf%a100,   a100_cccy, dm%dccc)
    call transpose_y_to_z(a100_cccy, a100_cccz, dm%dccc)
    call transpose_x_to_y(vf%a010,   a010_cccy, dm%dccc)
    call transpose_y_to_z(a010_cccy, a010_cccz, dm%dccc)
    call transpose_x_to_y(vf%a001,   a001_cccy, dm%dccc)
    call transpose_y_to_z(a001_cccy, a001_cccz, dm%dccc)
    call transpose_x_to_y(vf%a200,   a200_cccy, dm%dccc)
    call transpose_y_to_z(a200_cccy, a200_cccz, dm%dccc)
    call transpose_x_to_y(vf%a020,   a020_cccy, dm%dccc)
    call transpose_y_to_z(a020_cccy, a020_cccz, dm%dccc)
    call transpose_x_to_y(vf%a002,   a002_cccy, dm%dccc)
    call transpose_y_to_z(a002_cccy, a002_cccz, dm%dccc)
    call transpose_x_to_y(vf%a110,   a110_cccy, dm%dccc)
    call transpose_y_to_z(a110_cccy, a110_cccz, dm%dccc)
    call transpose_x_to_y(vf%a011,   a011_cccy, dm%dccc)
    call transpose_y_to_z(a011_cccy, a011_cccz, dm%dccc)
    call transpose_x_to_y(vf%a101,   a101_cccy, dm%dccc)
    call transpose_y_to_z(a101_cccy, a101_cccz, dm%dccc)
    call transpose_x_to_y(vf%dd,     dd_cccy,   dm%dccc)
    call transpose_y_to_z(dd_cccy,   dd_cccz,   dm%dccc)

    dtmp = dm%dccp
    do k=1, dtmp%zsz(3)
      do j=1, dtmp%zsz(2)
        do i=1, dtmp%zsz(1)
          w = qz_ccpz(i,j,k)
          lnx =  ZERO
          lny =  ZERO
          lnz =  ZERO
          a100 = ZERO
          a010 = ZERO
          a001 = ZERO
          a200 = ZERO
          a020 = ZERO
          a002 = ZERO
          a110 = ZERO
          a011 = ZERO
          a101 = ZERO
          d    = ZERO
          phi  = ZERO
          if (w>=ZERO) then
            coef = ONE
            if (k==1 .and. dm%ibcz(1,6)/=IBC_PERIODIC) then
              coef = ZERO
            else if (k==1 .and. dm%ibcz(1,6)==IBC_PERIODIC) then
              lnx  = lnx_cccz(i,j,dtmp%zsz(3))
              lny  = lny_cccz(i,j,dtmp%zsz(3))
              lnz  = lnz_cccz(i,j,dtmp%zsz(3))
              a100 = a100_cccz(i,j,dtmp%zsz(3))
              a010 = a010_cccz(i,j,dtmp%zsz(3))
              a001 = a001_cccz(i,j,dtmp%zsz(3))
              a200 = a200_cccz(i,j,dtmp%zsz(3))
              a020 = a020_cccz(i,j,dtmp%zsz(3))
              a002 = a002_cccz(i,j,dtmp%zsz(3))
              a110 = a110_cccz(i,j,dtmp%zsz(3))
              a011 = a011_cccz(i,j,dtmp%zsz(3))
              a101 = a101_cccz(i,j,dtmp%zsz(3))
              d    = dd_cccz  (i,j,dtmp%zsz(3))
              phi  = phi_cccz (i,j,dtmp%zsz(3))
            else
              lnx  = lnx_cccz(i,j,k-1)
              lny  = lny_cccz(i,j,k-1)
              lnz  = lnz_cccz(i,j,k-1)
              a100 = a100_cccz(i,j,k-1)
              a010 = a010_cccz(i,j,k-1)
              a001 = a001_cccz(i,j,k-1)
              a200 = a200_cccz(i,j,k-1)
              a020 = a020_cccz(i,j,k-1)
              a002 = a002_cccz(i,j,k-1)
              a110 = a110_cccz(i,j,k-1)
              a011 = a011_cccz(i,j,k-1)
              a101 = a101_cccz(i,j,k-1)
              d    = dd_cccz  (i,j,k-1)
              phi  = phi_cccz (i,j,k-1)
            end if

            xitg_a = ZERO
            xitg_b = ONE
            yitg_a = ZERO
            yitg_b = ONE
            zitg_a = ONE-dm%dt/dm%h(3)*w
            zitg_b = ONE

          else
            coef = -ONE
            if (k==dtmp%zsz(3) .and. dm%ibcz(2,6)/=IBC_PERIODIC) then
              coef = ZERO
            else
              lnx  = lnx_cccz(i,j,k)
              lny  = lny_cccz(i,j,k)
              lnz  = lnz_cccz(i,j,k)
              a100 = a100_cccz(i,j,k)
              a010 = a010_cccz(i,j,k)
              a001 = a001_cccz(i,j,k)
              a200 = a200_cccz(i,j,k)
              a020 = a020_cccz(i,j,k)
              a002 = a002_cccz(i,j,k)
              a110 = a110_cccz(i,j,k)
              a011 = a011_cccz(i,j,k)
              a101 = a101_cccz(i,j,k)
              d    = dd_cccz  (i,j,k)
              phi  = phi_cccz (i,j,k)
            end if

            xitg_a = ZERO
            xitg_b = ONE
            yitg_a = ZERO
            yitg_b = ONE
            zitg_a = ZERO
            zitg_b = -dm%dt/dm%h(3)*w

          end if

          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz))) then
            dirct = 1
          else if (abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
            dirct = 2
          else
            dirct = 3
          end if

          if (phi<vf%voflim .or. phi>(ONE-vf%voflim)) then
            hlz_ccpz(i,j,k) = (zitg_b-zitg_a)*phi
          else if (abs(lnx)<eps .and. abs(lny)<eps .and. abs(lnz)<eps) then
            hlz_ccpz(i,j,k) = (zitg_b-zitg_a)*phi
          else
            hlz_ccpz(i,j,k) = cal_flux(a100, a010, a001,                       &
                                       a200, a020, a002,                       &
                                       a110, a011, a101,                       &
                                       d, beta,                                &
                                       xitg_a, xitg_b,                         &
                                       yitg_a, yitg_b,                         &
                                       zitg_a, zitg_b,                         &
                                       dirct)
          end if

          hlz_ccpz(i,j,k) = hlz_ccpz(i,j,k)*dm%h(3)*coef

        enddo
      enddo
    enddo

    dtmp = dm%dccc
    do k=1, dtmp%zsz(3)
      do j=1, dtmp%zsz(2)
        do i=1, dtmp%zsz(1)
          if(k==dtmp%zsz(3) .and. dm%ibcz(2,6)==IBC_PERIODIC)then
            phi_cccz(i,j,k) = phi_cccz(i,j,k)-ONE/dm%h(3)*                     &
                              (hlz_ccpz(i,j,1)-hlz_ccpz(i,j,k))
            phi_cccz(i,j,k) = phi_cccz(i,j,k)/(ONE-dm%dt/dm%h(3)*              &
                              (qz_ccpz(i,j,1)-qz_ccpz(i,j,k)))
            tmp_cccz(i,j,k) = phi_cccz(i,j,k)*(qz_ccpz(i,j,1)-qz_ccpz(i,j,k))
          else
            phi_cccz(i,j,k) = phi_cccz(i,j,k)-ONE/dm%h(3)*                     &
                              (hlz_ccpz(i,j,k+1)-hlz_ccpz(i,j,k))
            phi_cccz(i,j,k) = phi_cccz(i,j,k)/(ONE-dm%dt/dm%h(3)*              &
                              (qz_ccpz(i,j,k+1)-qz_ccpz(i,j,k)))
            tmp_cccz(i,j,k) = phi_cccz(i,j,k)*(qz_ccpz(i,j,k+1)-qz_ccpz(i,j,k))
          end if
          phi3_cccz(i,j,k) = phi_cccz(i,j,k)
        enddo
      enddo
    enddo

    call transpose_z_to_y(phi_cccz,  phi_cccy, dm%dccc)
    call transpose_y_to_x(phi_cccy,  vf%phi,   dm%dccc)
    call transpose_z_to_y(tmp_cccz,  tmp_cccy, dm%dccc)
    call transpose_y_to_x(tmp_cccy,  rhsz,     dm%dccc)
    call transpose_z_to_y(phi3_cccz, tmp_cccy, dm%dccc)
    call transpose_y_to_x(tmp_cccy,  vf%phi3,  dm%dccc)

    call clip_vof(dm, vf)

    do k=1, dtmp%xsz(3)
      do j=1, dtmp%xsz(2)
        do i=1, dtmp%xsz(1)
          vf%phi(i,j,k) = vf%phi3(i,j,k)-dm%dt*(rhsx(i,j,k)/dm%h(1)            &
                                               +rhsy(i,j,k)/dm%h(2)            &
                                               +rhsz(i,j,k)/dm%h(3))
        enddo
      enddo
    enddo

    call clip_vof(dm, vf)

    !update kappa
    call cal_colour_function(dm, vf)

    !update properties
    call Update_vof_properties(dm, fl, vf)



!    !calculate the y-direction numerical flux over the y-pencil (serial)
!    call cal_colour_function(dm, vf)

!    dtmp = dm%dcpc
!    do k=1, dtmp%ysz(3)
!      do j=1, dtmp%ysz(2)
!        do i=1, dtmp%ysz(1)
!          v = fl%qy(i,j,k)
!          lnx =  ZERO
!          lny =  ZERO
!          lnz =  ZERO
!          a100 = ZERO
!          a010 = ZERO
!          a001 = ZERO
!          a200 = ZERO
!          a020 = ZERO
!          a002 = ZERO
!          a110 = ZERO
!          a011 = ZERO
!          a101 = ZEro
!          d    = ZERO
!          phi  = ZERO
!          if (v>=ZERO) then
!            coef = ONE
!            if (j==1 .and. dm%ibcy(1,6)/=IBC_PERIODIC) then
!              coef = ZERO
!            else if (j==1 .and. dm%ibcy(1,6)==IBC_PERIODIC) then
!              lnx  = vf%lnx(i,dtmp%ysz(2),k)
!              lny  = vf%lny(i,dtmp%ysz(2),k)
!              lnz  = vf%lnz(i,dtmp%ysz(2),k)
!              a100 = vf%a100(i,dtmp%ysz(2),k)
!              a010 = vf%a010(i,dtmp%ysz(2),k)
!              a001 = vf%a001(i,dtmp%ysz(2),k)
!              a200 = vf%a200(i,dtmp%ysz(2),k)
!              a020 = vf%a020(i,dtmp%ysz(2),k)
!              a002 = vf%a002(i,dtmp%ysz(2),k)
!              a110 = vf%a110(i,dtmp%ysz(2),k)
!              a011 = vf%a011(i,dtmp%ysz(2),k)
!              a101 = vf%a101(i,dtmp%ysz(2),k)
!              d    = vf%dd  (i,dtmp%ysz(2),k)
!              phi  = vf%phi (i,dtmp%ysz(2),k)
!            else
!              lnx  = vf%lnx(i,j-1,k)
!              lny  = vf%lny(i,j-1,k)
!              lnz  = vf%lnz(i,j-1,k)
!              a100 = vf%a100(i,j-1,k)
!              a010 = vf%a010(i,j-1,k)
!              a001 = vf%a001(i,j-1,k)
!              a200 = vf%a200(i,j-1,k)
!              a020 = vf%a020(i,j-1,k)
!              a002 = vf%a002(i,j-1,k)
!              a110 = vf%a110(i,j-1,k)
!              a011 = vf%a011(i,j-1,k)
!              a101 = vf%a101(i,j-1,k)
!              d    = vf%dd  (i,j-1,k)
!              phi  = vf%phi (i,j-1,k)
!            end if

!            xitg_a = ZERO
!            xitg_b = ONE
!            yitg_a = ONE-dm%dt/dm%h(2)*v
!            yitg_b = ONE
!            zitg_a = ZERO
!            zitg_b = ONE

!          else
!            coef = -ONE
!            if (j==dtmp%ysz(2) .and. dm%ibcy(2,6)/=IBC_PERIODIC) then
!              coef = ZERO
!            else
!              lnx  = vf%lnx(i,j,k)
!              lny  = vf%lny(i,j,k)
!              lnz  = vf%lnz(i,j,k)
!              a100 = vf%a100(i,j,k)
!              a010 = vf%a010(i,j,k)
!              a001 = vf%a001(i,j,k)
!              a200 = vf%a200(i,j,k)
!              a020 = vf%a020(i,j,k)
!              a002 = vf%a002(i,j,k)
!              a110 = vf%a110(i,j,k)
!              a011 = vf%a011(i,j,k)
!              a101 = vf%a101(i,j,k)
!              d    = vf%dd  (i,j,k)
!              phi  = vf%phi (i,j,k)
!            end if

!            xitg_a = ZERO
!            xitg_b = ONE
!            yitg_a = ZERO
!            yitg_b = -dm%dt/dm%h(2)*v
!            zitg_a = ZERO
!            zitg_b = ONE

!          end if

!          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz))) then
!            dirct = 1
!          else if (abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
!            dirct = 2
!          else
!            dirct = 3
!          end if

!          if (phi<vf%voflim .or. phi>(ONE-vf%voflim)) then
!            vf%gly(i,j,k) = (yitg_b-yitg_a)*phi
!          else
!            vf%gly(i,j,k) = cal_flux(a100, a010, a001,                         &
!                                     a200, a020, a002,                         &
!                                     a110, a011, a101,                         &
!                                     d, beta,                                  &
!                                     xitg_a, xitg_b,                           &
!                                     yitg_a, yitg_b,                           &
!                                     zitg_a, zitg_b,                           &
!                                     dirct)
!          end if

!          vf%gly(i,j,k) = vf%gly(i,j,k)*dm%h(2)*coef

!        enddo
!      enddo
!    enddo

!    dtmp = dm%dccc
!    do k=1, dtmp%ysz(3)
!      do j=1, dtmp%ysz(2)
!        do i=1, dtmp%ysz(1)
!          if(j==dtmp%ysz(2) .and. dm%ibcy(2,6)==IBC_PERIODIC)then
!            vf%phi(i,j,k) = vf%phi(i,j,k)-ONE/dm%h(2)                          &
!                          * (vf%gly(i,1,k)-vf%gly(i,j,k))
!            vf%phi(i,j,k) = vf%phi(i,j,k)/(ONE-dm%dt/dm%h(2)                   &
!                          * (fl%qy(i,1,k)-fl%qy(i,j,k)))
!          else
!            vf%phi(i,j,k) = vf%phi(i,j,k)-ONE/dm%h(2)                          &
!                          * (vf%gly(i,j+1,k)-vf%gly(i,j,k))
!            vf%phi(i,j,k) = vf%phi(i,j,k)/(ONE-dm%dt/dm%h(2)                   &
!                          * (fl%qy(i,j+1,k)-fl%qy(i,j,k)))
!          end if
!          vf%phi2(i,j,k) = vf%phi(i,j,k)
!        enddo
!      enddo
!    enddo

!    !calculate the z-direction numerical flux over the z-pencil (serial)
!    call cal_colour_function(dm, vf)

!    dtmp = dm%dccp
!    do k=1, dtmp%zsz(3)
!      do j=1, dtmp%zsz(2)
!        do i=1, dtmp%zsz(1)
!          w = fl%qz(i,j,k)
!          lnx =  ZERO
!          lny =  ZERO
!          lnz =  ZERO
!          a100 = ZERO
!          a010 = ZERO
!          a001 = ZERO
!          a200 = ZERO
!          a020 = ZERO
!          a002 = ZERO
!          a110 = ZERO
!          a011 = ZERO
!          a101 = ZERO
!          d    = ZERO
!          phi  = ZERO
!          if (w>=ZERO) then
!            coef = ONE
!            if (k==1 .and. dm%ibcz(1,6)/=IBC_PERIODIC) then
!              coef = ZERO
!            else if (k==1 .and. dm%ibcz(1,6)==IBC_PERIODIC) then
!              lnx  = vf%lnx(i,j,dtmp%zsz(3))
!              lny  = vf%lny(i,j,dtmp%zsz(3))
!              lnz  = vf%lnz(i,j,dtmp%zsz(3))
!              a100 = vf%a100(i,j,dtmp%zsz(3))
!              a010 = vf%a010(i,j,dtmp%zsz(3))
!              a001 = vf%a001(i,j,dtmp%zsz(3))
!              a200 = vf%a200(i,j,dtmp%zsz(3))
!              a020 = vf%a020(i,j,dtmp%zsz(3))
!              a002 = vf%a002(i,j,dtmp%zsz(3))
!              a110 = vf%a110(i,j,dtmp%zsz(3))
!              a011 = vf%a011(i,j,dtmp%zsz(3))
!              a101 = vf%a101(i,j,dtmp%zsz(3))
!              d    = vf%dd  (i,j,dtmp%zsz(3))
!              phi  = vf%phi (i,j,dtmp%zsz(3))
!            else
!              lnx  = vf%lnx(i,j,k-1)
!              lny  = vf%lny(i,j,k-1)
!              lnz  = vf%lnz(i,j,k-1)
!              a100 = vf%a100(i,j,k-1)
!              a010 = vf%a010(i,j,k-1)
!              a001 = vf%a001(i,j,k-1)
!              a200 = vf%a200(i,j,k-1)
!              a020 = vf%a020(i,j,k-1)
!              a002 = vf%a002(i,j,k-1)
!              a110 = vf%a110(i,j,k-1)
!              a011 = vf%a011(i,j,k-1)
!              a101 = vf%a101(i,j,k-1)
!              d    = vf%dd  (i,j,k-1)
!              phi  = vf%phi (i,j,k-1)
!            end if

!            xitg_a = ZERO
!            xitg_b = ONE
!            yitg_a = ZERO
!            yitg_b = ONE
!            zitg_a = ONE-dm%dt/dm%h(3)*w
!            zitg_b = ONE

!          else
!            coef = -ONE
!            if (k==dtmp%zsz(3) .and. dm%ibcz(2,6)/=IBC_PERIODIC) then
!              coef = ZERO
!            else
!              lnx  = vf%lnx(i,j,k)
!              lny  = vf%lny(i,j,k)
!              lnz  = vf%lnz(i,j,k)
!              a100 = vf%a100(i,j,k)
!              a010 = vf%a010(i,j,k)
!              a001 = vf%a001(i,j,k)
!              a200 = vf%a200(i,j,k)
!              a020 = vf%a020(i,j,k)
!              a002 = vf%a002(i,j,k)
!              a110 = vf%a110(i,j,k)
!              a011 = vf%a011(i,j,k)
!              a101 = vf%a101(i,j,k)
!              d    = vf%dd  (i,j,k)
!              phi  = vf%phi (i,j,k)
!            end if

!            xitg_a = ZERO
!            xitg_b = ONE
!            yitg_a = ZERO
!            yitg_b = ONE
!            zitg_a = ZERO
!            zitg_b = -dm%dt/dm%h(3)*w

!          end if

!          if (abs(lnx)==max(abs(lnx),abs(lny),abs(lnz))) then
!            dirct = 1
!          else if (abs(lny)==max(abs(lnx),abs(lny),abs(lnz))) then
!            dirct = 2
!          else
!            dirct = 3
!          end if

!          if (phi<vf%voflim .or. phi>(ONE-vf%voflim)) then
!            vf%hlz(i,j,k) = (zitg_b-zitg_a)*phi
!          else
!            vf%hlz(i,j,k) = cal_flux(a100, a010, a001,                         &
!                                     a200, a020, a002,                         &
!                                     a110, a011, a101,                         &
!                                     d, beta,                                  &
!                                     xitg_a, xitg_b,                           &
!                                     yitg_a, yitg_b,                           &
!                                     zitg_a, zitg_b,                           &
!                                     dirct)
!          end if

!          vf%hlz(i,j,k) = vf%hlz(i,j,k)*dm%h(3)*coef

!        enddo
!      enddo
!    enddo

!    dtmp = dm%dccc
!    do k=1, dtmp%zsz(3)
!      do j=1, dtmp%zsz(2)
!        do i=1, dtmp%zsz(1)
!          if(k==dtmp%zsz(3) .and. dm%ibcz(2,6)==IBC_PERIODIC)then
!            vf%phi(i,j,k) = vf%phi(i,j,k)-ONE/dm%h(3)                          &
!                          * (vf%hlz(i,j,1)-vf%hlz(i,j,k))
!            vf%phi(i,j,k) = vf%phi(i,j,k)/(ONE-dm%dt/dm%h(3)                   &
!                          * (fl%qz(i,j,1)-fl%qz(i,j,k)))
!          else
!            vf%phi(i,j,k) = vf%phi(i,j,k)-ONE/dm%h(3)                          &
!                          * (vf%hlz(i,j,k+1)-vf%hlz(i,j,k))
!            vf%phi(i,j,k) = vf%phi(i,j,k)/(ONE-dm%dt/dm%h(3)                   &
!                          * (fl%qz(i,j,k+1)-fl%qz(i,j,k)))
!          end if
!          vf%phi3(i,j,k) = vf%phi(i,j,k)
!        enddo
!      enddo
!    enddo

!    do k=1, dtmp%xsz(3)
!      do j=1, dtmp%xsz(2)
!        do i=1, dtmp%xsz(1)
!          vf%phi(i,j,k) = vf%phi3(i,j,k)
!          if(i==dtmp%xsz(1) .and. dm%ibcx(2,6)==IBC_PERIODIC)then
!            vf%phi(i,j,k) = vf%phi(i,j,k)-dm%dt*vf%phi1(i,j,k)*                &
!                            (fl%qx(1,j,k)-fl%qx(i,j,k))/dm%h(1)
!          else
!            vf%phi(i,j,k) = vf%phi(i,j,k)-dm%dt*vf%phi1(i,j,k)*                &
!                            (fl%qx(i+1,j,k)-fl%qx(i,j,k))/dm%h(1)          
!          end if
!          if(j==dtmp%ysz(2) .and. dm%ibcy(2,6)==IBC_PERIODIC)then
!            vf%phi(i,j,k) = vf%phi(i,j,k)-dm%dt*vf%phi2(i,j,k)*                &
!                            (fl%qy(i,1,k)-fl%qy(i,j,k))/dm%h(2)          
!          else
!            vf%phi(i,j,k) = vf%phi(i,j,k)-dm%dt*vf%phi2(i,j,k)*                &
!                            (fl%qy(i,j+1,k)-fl%qy(i,j,k))/dm%h(2)          
!          end if
!          if(k==dtmp%zsz(3) .and. dm%ibcz(2,6)==IBC_PERIODIC)then
!            vf%phi(i,j,k) = vf%phi(i,j,k)-dm%dt*vf%phi3(i,j,k)*                &
!                            (fl%qz(i,j,1)-fl%qz(i,j,k))/dm%h(3)           
!          else
!            vf%phi(i,j,k) = vf%phi(i,j,k)-dm%dt*vf%phi3(i,j,k)*                &
!                            (fl%qz(i,j,k+1)-fl%qz(i,j,k))/dm%h(3)           
!          end if
!          if(vf%phi(i,j,k)<ZERO) vf%phi(i,j,k) = ZERO
!          if(vf%phi(i,j,k)>ONE) vf%phi(i,j,k) = ONE
!        enddo
!      enddo
!    enddo

  end subroutine Solve_vof_eq

end module
