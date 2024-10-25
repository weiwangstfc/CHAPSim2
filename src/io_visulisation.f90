module io_visulisation_mod
  use io_tools_mod
  use parameters_constant_mod
  use print_msg_mod
  implicit none 

  character(len=*), parameter :: io_name = "solution-io"

  integer, parameter :: XDMF_HEADER = 1, &
                        XDMF_FOOTER = 2
  integer, parameter :: PLANE_AVERAGE = -1
  character(6), parameter :: SCALAR = "Scalar", &
                             VECTOR = "Vector", &
                             TENSOR = "Tensor"
  character(4), parameter :: CELL = "Cell", &
                             NODE = "Node"
  

  integer, allocatable :: nnd_visu(:, :)
  integer, allocatable :: ncl_visu(:, :)

  real(WP), allocatable :: xp(:), yp(:), zp(:)

  !character(6)  :: svisudim

  private :: write_visu_headerfooter
  private :: write_visu_field
  private :: visu_average_periodic_data
  private :: write_visu_profile

  public  :: write_visu_ini
  public  :: write_visu_flow
  public  :: write_visu_thermo
  public  :: write_visu_any3darray

  public  :: write_visu_stats_flow
  public  :: write_visu_stats_thermo
  
  
  
contains
!==========================================================================================================
! xszV means:
! x - xpencil, could also be y, z
! sz - size, could also be st, en
! V - visulisation
! xszV is the same as dppp%xsz/nskip, based on nodes, considering the skip nodes
!==========================================================================================================
  subroutine write_visu_ini(dm)
    use udf_type_mod
    use parameters_constant_mod, only: MAXP
    use decomp_2d, only: xszV, yszV, zszV
    use io_files_mod
    implicit none 
    type(t_domain), intent(in) :: dm

    integer :: i, j ,k
    
    character(120):: grid_flname
    character(120):: keyword
    integer :: iogrid

!----------------------------------------------------------------------------------------------------------
! allocate
!----------------------------------------------------------------------------------------------------------
    if(.not. allocated(nnd_visu)) allocate (nnd_visu(3, nxdomain))
    if(.not. allocated(ncl_visu)) allocate (ncl_visu(3, nxdomain))
    nnd_visu = 0
    ncl_visu = 0
!----------------------------------------------------------------------------------------------------------
! global size
!----------------------------------------------------------------------------------------------------------
    !svisudim = ''
    !if(dm%visu_idim == Ivisudim_3D) then
      !svisudim = "3d"
      nnd_visu(1, dm%idom) = xszV(1)
      nnd_visu(2, dm%idom) = yszV(2)
      nnd_visu(3, dm%idom) = zszV(3)
      do i = 1, 3
        if(dm%is_periodic(i)) then 
          ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
          nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
        else 
          ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
        end if
      end do
    ! !else if(dm%visu_idim == Ivisudim_2D_Xa) then
    !   svisudim = "2d_xa"
    !   nnd_visu(1, dm%idom) = 1
    !   nnd_visu(2, dm%idom) = yszV(2)
    !   nnd_visu(3, dm%idom) = zszV(3)
    !   do i = 1, 3
    !     if(dm%is_periodic(i)) then 
    !       ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
    !       nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
    !     else 
    !       ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
    !     end if
    !   end do
    ! else if(dm%visu_idim == Ivisudim_2D_Ya) then
    !   svisudim = "2d_ya"
    !   nnd_visu(1, dm%idom) = xszV(1)
    !   nnd_visu(2, dm%idom) = 1
    !   nnd_visu(3, dm%idom) = zszV(3)
    !   do i = 1, 3
    !     if(dm%is_periodic(i)) then 
    !       ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
    !       nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
    !     else 
    !       ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
    !     end if
    !   end do
    ! else if(dm%visu_idim == Ivisudim_2D_Za) then
    !   svisudim = "2d_za"
    !   nnd_visu(1, dm%idom) = xszV(1)
    !   nnd_visu(2, dm%idom) = yszV(2)
    !   nnd_visu(3, dm%idom) = 1
    !   do i = 1, 3
    !     if(dm%is_periodic(i)) then 
    !       ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
    !       nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
    !     else 
    !       ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
    !     end if
    !   end do
    ! else if(dm%visu_idim == Ivisudim_1D_XZa) then
    !   svisudim = "2d_xza"
    !   nnd_visu(1, dm%idom) = 1
    !   nnd_visu(2, dm%idom) = yszV(2)
    !   nnd_visu(3, dm%idom) = 1
    !   do i = 1, 3
    !     if(dm%is_periodic(i)) then 
    !       ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
    !       nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
    !     else 
    !       ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
    !     end if
    !   end do
    ! else
    !   svisudim = "3d"
    !   nnd_visu(1, dm%idom) = xszV(1)
    !   nnd_visu(2, dm%idom) = yszV(2)
    !   nnd_visu(3, dm%idom) = zszV(3)
    !   do i = 1, 3
    !     if(dm%is_periodic(i)) then 
    !       ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
    !       nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
    !     else 
    !       ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
    !     end if
    !   end do
    ! end if
!----------------------------------------------------------------------------------------------------------
! write grids
!----------------------------------------------------------------------------------------------------------    
    if(nrank == 0) then

      if(.not. allocated(xp)) allocate ( xp(nnd_visu(1, dm%idom)) )
      if(.not. allocated(yp)) allocate ( yp(nnd_visu(2, dm%idom)) )
      if(.not. allocated(zp)) allocate ( zp(nnd_visu(3, dm%idom)) )

      xp = MAXP
      yp = MAXP
      zp = MAXP

      do i = 1, nnd_visu(1, dm%idom)
        xp(i) = real(i-1, WP) * dm%h(1) * dm%visu_nskip(1)
      enddo
      do j = 1, nnd_visu(2, dm%idom)
        if(dm%is_stretching(2)) then 
        yp(j) = dm%yp(j)
        else 
          yp(j) = real(j-1, WP) * dm%h(2) * dm%visu_nskip(2)
        end if
      end do
      do k = 1, nnd_visu(3, dm%idom)
        zp(k) = real(k-1, WP) * dm%h(3) * dm%visu_nskip(3)
      enddo


      keyword = "grid_x"
      call generate_pathfile_name(grid_flname, dm%idom, keyword, dir_visu, 'dat')
      open(newunit = iogrid, file = trim(grid_flname), action = "write", status="replace")
      write(iogrid, *) xp
      close(iogrid)

      keyword = "grid_y"
      call generate_pathfile_name(grid_flname, dm%idom, keyword, dir_visu, 'dat')
      open(newunit = iogrid, file = trim(grid_flname), action = "write", status="replace")
      write(iogrid, *) yp
      close(iogrid)

      keyword = "grid_z"
      call generate_pathfile_name(grid_flname, dm%idom, keyword, dir_visu, 'dat')
      open(newunit = iogrid, file = trim(grid_flname), action = "write", status="replace")
      write(iogrid, *) zp
      close(iogrid)

    end if

    return
  end subroutine

!==========================================================================================================
! ref: https://www.xdmf.org/index.php/XDMF_Model_and_Format
!==========================================================================================================
  subroutine write_visu_headerfooter(dm, visuname, iheadfoot, iter)
    use precision_mod
    use parameters_constant_mod, only: MAXP
    use udf_type_mod, only: t_domain
    use decomp_2d, only: nrank, mytype, xszV, yszV, zszV
    use io_files_mod
    implicit none 
    integer, intent(in)        :: iheadfoot
    integer, intent(in)        :: iter
    type(t_domain), intent(in) :: dm
    character(*), intent(in)   :: visuname

    character(120):: keyword
    character(120):: visu_flname
    !character(120):: grid_flname(3)
    !character(1)  :: str(3)

    integer :: ioxdmf
    
    !integer :: i, j, k
    if(nrank /= 0) return
!----------------------------------------------------------------------------------------------------------
! visu file name
!----------------------------------------------------------------------------------------------------------
    keyword = trim(visuname)
    call generate_pathfile_name(visu_flname, dm%idom, keyword, dir_visu, 'xdmf', iter)
    open(newunit = ioxdmf, file = trim(visu_flname), action = "write", position="append")
!----------------------------------------------------------------------------------------------------------
! xdmf head
!----------------------------------------------------------------------------------------------------------
    if(iheadfoot == XDMF_HEADER) then
!----------------------------------------------------------------------------------------------------------
! grid file
! to do: check why the binary/ascii file does not work.
!----------------------------------------------------------------------------------------------------------
      ! str(1) = "x"
      ! str(2) = "y"
      ! str(3) = "z"
      ! do i = 1, 3
      !   keyword = trim(svisudim)//"_grid_"//trim(str(i))
      !   call generate_file_name(grid_flname(i), dm%idom, keyword, 'dat')
      !   if(.not.file_exists(trim(trim(dir_visu)//"/"//trim(grid_flname(i))))) then
      !     call Print_error_msg("Mesh file for visu does not exist. Filename = "//trim(trim(dir_visu)//"/"//trim(grid_flname(i))))
      !   end if
      ! end do

!----------------------------------------------------------------------------------------------------------
! write header
! geometry is based on node coordinates
! to do: write mesh into mesh.bin 
!----------------------------------------------------------------------------------------------------------
      write(ioxdmf, '(A)')'<?xml version="1.0" ?>'
      write(ioxdmf, '(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(ioxdmf, '(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      write(ioxdmf, *)' <Domain>'
      write(ioxdmf, *)'   <Grid Name="'//trim(keyword)//'" GridType="Uniform">'
      write(ioxdmf, *)'     <Topology name="topo" TopologyType="3DRectMesh"'
      write(ioxdmf, *)'        Dimensions="', nnd_visu(3, dm%idom), nnd_visu(2, dm%idom), nnd_visu(1, dm%idom),'">'
      write(ioxdmf, *)'     </Topology>'
      write(ioxdmf, *)'     <Geometry name="geo" GeometryType="VXVYVZ">'
      write(ioxdmf, *)'        <DataItem Format="XML"'
      write(ioxdmf, *)'            NumberType="Float" Precision="8" '
      write(ioxdmf, *)'            Dimensions="', nnd_visu(1, dm%idom), '">'
      write(ioxdmf, *)'              ', xp(1:nnd_visu(1, dm%idom))
      write(ioxdmf, *)'        </DataItem>'
      write(ioxdmf, *)'        <DataItem Format="XML"'
      write(ioxdmf, *)'            NumberType="Float" Precision="8" '
      write(ioxdmf, *)'            Dimensions="', nnd_visu(2, dm%idom), '">'
      write(ioxdmf, *)'              ', yp(1:nnd_visu(2, dm%idom))
      write(ioxdmf, *)'        </DataItem>'
      write(ioxdmf, *)'        <DataItem Format="XML"'
      write(ioxdmf, *)'            NumberType="Float" Precision="8" '
      write(ioxdmf, *)'            Dimensions="', nnd_visu(3, dm%idom), '">'
      write(ioxdmf, *)'              ', zp(1:nnd_visu(3, dm%idom))
      write(ioxdmf, *)'        </DataItem>'
      write(ioxdmf, *)'      </Geometry>'
    else if (iheadfoot == XDMF_FOOTER) then 
      write(ioxdmf, *)'   </Grid>'
      write(ioxdmf, *)' </Domain>'
      write(ioxdmf, '(A)')'</Xdmf>'
    else 
    end if
    close(ioxdmf)
    return
  end subroutine
!==========================================================================================================
! ref: https://www.xdmf.org/index.php/XDMF_Model_and_Format
!==========================================================================================================
  subroutine write_visu_field(dm, var, dtmp, varname, visuname, attributetype, centring, iter)
    use precision_mod
    use decomp_2d
    use decomp_2d_io
    use udf_type_mod, only: t_domain
    use io_files_mod
    use decomp_operation_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP), intent(in) :: var(:, :, :)
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: visuname
    character(*), intent(in) :: attributetype
    character(*), intent(in) :: centring
    type(DECOMP_INFO), intent(in) :: dtmp
    integer, intent(in), optional :: iter

    character(120):: data_flname
    character(120):: data_flname_path
    character(120):: visu_flname_path
    character(120):: keyword
    integer :: nsz(3)
    integer :: ioxdmf

    if((.not. is_same_decomp(dtmp, dm%dccc))) then
      if(nrank == 0) call Print_error_msg("Data is not stored at cell centre. varname = " // trim(varname))
    end if
!----------------------------------------------------------------------------------------------------------
! xmdf file name
!----------------------------------------------------------------------------------------------------------
    keyword = trim(visuname)
    call generate_pathfile_name(visu_flname_path, dm%idom, keyword, dir_visu, 'xdmf', iter)
!----------------------------------------------------------------------------------------------------------
! write data into binary file
!----------------------------------------------------------------------------------------------------------
    if(dm%visu_idim == Ivisudim_3D) then
      keyword = trim(varname)
      call generate_pathfile_name(data_flname_path, dm%idom, keyword, dir_data, 'bin', iter)
      if(.not.file_exists(data_flname_path)) &
      call decomp_2d_write_one(X_PENCIL, var, trim(data_flname_path), dtmp)

    else if(dm%visu_idim == Ivisudim_1D_XZa) then
      !to add 1D profile
    else 
      keyword = trim(varname)
      call generate_file_name(data_flname, dm%idom, keyword, 'bin', iter)
      call generate_pathfile_name(data_flname_path, dm%idom, keyword, dir_data, 'bin', iter)
      call decomp_2d_write_plane(X_PENCIL, var, dm%visu_idim, PLANE_AVERAGE, trim(dir_data), trim(data_flname), io_name, dtmp)
    end if
!----------------------------------------------------------------------------------------------------------
! dataitem for xdmf file
!----------------------------------------------------------------------------------------------------------
    if (nrank == 0) then
      open(newunit = ioxdmf, file = trim(visu_flname_path), action = "write", status = "old", position = "append")

      if(trim(centring) == TRIM(CELL)) then
        nsz(1:3) = ncl_visu(1:3, dm%idom)
      else if (trim(centring) == TRIM(NODE)) then
        nsz(1:3) = nnd_visu(1:3, dm%idom)
      else
      end if
      write(ioxdmf, *)'      <Attribute Name="'//trim(varname)// &
                            '" AttributeType="'//trim(attributetype)// &
                            '" Center="'//trim(centring)//'">'
      write(ioxdmf, *)'           <DataItem Format="Binary"'
      write(ioxdmf, *)'            NumberType="Float" Precision="8" Endian="little" Seek="0"'
      write(ioxdmf, *)'            Dimensions="', nsz(3), nsz(2), nsz(1), '">'
      write(ioxdmf, *)'              '//"../"//trim(data_flname_path)
      write(ioxdmf, *)'           </DataItem>'
      write(ioxdmf, *)'        </Attribute>'
      close(ioxdmf)
    end if

    return
  end subroutine 
  !==========================================================================================================
! ref: https://www.xdmf.org/index.php/XDMF_Model_and_Format
!==========================================================================================================
  subroutine write_visu_profile(dm, var, dtmp, varname, visuname, attributetype, centring, idim, iter)
    use precision_mod
    use decomp_2d
    use decomp_2d_io
    use udf_type_mod, only: t_domain
    use io_files_mod
    use decomp_operation_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP), intent(in) :: var(:)
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: visuname
    character(*), intent(in) :: attributetype
    character(*), intent(in) :: centring
    type(DECOMP_INFO), intent(in) :: dtmp
    integer, intent(in) :: idim
    integer, intent(in), optional :: iter

    character(120):: data_flname
    character(120):: data_flname_path
    character(120):: visu_flname_path
    character(120):: keyword
    integer :: nsz(3)
    integer :: ioxdmf, iofl

    integer :: j

    if((.not. is_same_decomp(dtmp, dm%dccc))) then
      if(nrank == 0) call Print_error_msg("Data is not stored at cell centre. varname = " // trim(varname))
    end if
!----------------------------------------------------------------------------------------------------------
! write data 
!----------------------------------------------------------------------------------------------------------
    keyword = trim(varname)
    call generate_pathfile_name(data_flname_path, dm%idom, keyword, dir_data, 'dat', iter)
    open(newunit = iofl, file = data_flname_path, action = "write", status="replace")
    if(idim /= 2) call Print_error_msg('Error in direction')
    do j = 1, dtmp%ysz(2)
      write(iofl, *) j, dm%yc(j), var(j) 
    end do

    return
  end subroutine 
!==========================================================================================================
  subroutine write_visu_flow(fl, dm, str)
    use udf_type_mod
    use precision_mod
    use operations
    use parameters_constant_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_flow), intent(in) :: fl
    character(4), intent(in), optional :: str

    integer :: iter 
    character(120) :: visuname
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil

    iter = fl%iteration
    visuname = 'flow'
    if(present(str)) visuname = trim(visuname)//'_'//trim(str)
!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_HEADER, iter)
!----------------------------------------------------------------------------------------------------------
! write data, pressure, to cell centre
!----------------------------------------------------------------------------------------------------------
    call write_visu_field(dm, fl%pres, dm%dccc, "pr_visu", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! qx, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(fl%qx, accc, dm, dm%iAccuracy, dm%ibcx_qx(:), dm%fbcx_qx)
    call write_visu_field(dm, accc, dm%dccc, "qx_visu", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! qy, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy(:), dm%fbcy_qy)
    call transpose_y_to_x(accc_ypencil, accc, dm%dccc)
    call write_visu_field(dm, accc, dm%dccc, "qy_visu", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! qz, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qz, accp_ypencil, dm%dccp)
    call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
    call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz(:))
    call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
    call transpose_y_to_x(accc_ypencil, accc, dm%dccc)
    call write_visu_field(dm, accc, dm%dccc, "qz_visu", trim(visuname), SCALAR, CELL, iter)

    if(dm%is_thermo) then
!----------------------------------------------------------------------------------------------------------
! gx, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
      call Get_x_midp_P2C_3D(fl%gx, accc, dm, dm%iAccuracy, dm%ibcx_qx(:), dm%fbcx_gx)
      call write_visu_field(dm, accc, dm%dccc, "gx_visu", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! gy, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
      call transpose_x_to_y(fl%gy, acpc_ypencil, dm%dcpc)
      call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy(:), dm%fbcy_gy)
      call transpose_y_to_x(accc_ypencil, accc, dm%dccc)
      call write_visu_field(dm, accc, dm%dccc, "gy_visu", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! gz, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
      call transpose_x_to_y(fl%gz, accp_ypencil, dm%dccp)
      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
      call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz(:), dm%fbcz_gz)
      call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
      call transpose_y_to_x(accc_ypencil, accc, dm%dccc)
      call write_visu_field(dm, accc, dm%dccc, "gz_visu", trim(visuname), SCALAR, CELL, iter)
    end if
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_FOOTER, iter)

    if(nrank == 0) call Print_debug_mid_msg("Write out visulisation for flow field.")
    
    return
  end subroutine

  !==========================================================================================================
  subroutine write_visu_thermo(tm, fl, dm, str)
    use udf_type_mod
    use precision_mod
    use operations
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm
    type(t_flow),   intent(in) :: fl
    character(4), intent(in), optional :: str

    integer :: iter 
    character(120) :: visuname

    iter = tm%iteration
    visuname = 'thermo'
    if(present(str)) visuname = trim(visuname)//'_'//trim(str)
!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_HEADER, iter)
!----------------------------------------------------------------------------------------------------------
! write data, temperature, to cell centre
!----------------------------------------------------------------------------------------------------------
    call write_visu_field(dm, tm%tTemp, dm%dccc, "Temp_visu", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%dDens, dm%dccc, "Dens_visu", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%mVisc, dm%dccc, "Visc_visu", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, tm%kCond, dm%dccc, "Cond_visu", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, tm%hEnth, dm%dccc, "Enth_visu", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_FOOTER, iter)

    if(nrank == 0) call Print_debug_mid_msg("Write out visulisation for thermal field.")
    
    return
  end subroutine

!==========================================================================================================
  subroutine write_visu_stats_flow(fl, dm)
    use udf_type_mod
    use precision_mod
    use operations
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(in) :: fl

    integer :: iter 
    character(120) :: visuname

!==========================================================================================================
! write time averaged 3d data
!==========================================================================================================
    iter = fl%iteration
    visuname = 'time_averaged_flow'
!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_HEADER, iter)
!----------------------------------------------------------------------------------------------------------
! write data, 
!----------------------------------------------------------------------------------------------------------
    call write_visu_field(dm, fl%pr_mean,                     dm%dccc, "time_averaged_pr", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%u_vector_mean  (:, :, :, 1), dm%dccc, "time_averaged_ux", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%u_vector_mean  (:, :, :, 2), dm%dccc, "time_averaged_uy", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%u_vector_mean  (:, :, :, 3), dm%dccc, "time_averaged_uz", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%uu_tensor6_mean(:, :, :, 1), dm%dccc, "time_averaged_uu", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%uu_tensor6_mean(:, :, :, 2), dm%dccc, "time_averaged_vv", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%uu_tensor6_mean(:, :, :, 3), dm%dccc, "time_averaged_ww", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%uu_tensor6_mean(:, :, :, 4), dm%dccc, "time_averaged_uv", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%uu_tensor6_mean(:, :, :, 5), dm%dccc, "time_averaged_uw", trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, fl%uu_tensor6_mean(:, :, :, 6), dm%dccc, "time_averaged_vw", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_FOOTER, iter)
!==========================================================================================================
! write time averaged and space averaged 3d data (stored 2d or 1d data)
!==========================================================================================================
    if( ANY(dm%is_periodic(:))) then

    iter = fl%iteration
    visuname = 'time_space_averaged_flow'
!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    if(.not. (dm%is_periodic(1) .and. dm%is_periodic(3))) &
    call write_visu_headerfooter(dm, trim(visuname), XDMF_HEADER, iter)
!----------------------------------------------------------------------------------------------------------
! write data, 
!----------------------------------------------------------------------------------------------------------
    call visu_average_periodic_data(                    fl%pr_mean, dm%dccc, dm, "time_space_averaged_pr", trim(visuname), iter)
    call visu_average_periodic_data(  fl%u_vector_mean(:, :, :, 1), dm%dccc, dm, "time_space_averaged_ux", trim(visuname), iter)
    call visu_average_periodic_data(  fl%u_vector_mean(:, :, :, 2), dm%dccc, dm, "time_space_averaged_uy", trim(visuname), iter)
    call visu_average_periodic_data(  fl%u_vector_mean(:, :, :, 3), dm%dccc, dm, "time_space_averaged_uz", trim(visuname), iter)
    call visu_average_periodic_data(fl%uu_tensor6_mean(:, :, :, 1), dm%dccc, dm, "time_space_averaged_uu", trim(visuname), iter)
    call visu_average_periodic_data(fl%uu_tensor6_mean(:, :, :, 2), dm%dccc, dm, "time_space_averaged_vv", trim(visuname), iter)
    call visu_average_periodic_data(fl%uu_tensor6_mean(:, :, :, 3), dm%dccc, dm, "time_space_averaged_ww", trim(visuname), iter)
    call visu_average_periodic_data(fl%uu_tensor6_mean(:, :, :, 4), dm%dccc, dm, "time_space_averaged_uv", trim(visuname), iter)
    call visu_average_periodic_data(fl%uu_tensor6_mean(:, :, :, 5), dm%dccc, dm, "time_space_averaged_uw", trim(visuname), iter)
    call visu_average_periodic_data(fl%uu_tensor6_mean(:, :, :, 6), dm%dccc, dm, "time_space_averaged_vw", trim(visuname), iter)
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    if(.not. (dm%is_periodic(1) .and. dm%is_periodic(3))) &
    call write_visu_headerfooter(dm, trim(visuname), XDMF_FOOTER, iter)

    end if

    return
  end subroutine

  !==========================================================================================================
  subroutine write_visu_stats_thermo(tm, dm)
    use udf_type_mod
    use precision_mod
    use operations
    implicit none 
    
    type(t_domain), intent(in) :: dm
    type(t_thermo), intent(in) :: tm

    integer :: iter 
    character(120) :: visuname
    
!==========================================================================================================
! write time averaged 3d data
!==========================================================================================================
    iter = tm%iteration
    visuname = 'time_averaged_thermo'
!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_HEADER, iter)
!----------------------------------------------------------------------------------------------------------
! write data, 
!----------------------------------------------------------------------------------------------------------
    call write_visu_field(dm, tm%t_mean,  dm%dccc, "time_averaged_T",  trim(visuname), SCALAR, CELL, iter)
    call write_visu_field(dm, tm%tt_mean, dm%dccc, "time_averaged_TT", trim(visuname), SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(visuname), XDMF_FOOTER, iter)
!==========================================================================================================
! write time averaged and space averaged 3d data (stored 2d or 1d data)
!==========================================================================================================
    if( ANY(dm%is_periodic(:))) then

    iter = tm%iteration
    visuname = 'time_space_averaged_thermo'
!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    if(.not. (dm%is_periodic(1) .and. dm%is_periodic(3))) &
    call write_visu_headerfooter(dm, trim(visuname), XDMF_HEADER, iter)
!----------------------------------------------------------------------------------------------------------
! write data, 
!----------------------------------------------------------------------------------------------------------
    call visu_average_periodic_data(tm%t_mean,   dm%dccc, dm, "time_space_averaged_T",   trim(visuname), iter)
    call visu_average_periodic_data(tm%tt_mean,  dm%dccc, dm, "time_space_averaged_TT",  trim(visuname), iter)
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    if(.not. (dm%is_periodic(1) .and. dm%is_periodic(3))) &
    call write_visu_headerfooter(dm, trim(visuname), XDMF_FOOTER, iter)

    end if
    
    return
  end subroutine
  
  !==========================================================================================================
  subroutine write_visu_any3darray(var, varname, visuname, dtmp, dm, iter)
    use udf_type_mod
    use precision_mod
    use operations
    use decomp_operation_mod
    use io_files_mod
    use parameters_constant_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    type(DECOMP_INFO), intent(in) :: dtmp
    character(*), intent(in) :: varname
    real(WP), dimension( dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3) ), intent(in) :: var
    integer, intent(in) :: iter 
    character(*), intent(in) :: visuname

    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil

    character(120) :: keyword

!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    keyword = trim(visuname)//"_"//trim(varname)//'_visu'
    call write_visu_headerfooter(dm, trim(keyword), XDMF_HEADER, iter)
    !if(nrank==0) write(*,*) keyword, iter
!----------------------------------------------------------------------------------------------------------
! write data, temperature, to cell centre
!----------------------------------------------------------------------------------------------------------
    if (is_same_decomp(dtmp, dm%dccc)) then
      call write_visu_field(dm, var, dm%dccc, trim(varname), trim(keyword), SCALAR, CELL, iter)

    else if (is_same_decomp(dtmp, dm%dpcc)) then
      call Get_x_midp_P2C_3D(var, accc, dm, dm%iAccuracy, dm%ibcx_qx(:), dm%fbcx_qx)
      call write_visu_field(dm, accc, dm%dccc, trim(varname), trim(keyword), SCALAR, CELL, iter)

    else if (is_same_decomp(dtmp, dm%dcpc)) then
      call transpose_x_to_y(var, acpc_ypencil, dm%dcpc)
      call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy(:), dm%fbcy_qy)
      call transpose_y_to_x(accc_ypencil, accc, dm%dccc)
      call write_visu_field(dm, accc, dm%dccc, trim(varname), trim(keyword), SCALAR, CELL, iter)

    else if (is_same_decomp(dtmp, dm%dccp)) then
      call transpose_x_to_y(var, accp_ypencil, dm%dccp)
      call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
      call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz(:), dm%fbcz_qz)
      call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
      call transpose_y_to_x(accc_ypencil, accc, dm%dccc)
      call write_visu_field(dm, accc, dm%dccc, trim(varname), trim(keyword), SCALAR, CELL, iter)

    else
      call Print_error_msg ("Given decomp_into is not supported "//trim(varname))
    end if
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    call write_visu_headerfooter(dm, trim(keyword), XDMF_FOOTER, iter)
    
    return
  end subroutine


!==========================================================================================================
!==========================================================================================================
  subroutine visu_average_periodic_data(data_in, dtmp, dm, str1, str2, iter)
    use udf_type_mod
    use parameters_constant_mod
    implicit none
    type(DECOMP_INFO), intent(in) :: dtmp
    type(t_domain), intent(in) :: dm
    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)), intent(in)  :: data_in
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    integer, intent(in) :: iter
    

    real(WP), dimension(dtmp%xsz(1), dtmp%xsz(2), dtmp%xsz(3)) :: a_xpencil
    real(WP), dimension(dtmp%ysz(1), dtmp%ysz(2), dtmp%ysz(3)) :: a_ypencil
    real(WP), dimension(dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3)) :: a_zpencil
    real(WP), dimension(dtmp%zsz(1), dtmp%zsz(2), dtmp%zsz(3)) :: b_zpencil

    real(WP), dimension( dtmp%ysz(2)) :: var

    integer :: i, j, k
    real(WP) :: sum

    if(dm%is_periodic(1) .and. &
       dm%is_periodic(3) .and. &
       dm%is_periodic(2)) then

    ! do nothing here, but bulk value output

    else if(dm%is_periodic(1) .and. &
            dm%is_periodic(3) .and. &
      .not. dm%is_periodic(2)) then

      do j = 1, dtmp%xsz(2)
        do k = 1, dtmp%xsz(3)
          sum =  ZERO
          do i = 1, dtmp%xsz(1)
            sum = sum + data_in(i, j, k) 
          end do
          sum =  sum/real(dtmp%xsz(1), WP)
          a_xpencil(:, j, k) = sum
        end do
      end do

      call transpose_x_to_y(a_xpencil, a_ypencil, dtmp)
      call transpose_y_to_z(a_ypencil, a_zpencil, dtmp)
      do i = 1, dtmp%zsz(1)
        do j = 1, dtmp%zsz(2)
          sum =  ZERO
          do k = 1, dtmp%zsz(3)
            sum = sum + a_zpencil(i, j, k) 
          end do
          b_zpencil(i, j, :) =  sum/real(dtmp%zsz(3), WP)
        end do
      end do
      call transpose_z_to_y(b_zpencil, a_ypencil, dtmp)
      var = a_ypencil(1,:,1)
      call write_visu_profile(dm, var, dm%dccc, trim(str1), trim(str2), SCALAR, CELL, 2, iter)
  
    else if(dm%is_periodic(1) .and. &
      .not. dm%is_periodic(3) .and. &
      .not. dm%is_periodic(2)) then

      do j = 1, dtmp%xsz(2)
        do k = 1, dtmp%xsz(3)
          sum =  ZERO
          do i = 1, dtmp%xsz(1)
            sum = sum + data_in(i, j, k) 
          end do
          sum =  sum/real(dtmp%xsz(1), WP)
          a_xpencil(:, j, k) = sum
        end do
      end do

      call write_visu_field(dm, a_xpencil, dm%dccc, trim(str1), trim(str2), SCALAR, CELL, iter)

    else if( &
      .not. dm%is_periodic(1) .and. &
            dm%is_periodic(3) .and. &
      .not. dm%is_periodic(2)) then

      call transpose_x_to_y(data_in,   a_ypencil, dtmp)
      call transpose_y_to_z(a_ypencil, a_zpencil, dtmp)
      do j = 1, dtmp%zsz(2)
        do i = 1, dtmp%zsz(1)
          sum =  ZERO
          do k = 1, dtmp%zsz(3)
            sum = sum + a_zpencil(i, j, k) 
          end do
          sum =  sum /real(dtmp%zsz(3), WP)
          b_zpencil(i, j, :) = sum
        end do
      end do

      call transpose_z_to_y(b_zpencil, a_ypencil, dtmp)
      call transpose_y_to_x(a_ypencil, a_xpencil, dtmp)

      call write_visu_field(dm, a_xpencil, dm%dccc, trim(str1), trim(str2), SCALAR, CELL, iter)

    else
      ! do nothing here
      !data_out = data_in
    end if


    return
  end subroutine



end module