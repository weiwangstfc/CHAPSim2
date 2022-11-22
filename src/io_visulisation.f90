module io_visulisation_mod
  use io_tools_mod
  implicit none 

  integer, parameter :: XDMF_HEADER = 1, &
                        XDMF_FOOTER = 2
  integer, parameter :: PLANE_AVERAGE = -1
  character(6), parameter :: SCALAR = "Scalar", &
                             VECTOR = "Vector", &
                             TENSOR = "Tensor"
  character(4), parameter :: CELL = "Cell", &
                             NODE = "Node"
  character(len=*), parameter :: io_name = "solution-io"

  integer, allocatable :: nnd_visu(:, :)
  integer, allocatable :: ncl_visu(:, :)

  character(6)  :: svisudim

  private :: write_snapshot_headerfooter
  private :: write_field

  public  :: write_snapshot_ini
  public  :: write_snapshot_flow


contains
!==========================================================================================================
! xszV means:
! x - xpencil, could also be y, z
! sz - size, could also be st, en
! V - visulisation
! xszV is the same as dppp%xsz/nskip, based on nodes, considering the skip nodes
!==========================================================================================================
  subroutine write_snapshot_ini(dm)
    use udf_type_mod
    use parameters_constant_mod, only: MAXP
    use decomp_2d, only: xszV, yszV, zszV
    use files_io_mod
    implicit none 
    type(t_domain), intent(in) :: dm

    integer :: i, j ,k
    real(WP), allocatable :: xp(:), yp(:), zp(:)
    character(120):: grid_flname
    character(120):: keyword
    integer :: iogrid

!----------------------------------------------------------------------------------------------------------
! allocate
!----------------------------------------------------------------------------------------------------------
    allocate (nnd_visu(3, nxdomain))
    allocate (ncl_visu(3, nxdomain))
    nnd_visu = 0
    ncl_visu = 0
!----------------------------------------------------------------------------------------------------------
! global size
!----------------------------------------------------------------------------------------------------------
    if(dm%visu_idim == Ivisudim_3D) then
      svisudim = "3d"
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
    else if(dm%visu_idim == Ivisudim_2D_Xa) then
      svisudim = "2d_xa"
      nnd_visu(1, dm%idom) = 1
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
    else if(dm%visu_idim == Ivisudim_2D_Ya) then
      svisudim = "2d_ya"
      nnd_visu(1, dm%idom) = xszV(1)
      nnd_visu(2, dm%idom) = 1
      nnd_visu(3, dm%idom) = zszV(3)
      do i = 1, 3
        if(dm%is_periodic(i)) then 
          ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
          nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
        else 
          ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
        end if
      end do
    else if(dm%visu_idim == Ivisudim_2D_Za) then
      svisudim = "2d_za"
      nnd_visu(1, dm%idom) = xszV(1)
      nnd_visu(2, dm%idom) = yszV(2)
      nnd_visu(3, dm%idom) = 1
      do i = 1, 3
        if(dm%is_periodic(i)) then 
          ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
          nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
        else 
          ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
        end if
      end do
    else if(dm%visu_idim == Ivisudim_1D_XZa) then
      svisudim = "2d_xza"
      nnd_visu(1, dm%idom) = 1
      nnd_visu(2, dm%idom) = yszV(2)
      nnd_visu(3, dm%idom) = 1
      do i = 1, 3
        if(dm%is_periodic(i)) then 
          ncl_visu(i, dm%idom) = nnd_visu(i, dm%idom)
          nnd_visu(i, dm%idom) = nnd_visu(i, dm%idom) + 1
        else 
          ncl_visu(i, dm%idom) = MAX(nnd_visu(i, dm%idom) - 1, 1)
        end if
      end do
    else
      call Print_error_msg("Input visu_idim does not support.") 
    end if
!----------------------------------------------------------------------------------------------------------
! write grids
!----------------------------------------------------------------------------------------------------------    
    if(nrank == 0) then

      allocate ( xp(nnd_visu(1, dm%idom)) )
      allocate ( yp(nnd_visu(2, dm%idom)) )
      allocate ( zp(nnd_visu(3, dm%idom)) )

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
        zp(k) = real(k-1, WP) * dm%h(2) * dm%visu_nskip(3)
      enddo


      keyword = trim(svisudim)//"_grid_x"
      call generate_pathfile_name(grid_flname, dm%idom, keyword, dir_visu, 'xml')
      open(newunit = iogrid, file = trim(grid_flname), action = "write", status="replace")
      write(iogrid, *) xp
      close(iogrid)

      keyword = trim(svisudim)//"_grid_y"
      call generate_pathfile_name(grid_flname, dm%idom, keyword, dir_visu, 'xml')
      open(newunit = iogrid, file = trim(grid_flname), action = "write", status="replace")
      write(iogrid, *) yp
      close(iogrid)

      keyword = trim(svisudim)//"_grid_z"
      call generate_pathfile_name(grid_flname, dm%idom, keyword, dir_visu, 'xml')
      open(newunit = iogrid, file = trim(grid_flname), action = "write", status="replace")
      write(iogrid, *) zp
      close(iogrid)

      deallocate(xp, yp, zp)
    end if

    return
  end subroutine

!==========================================================================================================
! ref: https://www.xdmf.org/index.php/XDMF_Model_and_Format
!==========================================================================================================
  subroutine write_snapshot_headerfooter(dm, iheadfoot, iter)
    use precision_mod
    use parameters_constant_mod, only: MAXP
    use udf_type_mod, only: t_domain
    use decomp_2d, only: nrank, mytype, xszV, yszV, zszV
    use files_io_mod
    implicit none 
    integer, intent(in)        :: iheadfoot
    integer, intent(in)        :: iter
    type(t_domain), intent(in) :: dm

    character(120):: keyword
    character(120):: visu_flname
    character(120):: grid_flname(3)
    character(1)  :: str(3)

    integer :: ioxdmf
    logical :: file_exists
    real(WP), allocatable :: xp(:), yp(:), zp(:)
    
    
    integer :: i, j, k
    if(nrank /= 0) return
!----------------------------------------------------------------------------------------------------------
! file name
!----------------------------------------------------------------------------------------------------------
    keyword = "snapshot_"//trim(svisudim)
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
      str(1) = "x"
      str(2) = "y"
      str(3) = "z"
      do i = 1, 3
        keyword = trim(svisudim)//"_grid_"//trim(str(i))
        call generate_file_name(grid_flname(i), dm%idom, keyword, 'xml')
        INQUIRE(FILE = trim(trim(dir_visu)//"/"//trim(grid_flname(i))), exist = file_exists)
        if(.not.file_exists) then
          call Print_error_msg("Mesh file for visu does not exist. Filename = "//trim(trim(dir_visu)//"/"//trim(grid_flname(i))))
        end if
      end do

      allocate ( xp(nnd_visu(1, dm%idom)) )
      allocate ( yp(nnd_visu(2, dm%idom)) )
      allocate ( zp(nnd_visu(3, dm%idom)) )

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
        zp(k) = real(k-1, WP) * dm%h(2) * dm%visu_nskip(3)
      enddo

!----------------------------------------------------------------------------------------------------------
! write header
! geometry is based on node coordinates
! to do: write mesh into mesh.bin 
!----------------------------------------------------------------------------------------------------------
      write(ioxdmf, '(A)')'<?xml version="1.0" ?>'
      write(ioxdmf, '(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(ioxdmf, '(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      write(ioxdmf, *)' <Domain>'
      write(ioxdmf, *)'   <Grid Name="'//trim(svisudim)//'" GridType="Uniform">'
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
      deallocate(xp, yp, zp)
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
  subroutine write_field(dm, var, dtmp, varname, attributetype, centring, iter)
    use precision_mod
    use decomp_2d
    use decomp_2d_io
    use udf_type_mod, only: t_domain
    use files_io_mod
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP), intent(in) :: var(:, :, :)
    character(len=*), intent(in) :: varname
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

!----------------------------------------------------------------------------------------------------------
! file name
!----------------------------------------------------------------------------------------------------------
    keyword = "snapshot_"//trim(svisudim)
    call generate_pathfile_name(visu_flname_path, dm%idom, keyword, dir_visu, 'xdmf', iter)
!----------------------------------------------------------------------------------------------------------
! write data into binary file
!----------------------------------------------------------------------------------------------------------
    if(dm%visu_idim == Ivisudim_3D) then
      keyword = trim(svisudim)//"_"//trim(varname)
      call generate_pathfile_name(data_flname_path, dm%idom, keyword, dir_data, 'bin', iter)
      call decomp_2d_write_one(X_PENCIL, var, trim(data_flname_path), dtmp)
    else if(dm%visu_idim == Ivisudim_1D_XZa) then
      !to add 1D profile
    else 
      keyword = trim(svisudim)//"_"//trim(varname)
      call generate_file_name(data_flname, dm%idom, keyword, 'bin', iter)
      call generate_pathfile_name(data_flname_path, dm%idom, keyword, dir_data, 'bin', iter)
      call decomp_2d_write_plane(X_PENCIL, var, dm%visu_idim, PLANE_AVERAGE, dir_data, trim(data_flname), io_name, dtmp)
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
  subroutine write_snapshot_flow(dm, fl)
    use udf_type_mod
    use precision_mod
    use operations
    implicit none 
    type(t_domain), intent(in) :: dm
    type(t_flow), intent(in) :: fl

    integer :: iter 
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
    real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
    real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
    real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil


    iter = fl%iteration
!----------------------------------------------------------------------------------------------------------
! write xdmf header
!----------------------------------------------------------------------------------------------------------
    call write_snapshot_headerfooter(dm, XDMF_HEADER, iter)
!----------------------------------------------------------------------------------------------------------
! write data, pressure, to cell centre
!----------------------------------------------------------------------------------------------------------
    call write_field(dm, fl%pres, dm%dccc, "pr", SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! qx, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
    call Get_x_midp_P2C_3D(fl%qx, accc, dm, dm%ibcx(:, 1))
    call write_field(dm, accc, dm%dccc, "ux", SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! qy, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%ibcy(:, 2))
    call transpose_y_to_x(accc_ypencil, accc, dm%dccc)
    call write_field(dm, accc, dm%dccc, "uy", SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! qz, default x-pencil, staggered to cell centre
!----------------------------------------------------------------------------------------------------------
    call transpose_x_to_y(fl%qz, accp_ypencil, dm%dccp)
    call transpose_y_to_z(accp_ypencil, accp_zpencil, dm%dccp)
    call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%ibcz(:, 3))
    call transpose_z_to_y(accc_zpencil, accc_ypencil, dm%dccc)
    call transpose_y_to_z(accc_ypencil, accc, dm%dccc)
    call write_field(dm, accc, dm%dccc, "uz", SCALAR, CELL, iter)
!----------------------------------------------------------------------------------------------------------
! write xdmf footer
!----------------------------------------------------------------------------------------------------------
    call write_snapshot_headerfooter(dm, XDMF_FOOTER, iter)
    
    return
  end subroutine



end module