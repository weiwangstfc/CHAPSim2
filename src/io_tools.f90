module io_tools_mod
  implicit none

  !----------------------------------------------------------------------------------------------------------
  ! io parameters
  !----------------------------------------------------------------------------------------------------------
  integer, parameter :: Ivisudim_3D    = 0, &
                        Ivisudim_2D_Xa = 1, & ! x averaged, should not change this value. 
                        Ivisudim_2D_Ya = 2, & ! y averaged
                        Ivisudim_2D_Za = 3, & ! z averaged
                        Ivisudim_1D_XZa= 4    ! xz averaged
  integer, parameter :: X_PENCIL = 1, & ! x-pencil
                        Y_PENCIL = 2, & ! y-pencil
                        Z_PENCIL = 3    ! z-pencil

  public :: initialize_decomp_io
  public :: generate_file_name
  public :: generate_pathfile_name
  public :: mesh_output
  
contains
!==========================================================================================================
  subroutine initialize_decomp_io(dm)
    use udf_type_mod
    use decomp_2d_io
    implicit none 
    type(t_domain), intent(in) :: dm
    
    logical is_start1 ! is index starting from 1.
!----------------------------------------------------------------------------------------------------------
! if not #ifdef ADIOS2, do nothing below.
!----------------------------------------------------------------------------------------------------------
    call decomp_2d_io_init()
!----------------------------------------------------------------------------------------------------------
! re-define the grid mesh size, considering the nskip
! based on decomp_info of dppp (default one defined)
!---------------------------------------------------------------------------------------------------------- 
    is_start1 = .true.
    call init_coarser_mesh_statV(dm%visu_nskip(1), dm%visu_nskip(2), dm%visu_nskip(3), is_start1)
    call init_coarser_mesh_statS(dm%stat_nskip(1), dm%stat_nskip(2), dm%stat_nskip(3), is_start1)

  end subroutine 
!==========================================================================================================
  subroutine generate_pathfile_name(flname, dmtag, keyword, path, extension, timetag)
    use typeconvert_mod
    implicit none 
    integer, intent(in)      :: dmtag
    character(*), intent(in) :: keyword
    character(*), intent(in) :: path
    character(*), intent(in) :: extension
    character(120), intent(out) :: flname
    integer, intent(in), optional     :: timetag

    if(present(timetag)) then
      flname = trim(path)//"/domain"//trim(int2str(dmtag))//'_'//trim(keyword)//'_'//trim(int2str(timetag))//"."//trim(extension)
    else 
      flname = trim(path)//"/domain"//trim(int2str(dmtag))//'_'//trim(keyword)//"."//trim(extension)
    end if

    return
  end subroutine
!==========================================================================================================
  subroutine generate_file_name(flname, dmtag, keyword, extension, timetag)
    use typeconvert_mod
    implicit none 
    integer, intent(in)      :: dmtag
    
    character(*), intent(in) :: keyword
    character(*), intent(in) :: extension
    character(120), intent(out) :: flname
    integer, intent(in), optional      :: timetag

    if(present(timetag)) then
      flname = "domain"//trim(int2str(dmtag))//'_'//trim(keyword)//'_'//trim(int2str(timetag))//"."//trim(extension)
    else 
      flname = "domain"//trim(int2str(dmtag))//'_'//trim(keyword)//"."//trim(extension)
    end if

    
    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine mesh_output(dm)
    use udf_type_mod
    implicit none

    type(t_domain), intent(in) :: dm


    real(WP), allocatable :: id(:, :, :)

    allocate( id ( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) )
    id(:, :, :) = real(nrank, WP)
    call write_snapshot_any3darray(id, 'rank', 'mesh', dm%dpcc, dm, 0)
    deallocate(id)

    return
  end subroutine mesh_output
end module