module decomp_extended_mod
  use parameters_constant_mod
  implicit none


  public :: ypencil_index_lgl2ggl
  public :: zpencil_index_llg2ggg
  public :: zpencil_index_ggg2llg

  contains
!==========================================================================================================
  subroutine ypencil_index_lgl2ggl(vin, vou, dtmp)
    use decomp_2d
    implicit none

    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(dtmp%ysz(1),               dtmp%ysz(2), dtmp%ysz(3)), intent(in)  :: vin
    real(WP), dimension(dtmp%yst(1) : dtmp%yen(2), dtmp%ysz(2), dtmp%zsz(3)), intent(out) :: vou

    integer :: i, j, k, ii
    vou = ZERO
    do k = 1, dtmp%ysz(3)
      do j = 1, dtmp%ysz(2)
        do i = 1, dtmp%ysz(1)
          ii = dtmp%yst(1) + i - 1
          vou(ii, j, k) = vin(i, j, k)
        end do
      end do
    end do
    return
  end subroutine 
!==========================================================================================================
  subroutine zpencil_index_llg2ggg(vin, vou, dtmp)
    use decomp_2d
    use index_mod
    implicit none

    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(dtmp%zsz(1),               dtmp%zsz(2),               dtmp%zsz(3)), intent(in)  :: vin
    real(WP), dimension(dtmp%zst(1) : dtmp%zen(1), dtmp%zst(2) : dtmp%zen(2), dtmp%zsz(3)), intent(out) :: vou

    integer :: i, j, k, jj, ii

    vou = ZERO
    do k = 1, dtmp%zsz(3)
      do j = 1, dtmp%zsz(2)
        jj = dtmp%zst(2) + j - 1 !local2global_yid(j, dtmp)
        do i = 1, dtmp%zsz(1)
          ii = dtmp%zst(1) + i - 1
          vou(ii, jj, k) = vin(i, j, k)
        end do
      end do
    end do
    return
  end subroutine
!==========================================================================================================  
  subroutine zpencil_index_ggg2llg(vin, vou, dtmp)
    use decomp_2d
    use index_mod
    implicit none

    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(dtmp%zst(1) : dtmp%zen(1), dtmp%zst(2) : dtmp%zen(2), dtmp%zsz(3)), intent(in)   :: vin
    real(WP), dimension(dtmp%zsz(1),               dtmp%zsz(2),               dtmp%zsz(3)), intent(out)  :: vou
    

    integer :: i, j, k, jj, ii
!write(*,*) 'vin', nrank, size(vin, 1), size(vin, 2),size(vin, 3)
!write(*,*) 'vou', nrank, size(vou, 1), size(vou, 2),size(vou, 3)

    vou = ZERO
    do k = 1, dtmp%zsz(3)
      do j = 1, dtmp%zsz(2)
        jj = dtmp%zst(2) + j - 1 !local2global_yid(j, dtmp)
        do i = 1, dtmp%zsz(1)
          ii = dtmp%zst(1) + i - 1
          vou(i, j, k) = vin(ii, jj, k)
        end do
      end do
    end do
    return
  end subroutine
end module 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================================================================================================
module poisson_interface_mod
  use parameters_constant_mod
  use fft2decomp_interface_mod
  use decomp_2d_poisson
  use decomp_extended_mod
  use fishpack_fft
  implicit none

  public :: initialise_fft
  public :: solve_fft_poisson

contains
!==========================================================================================================
!==========================================================================================================
  subroutine initialise_fft(dm)
    use udf_type_mod
    implicit none 
    type(t_domain), intent(in) :: dm

    if(fft_lib == FFT_2DECOMP ) then 
      call build_up_fft2decomp_interface(dm)
      call decomp_2d_poisson_init()
    else if(fft_lib == FFT_FISHPACK) then 
      call fishpack_fft_init(dm)
    else 
      write(*, *) 'Error in selecting FFT libs'
      STOP
    end if
  return 
  end subroutine 
!==========================================================================================================
!==========================================================================================================
  subroutine solve_fft_poisson(rhs_xpencil, dm)
    use udf_type_mod
    implicit none 
    type(t_domain), intent(in) :: dm
    integer :: i, j, k
    real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ), intent(INOUT) :: rhs_xpencil
    real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: rhs_ypencil
    real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: rhs_zpencil
    real(WP), dimension( dm%dccc%zst(1) : dm%dccc%zen(1), &
                         dm%dccc%zst(2) : dm%dccc%zen(2), &
                         dm%dccc%zst(3) : dm%dccc%zen(3) ) :: rhs_zpencil_ggg

    if(fft_lib == FFT_2DECOMP ) then 
      call transpose_x_to_y (rhs_xpencil, rhs_ypencil, dm%dccc)
      call transpose_y_to_z (rhs_ypencil, rhs_zpencil, dm%dccc)
      call zpencil_index_llg2ggg(rhs_zpencil, rhs_zpencil_ggg, dm%dccc)

      call poisson(rhs_zpencil_ggg)

      call zpencil_index_ggg2llg(rhs_zpencil_ggg, rhs_zpencil, dm%dccc)
      call transpose_z_to_y (rhs_zpencil, rhs_ypencil, dm%dccc)
      call transpose_y_to_x (rhs_ypencil, rhs_xpencil, dm%dccc)
    else if(fft_lib == FFT_FISHPACK) then 
      call fishpack_fft_simple(rhs_xpencil, dm)
    else 
      write(*, *) 'Error in selecting FFT libs'
      STOP
    end if

    write(*, '(3I3, E13.5)') (i, j, k, rhs_xpencil(i, j, k), i = 1, dm%dccc%xsz(1), j = 1, dm%dccc%xsz(2), k = 1, dm%dccc%xsz(3))
    
  return 
  end subroutine 


end module 
