!-------------------------------------------------------------------------------
!                      CHAPSim version 2.0.0
!                      --------------------------
! This file is part of CHAPSim, a general-purpose CFD tool.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 3 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------
!===============================================================================
!> \file input_general.f90
!>
!> \brief Reading the input parameters from the given file.
!>
!===============================================================================
module input_general_mod
  use precision_mod
  use parameters_constant_mod, only: ZERO
  implicit none

  character(len = 9), parameter :: INPUT_FILE = 'input.ini'

  integer, parameter :: ICASE_CHANNEL = 1, &
                        ICASE_PIPE    = 2, &
                        ICASE_ANNUAL  = 3, &
                        ICASE_TGV     = 4, &
                        ICASE_SINETEST= 5

  integer, parameter :: ICARTESIAN   = 1, &
                        ICYLINDRICAL = 2
                        
  integer, parameter :: ISTRET_NO     = 0, &
                        ISTRET_2SIDES = 1, &
                        ISTRET_BOTTOM = 2, &
                        ISTRET_TOP    = 3, &
                        ISTRET_CENTRE = 4

  integer, parameter :: ITIME_RK3    = 3, &
                        ITIME_RK3_CN = 2, &
                        ITIME_AB2    = 1

  integer, parameter :: IBC_INTERIOR    = 0, &
                        IBC_PERIODIC    = 1, &
                        IBC_UDIRICHLET  = 2, &
                        IBC_SYMMETRIC   = 3, &
                        IBC_ASYMMETRIC  = 4

!                        IBC_INLET_MEAN  = 4, &
!                        IBC_INLET_TG    = 5, &
!                        IBC_INLET_MAP   = 6, &
!                        IBC_INLET_DB    = 7, &
!                        IBC_OUTLET_EXPO = 8, &
!                        IBC_OUTLET_CONV = 9, &
!                        IBC_INTERIOR    = 0, &
                        
  integer, parameter :: IACCU_CD2 = 1, &
                        IACCU_CD4 = 2, &
                        IACCU_CP4 = 3, &
                        IACCU_CP6 = 4

  integer, parameter :: NDIM = 3

  integer, parameter :: INITIAL_RANDOM  = 1, &
                        INITIAL_RESTART = 2, &
                        INITIAL_INTERPL = 3

  integer, parameter :: IVIS_EXPLICIT   = 1, &
                        IVIS_SEMIMPLT   = 2

  integer, parameter :: IDRVF_NO        = 0, &
                        IDRVF_MASSFLUX  = 1, &
                        IDRVF_SKINFRIC  = 2, &
                        IDRVF_PRESLOSS  = 3


  ! flow type
  integer :: icase
  integer :: ithermo
  integer :: icht

  ! domain decomposition
  integer :: p_row
  integer :: p_col

  ! domain geometry
  real(WP) :: lxx, lzz, lyt, lyb

  ! domain mesh
  integer :: ncx, ncy, ncz
  integer :: istret
  real(WP) :: rstret

  ! flow parameter
  real(WP) :: ren

  ! time stepping
  real(WP) :: dt
  integer :: nIterFlow0
  integer :: nIterFlow1

  ! boundary condition
  integer :: ifbcx(1:2)
  integer :: ifbcy(1:2)
  integer :: ifbcz(1:2)
  real(WP) :: uxinf(2)
  real(WP) :: uyinf(2)
  real(WP) :: uzinf(2)


  ! InOutParam
  integer :: irestart
  integer :: ncheckpoint
  integer :: nvisu
  integer :: iterStatsFirst
  integer :: nstats

  ! NumOption
  integer :: iAccuracy
  integer :: iTimeScheme
  integer :: iviscous
  integer :: ipressure

  ! initial fields
  real(WP) :: renIni
  integer  :: nIterIniRen
  real(WP) :: initNoise

  ! PeriodicDrv
  integer :: idriven
  real(WP) :: drvf
  
  ! ThermoParam
  integer :: ifluid
  integer :: igravity
  real(WP) :: lenRef
  real(WP) :: t0Ref
  real(WP) :: tiRef
  integer :: itbcy(1:2)
  real(WP) :: tbcy(1:2)
  integer :: nIterThermo0
  integer :: nIterThermo1

  ! parameters from restart
  integer :: iterchkpt = 0       ! iteration number from restart/checkpoint
  real(WP) :: tThermo  = ZERO
  real(WP) :: tFlow    = ZERO

  ! derived parameters
  logical :: is_periodic(3)
  integer :: icoordinate

  integer :: nsubitr
  real(WP) :: tGamma(0 : 3)
  real(WP) :: tZeta (0 : 3)
  real(WP) :: tAlpha(0 : 3)

  ! procedure
  public  :: Initialize_general_input
  private :: Set_periodic_bc
  private :: Set_timestepping_coefficients
  
  
contains
!===============================================================================
!===============================================================================
!> \brief Reading the input parameters from the given file. The file name could
!> be changed in the above module.     
!>
!> This subroutine is called at beginning of solver.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
  subroutine Initialize_general_input ()
!===============================================================================
! Module files
!===============================================================================
    use iso_fortran_env, only : ERROR_UNIT, IOSTAT_END
    use parameters_constant_mod, only: ZERO, ONE, TWO, PI
    use mpi_mod, only: ncol, nrow
    implicit none
!===============================================================================
! Local arguments
!===============================================================================
    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit

    character(len = 80) :: section_name
    character(len = 80) :: variableName
    integer :: slen
!===============================================================================
    open ( newunit = inputUnit, &
           file    = INPUT_FILE, &
           status  = 'old', &
           action  = 'read', &
           iostat  = ioerr, &
           iomsg   = iotxt )
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', INPUT_FILE, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 1
    end if
    
    do 
      read(inputUnit, '(a)', iostat = ioerr) section_name
      slen = len_trim(section_name)
      if (ioerr /=0 ) exit
      if ( (section_name(1:1) == ';') .or. &
          (section_name(1:1) == '#') .or. &
          (section_name(1:1) == ' ') .or. &
          (slen == 0) ) then
        cycle
      end if

      block_section: if ( section_name(1:slen) == '[mpi]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, nrow
        read(inputUnit, *, iostat = ioerr) variableName, ncol

      else if ( section_name(1:slen) == '[flowtype]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, icase
        read(inputUnit, *, iostat = ioerr) variableName, ithermo
        read(inputUnit, *, iostat = ioerr) variableName, icht

      else if ( section_name(1:slen) == '[decomposition]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, p_row
        read(inputUnit, *, iostat = ioerr) variableName, p_col

      else if ( section_name(1:slen) == '[geometry]' )  then 

        read(inputUnit, *, iostat = ioerr) variableName, lxx
        read(inputUnit, *, iostat = ioerr) variableName, lyt
        read(inputUnit, *, iostat = ioerr) variableName, lyb
        read(inputUnit, *, iostat = ioerr) variableName, lzz

      else if ( section_name(1:slen) == '[mesh]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ncx
        read(inputUnit, *, iostat = ioerr) variableName, ncy
        read(inputUnit, *, iostat = ioerr) variableName, ncz
        read(inputUnit, *, iostat = ioerr) variableName, istret
        read(inputUnit, *, iostat = ioerr) variableName, rstret

      else if ( section_name(1:slen) == '[flowparams]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ren

      else if ( section_name(1:slen) == '[timestepping]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, dt
        read(inputUnit, *, iostat = ioerr) variableName, nIterFlow0
        read(inputUnit, *, iostat = ioerr) variableName, nIterFlow1

      else if ( section_name(1:slen) == '[boundary]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, &
            ifbcx(1), ifbcx(2), uxinf(1), uxinf(2)
        read(inputUnit, *, iostat = ioerr) variableName, &
            ifbcy(1), ifbcy(2), uyinf(1), uyinf(2)
        read(inputUnit, *, iostat = ioerr) variableName, &
            ifbcz(1), ifbcz(2), uzinf(1), uzinf(2)

      else if ( section_name(1:slen) == '[ioparams]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, irestart
        read(inputUnit, *, iostat = ioerr) variableName, ncheckpoint
        read(inputUnit, *, iostat = ioerr) variableName, nvisu
        read(inputUnit, *, iostat = ioerr) variableName, iterStatsFirst
        read(inputUnit, *, iostat = ioerr) variableName, nstats

      else if ( section_name(1:slen) == '[schemes]' ) then
        read(inputUnit, *, iostat = ioerr) variableName, iAccuracy
        read(inputUnit, *, iostat = ioerr) variableName, iTimeScheme
        read(inputUnit, *, iostat = ioerr) variableName, iviscous
        read(inputUnit, *, iostat = ioerr) variableName, ipressure

      else if ( section_name(1:slen) == '[initialization]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, renIni
        read(inputUnit, *, iostat = ioerr) variableName, nIterIniRen
        read(inputUnit, *, iostat = ioerr) variableName, initNoise

      else if ( section_name(1:slen) == '[periodicdriven]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, idriven
        read(inputUnit, *, iostat = ioerr) variableName, drvf

      else if ( section_name(1:slen) == '[thermohydraulics]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ifluid
        read(inputUnit, *, iostat = ioerr) variableName, igravity
        read(inputUnit, *, iostat = ioerr) variableName, lenRef
        read(inputUnit, *, iostat = ioerr) variableName, t0Ref
        read(inputUnit, *, iostat = ioerr) variableName, tiRef
        read(inputUnit, *, iostat = ioerr) variableName, itbcy(1), itbcy(2)
        read(inputUnit, *, iostat = ioerr) variableName, tbcy(1), tbcy(2)
        read(inputUnit, *, iostat = ioerr) variableName, nIterThermo0
        read(inputUnit, *, iostat = ioerr) variableName, nIterThermo1

      else
        exit
      end if block_section
    end do

    if(ioerr /= IOSTAT_END) &
    call Print_error_msg( 'Problem reading '//INPUT_FILE // &
    'in Subroutine: '// "Initialize_general_input")

    close(inputUnit)

    ! set up some default values to overcome wrong input
    if (icase == ICASE_CHANNEL) then
      if(istret /= ISTRET_2SIDES) &
      call Print_warning_msg ("Grids are not two-side clustered.")
      lyb = - ONE
      lyt = ONE
    else if (icase == ICASE_PIPE) then
      if(istret /= ISTRET_TOP)    &
      call Print_warning_msg ("Grids are not near-wall clustered.")
      lyb = ZERO
      lyt = ONE
    else if (icase == ICASE_ANNUAL) then
      if(istret /= ISTRET_2SIDES) &
      call Print_warning_msg ("Grids are not two-side clustered.")
      lyt = ONE
    else if (icase == ICASE_TGV) then
      if(istret /= ISTRET_NO) &
      call Print_warning_msg ("Grids are clustered.")
      lxx = TWO * PI
      lzz = TWO * PI
      lyt = PI
      lyb = -PI
    else if (icase == ICASE_SINETEST) then
      if(istret /= ISTRET_NO) &
      call Print_warning_msg ("Grids are clustered.")
      lxx = TWO * PI
      lzz = TWO * PI
      lyt = PI
      lyb = -PI
    else 
      ! do nothing...
    end if

    ! to set up cooridnates
    if (icase == ICASE_CHANNEL) then
      icoordinate = ICARTESIAN
    else if (icase == ICASE_PIPE) then
      icoordinate = ICYLINDRICAL
      ifbcy(1) = IBC_INTERIOR
    else if (icase == ICASE_ANNUAL) then
      icoordinate = ICYLINDRICAL
    else if (icase == ICASE_TGV) then
      icoordinate = ICARTESIAN
    else 
      icoordinate = ICARTESIAN
    end if

    ! to set up periodic boundary conditions
    is_periodic(:) = .false.
    call Set_periodic_bc ( ifbcx, is_periodic(1) )
    call Set_periodic_bc ( ifbcy, is_periodic(2) )
    call Set_periodic_bc ( ifbcz, is_periodic(3) )

    ! to set up parameters for time stepping
    call Set_timestepping_coefficients ( )

    return
  end subroutine Initialize_general_input
!===============================================================================
!===============================================================================
!> \brief Periodic B.C. configuration if one side of periodic bc is detected.     
!>
!> This subroutine is locally called once by \ref Initialize_general_input.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  bc            boundary condition index
!> \param[out]    flg           logical flag for periodic b.c.
!_______________________________________________________________________________
  subroutine Set_periodic_bc( bc, flg )
    integer, intent(inout) :: bc(1:2)
    logical, intent(out) :: flg

    if ( (bc(1) == IBC_PERIODIC) .or. (bc(2) == IBC_PERIODIC) ) then
      bc(1) = IBC_PERIODIC
      bc(2) = IBC_PERIODIC
      flg = .true.
    else 
      flg = .false.
    end if

  end subroutine Set_periodic_bc
!===============================================================================
!===============================================================================
!> \brief Define parameters for time stepping.     
!>
!> This subroutine is locally called once by \ref Initialize_general_input.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
  subroutine Set_timestepping_coefficients()
    use parameters_constant_mod
    implicit none

    if(iTimeScheme == ITIME_RK3 .or. &
       iTimeScheme == ITIME_RK3_CN) then
      
      nsubitr = 3
      tGamma(0) = ONE
      tGamma(1) = EIGHT / FIFTEEN
      tGamma(2) = FIVE / TWELVE
      tGamma(3) = THREE / FOUR

      tZeta (0) = ZERO
      tZeta (1) = ZERO
      tZeta (2) = -SEVENTEEN / SIXTY
      tZeta (3) = -FIVE / TWELVE

    else if (iTimeScheme == ITIME_AB2) then

      nsubitr = 1
      tGamma(0) = ONE
      tGamma(1) = THREE / TWO
      tGamma(2) = ZERO
      tGamma(3) = ZERO

      tZeta (0) = ZERO
      tZeta (1) = -HALF
      tZeta (2) = ZERO
      tZeta (3) = ZERO

    else 

      nsubitr = 0
      tGamma(:) = ZERO
      tZeta (:) = ZERO

    end if 
    
    tAlpha(:) = tGamma(:) + tZeta(:)

  end subroutine Set_timestepping_coefficients

end module 



