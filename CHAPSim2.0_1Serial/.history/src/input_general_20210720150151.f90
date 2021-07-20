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
  use precision_mod, only : WP
  use parameters_constant_mod, only : ZERO
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

  integer, parameter :: IBC_INTERIOR    = 9, &
                        IBC_PERIODIC    = 1, &
                        IBC_UDIRICHLET  = 2, &
                        IBC_SYMMETRIC   = 3, &
                        IBC_ASYMMETRIC  = 4, &
                        IBC_NEUMANN     = 5, &
                        IBC_CONVECTIVE  = 6, &
                        IBC_TURBGEN     = 7, &
                        IBC_DATABASE    = 8
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

  integer, parameter :: INITIAL_RANDOM  = 0, &
                        INITIAL_RESTART = 1, &
                        INITIAL_INTERPL = 2

  integer, parameter :: IVIS_EXPLICIT   = 1, &
                        IVIS_SEMIMPLT   = 2

  integer, parameter :: IDRVF_NO        = 0, &
                        IDRVF_MASSFLUX  = 1, &
                        IDRVF_SKINFRIC  = 2, &
                        IDRVF_PRESLOSS  = 3

  integer, parameter :: THERMAL_BC_CONST_T  = 0, &
                        THERMAL_BC_CONST_H  = 1

!-------------------------------------------------------------------------------
! domain decomposition
!-------------------------------------------------------------------------------
  integer :: p_row
  integer :: p_col
!-------------------------------------------------------------------------------
! flow type
!-------------------------------------------------------------------------------
  integer :: icase
  integer :: ithermo
  integer :: icht
!-------------------------------------------------------------------------------
! boundary condition
!-------------------------------------------------------------------------------
  integer  :: ifbcx(1:2)
  integer  :: ifbcy(1:2)
  integer  :: ifbcz(1:2)
  real(WP) :: uxinf(2)
  real(WP) :: uyinf(2)
  real(WP) :: uzinf(2)
!-------------------------------------------------------------------------------
! flow parameter
!-------------------------------------------------------------------------------
  real(WP) :: ren
  integer  :: idriven
  real(WP) :: drvf
!-------------------------------------------------------------------------------
! domain geometry
!-------------------------------------------------------------------------------
  real(WP) :: lxx, lzz, lyt, lyb
!-------------------------------------------------------------------------------
! domain mesh
!-------------------------------------------------------------------------------
  integer  :: ncx, ncy, ncz
  integer  :: istret
  real(WP) :: rstret
!-------------------------------------------------------------------------------
! initialization
!-------------------------------------------------------------------------------
  integer :: irestart
  integer :: nrsttckpt
  real(WP) :: renIni
  integer  :: nIterIniRen
  real(WP) :: initNoise
!-------------------------------------------------------------------------------
! time stepping
!-------------------------------------------------------------------------------
  real(WP) :: dt
!-------------------------------------------------------------------------------
! schemes
!-------------------------------------------------------------------------------
  integer :: iAccuracy
  integer :: iTimeScheme
  integer :: iviscous
!-------------------------------------------------------------------------------
! simulation control
!-------------------------------------------------------------------------------
  integer  :: nIterFlowStart
  integer  :: nIterFlowEnd
  integer  :: nIterThermoStart
  integer  :: nIterThermoEnd
!-------------------------------------------------------------------------------
! InOutParam
!-------------------------------------------------------------------------------
  integer :: nfreqckpt
  integer :: nvisu
  integer :: nIterStatsStart
  integer :: nfreqStats
!-------------------------------------------------------------------------------
! ThermoParam
!-------------------------------------------------------------------------------
  integer  :: ifluid
  integer  :: igravity
  real(WP) :: lenRef
  real(WP) :: t0Ref
  real(WP) :: tiRef
  integer  :: itbcy(1:2)
  real(WP) :: tbcy(1:2)

!-------------------------------------------------------------------------------
! ThermoParam
!-------------------------------------------------------------------------------
  ! parameters from restart
  real(WP) :: tThermo  = ZERO
  real(WP) :: tFlow    = ZERO

  ! derived parameters
  logical :: is_periodic(3)
  integer :: icoordinate

  integer  :: nsubitr
  real(WP) :: tGamma(0 : 3)
  real(WP) :: tZeta (0 : 3)
  real(WP) :: tAlpha(0 : 3)

  real(WP) :: sigma1p
  real(WP) :: sigma2p

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
    use iso_fortran_env,         only : ERROR_UNIT, IOSTAT_END
    use parameters_constant_mod, only : ZERO, ONE, TWO, PI
    use mpi_mod,                 only : ncol, nrow
    implicit none
!===============================================================================
! Local arguments
!===============================================================================
    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit

    character(len = 80) :: section_name
    character(len = 80) :: variableName
    character(len = 20) :: formati='(2X, A32, I20.1)'
    character(len = 20) :: format2i='(2X, A32, 2I10.1)'
    character(len = 20) :: formatr='(2X, A32, F20.4)'
    character(len = 20) :: formate='(2X, A32, E20.4)'
    character(len = 20) :: format2r='(2X, A32, 2F10.4)'
    character(len = 20) :: formats='(2X, A32, A20)'
    integer :: slen
!===============================================================================
    call Print_debug_start_msg("CHAPSim2.0 Starts ...")


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

    call Print_debug_start_msg("Reading General Parameters from "//INPUT_FILE//" ...")
    
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
      call Print_debug_start_msg("Reading "//section_name(1:slen))

      block_section: if ( section_name(1:slen) == '[flowtype]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, icase
        read(inputUnit, *, iostat = ioerr) variableName, ithermo
        read(inputUnit, *, iostat = ioerr) variableName, icht

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

        if(icase == ICASE_CHANNEL) write(*, formats) ' Case : ', "Channel flow" 
        if(icase == ICASE_PIPE)    write(*, formats) ' Case : ', "Pipe flow"
        if(icase == ICASE_ANNUAL)  write(*, formats) ' Case : ', "Annual flow"
        if(icase == ICASE_TGV)     write(*, formats) ' Case : ', "Taylor Green Vortex flow"
        if(ithermo == 0)           write(*, formats) ' Thermal field : ', 'No' 
        if(ithermo == 1)           write(*, formats) ' Thermal field : ', 'Yes' 
        if(icht    == 0)           write(*, formats) ' Conjugate Heat Transfer : ', 'No' 
        if(icht    == 1)           write(*, formats) ' Conjugate Heat Transfer : ', 'Yes'

      else if ( section_name(1:slen) == '[boundary]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, &
            ifbcx(1), ifbcx(2), uxinf(1), uxinf(2)
        read(inputUnit, *, iostat = ioerr) variableName, &
            ifbcy(1), ifbcy(2), uyinf(1), uyinf(2)
        read(inputUnit, *, iostat = ioerr) variableName, &
            ifbcz(1), ifbcz(2), uzinf(1), uzinf(2)

        ! to set up periodic boundary conditions
        is_periodic(:) = .false.
        call Set_periodic_bc ( ifbcx, is_periodic(1) )
        call Set_periodic_bc ( ifbcy, is_periodic(2) )
        call Set_periodic_bc ( ifbcz, is_periodic(3) )

        if(is_periodic(1)) then
          write(*, formats) ' BC in x : ', "Periodic" 
        else
          write(*, format2i) ' BC in x : ', ifbcx(1:2)
          write(*, format2r) ' xBC input : ', uxinf(1:2)
        end if

        if(is_periodic(2)) then
          write(*, formats) ' BC in y : ', "Periodic" 
        else
          write(*, format2i) ' BC in y : ', ifbcy(1:2)
          write(*, format2r) ' yBC input : ', uyinf(1:2)
        end if

        if(is_periodic(3)) then
          write(*, formats) ' BC in z : ', "Periodic" 
        else
          write(*, format2i) ' BC in z : ', ifbcz(1:2)
          write(*, format2r) ' zBC input : ', uzinf(1:2)
        end if

      else if ( section_name(1:slen) == '[flowparams]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ren
        read(inputUnit, *, iostat = ioerr) variableName, idriven
        read(inputUnit, *, iostat = ioerr) variableName, drvf

        write(*, formatr) ' Reynolds number : ', ren
        if(is_periodic(1)) then
          if(idriven == IDRVF_MASSFLUX) write(*, formats) ' Flow driven by : ', "Constant Mass Flux" 
          if(idriven == IDRVF_SKINFRIC) then
            write(*, formats) ' Flow driven by : ', "Provided Skin Friction" 
            write(*, formatr) ' Skin Friction : ', drvf
          end if
        end if

      else if ( section_name(1:slen) == '[decomposition]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, p_row
        read(inputUnit, *, iostat = ioerr) variableName, p_col
        write(*, formati) ' Processor in Row :',     nrow
        write(*, formati) ' Processor in Column :',  ncol

      else if ( section_name(1:slen) == '[geometry]' )  then 

        read(inputUnit, *, iostat = ioerr) variableName, lxx
        read(inputUnit, *, iostat = ioerr) variableName, lyt
        read(inputUnit, *, iostat = ioerr) variableName, lyb
        read(inputUnit, *, iostat = ioerr) variableName, lzz

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
          lyt =   PI
          lyb = - PI
        else if (icase == ICASE_SINETEST) then
          if(istret /= ISTRET_NO) &
          call Print_warning_msg ("Grids are clustered.")
          lxx = TWO * PI
          lzz = TWO * PI
          lyt =   PI
          lyb = - PI
        else 
          ! do nothing...
        end if

        write(*, formatr) ' length in x : ', lxx
        write(*, formatr) ' length in z : ', lzz
        write(*, formatr) ' length in y : ', lyt - lyb
        write(*, formatr) ' bottom in y : ', lyb
        write(*, formatr) '    top in y : ', lyt

      else if ( section_name(1:slen) == '[mesh]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ncx
        read(inputUnit, *, iostat = ioerr) variableName, ncy
        read(inputUnit, *, iostat = ioerr) variableName, ncz
        read(inputUnit, *, iostat = ioerr) variableName, istret
        read(inputUnit, *, iostat = ioerr) variableName, rstret

        write(*, formati) ' Mesh Cell Number in x : ', ncx
        write(*, formati) ' Mesh Cell Number in y : ', ncy
        write(*, formati) ' Mesh Cell Number in z : ', ncz
        if(istret == ISTRET_NO) write(*, formats) ' Y mesh stretching : ', 'No' 
        if(istret /= ISTRET_NO) then
          write(*, formats) ' Y mesh stretching : ', 'Yes' 
          write(*, formatr) ' Stretching factor beta : ', rstret
        end if

      else if ( section_name(1:slen) == '[timestepping]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, dt
        write(*, formate) ' Physical Time Step : ', dt
      
      else if ( section_name(1:slen) == '[initialization]' ) then
        read(inputUnit, *, iostat = ioerr) variableName, irestart
        read(inputUnit, *, iostat = ioerr) variableName, nrsttckpt
        read(inputUnit, *, iostat = ioerr) variableName, nfreqckpt
        read(inputUnit, *, iostat = ioerr) variableName, renIni
        read(inputUnit, *, iostat = ioerr) variableName, nIterIniRen
        read(inputUnit, *, iostat = ioerr) variableName, initNoise

        if(irestart == 0) then
          nrsttckpt = 0
          write(*, formats) ' Start from : ', 'Scratch'
          write(*, formatr) ' Starting Reynolds No : ', renIni
          write(*, formati) ' Starting Re lasts until : ', nIterIniRen
          write(*, formatr) ' Starting velocity perturbation : ', initNoise
        else 
          write(*, formats) ' Start from : ', 'Restart' 
          write(*, formati) ' Restart iteration : ', nrsttckpt
        end if

      else if ( section_name(1:slen) == '[simcontrol]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, nIterFlowStart
        read(inputUnit, *, iostat = ioerr) variableName, nIterFlowEnd
        read(inputUnit, *, iostat = ioerr) variableName, nIterThermoStart
        read(inputUnit, *, iostat = ioerr) variableName, nIterThermoEnd

        write(*, format2i) ' Flow simulation lasts between : ', nIterFlowStart, nIterFlowEnd
        if(ithermo == 1) write(*, format2i) ' Heat simulation lasts between : ', nIterThermoStart, nIterThermoEnd

      else if ( section_name(1:slen) == '[ioparams]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, nfreqckpt
        read(inputUnit, *, iostat = ioerr) variableName, nvisu
        read(inputUnit, *, iostat = ioerr) variableName, nIterStatsStart
        read(inputUnit, *, iostat = ioerr) variableName, nfreqStats

        write(*, formati) ' Raw  data written at every : ', nfreqckpt
        write(*, formati) ' Visualisation data written at every : ', nvisu
        write(*, formati) ' Statistic Data written from : ', nIterStatsStart
        write(*, formati) ' Statistic Data written at every : ', nfreqStats

      else if ( section_name(1:slen) == '[schemes]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, iAccuracy
        read(inputUnit, *, iostat = ioerr) variableName, iTimeScheme
        read(inputUnit, *, iostat = ioerr) variableName, iviscous

        write(*, formati) ' Spatial Accuracy Order: ', iAccuracy
        write(*, formatr) ' Time Stepping Scheme: ', dt
        if(iviscous == IVIS_EXPLICIT) write(*, formats) ' Viscous Term : ', 'Explicit Scheme'
        if(iviscous == IVIS_SEMIMPLT) write(*, formats) ' Viscous Term : ', 'Semi-implicit Scheme'

       
      else if ( section_name(1:slen) == '[thermohydraulics]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ifluid
        read(inputUnit, *, iostat = ioerr) variableName, igravity
        read(inputUnit, *, iostat = ioerr) variableName, lenRef
        read(inputUnit, *, iostat = ioerr) variableName, t0Ref
        read(inputUnit, *, iostat = ioerr) variableName, tiRef
        read(inputUnit, *, iostat = ioerr) variableName, itbcy(1), itbcy(2)
        read(inputUnit, *, iostat = ioerr) variableName, tbcy(1), tbcy(2)
        
        write(*, formati) ' Fluid type : ', ifluid
        write(*, formati) ' Gravity force direction : ', igravity
        write(*, formatr) ' Reference length for normalisation : ', lenRef
        write(*, formatr) ' Reference temperature for normalisation : ', t0Ref
        write(*, formatr) ' Initialisation temperature : ', tiRef
        write(*, format2i) ' Thermal BC type: ', itbcy(1), itbcy(2)
        write(*, format2r) ' Thermal BC value: ', tbcy(1), tbcy(2)

      else
        exit
      end if block_section
    end do

    if(ioerr /= IOSTAT_END) &
    call Print_error_msg( 'Problem reading '//INPUT_FILE // &
    'in Subroutine: '// "Initialize_general_input")

    close(inputUnit)

    call Set_timestepping_coefficients ( )

    call Print_debug_end_msg

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
    logical, intent(out  ) :: flg

    if ( (bc(1) == IBC_PERIODIC) .or. (bc(2) == IBC_PERIODIC) ) then
      bc(1) = IBC_PERIODIC
      bc(2) = IBC_PERIODIC
      flg = .true.
    else 
      flg = .false.
    end if
    return
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

    !option 1: to set up pressure treatment, for O(dt^2)
    sigma1p = ONE
    sigma2p = HALF

    !option 2: to set up pressure treatment, for O(dt)
    !sigma1p = ONE
    !sigma2p = ONE

    if(iTimeScheme == ITIME_RK3     .or. &
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



