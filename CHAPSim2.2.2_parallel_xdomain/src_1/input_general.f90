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
  use parameters_constant_mod
  use type_vars_mod
  implicit none
!-------------------------------------------------------------------------------
! global variable, time related
!-------------------------------------------------------------------------------
  integer :: ndomain
  logical :: is_any_energyeq

  integer :: ifluid
  real(WP) :: lenRef
  real(WP) :: T0Ref

  integer  :: niter

  integer  :: nsubitr

  integer  :: iTimeScheme
  real(WP) :: tGamma(0 : 3)
  real(WP) :: tZeta (0 : 3)
  real(WP) :: tAlpha(0 : 3)

  real(WP) :: sigma1p
  real(WP) :: sigma2p

!-------------------------------------------------------------------------------
! procedure
!-------------------------------------------------------------------------------
  public  :: Read_general_input_common
  public  :: Read_general_input_subdomain
  private :: Set_periodic_bc
  private :: Set_timestepping_coefficients
  
  
contains
!===============================================================================
!===============================================================================
!> \brief Reading the input parameters from the given file. The file name could
!> be changed in the above module.     
!>
!> This subroutine is called at beginning of solver.
!> [mpi] all ranks
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
  subroutine Read_general_input_common (flname)
!===============================================================================
! Module files
!===============================================================================
    use iso_fortran_env,         only : ERROR_UNIT, IOSTAT_END
    use parameters_constant_mod, only : ZERO, ONE, TWO, PI
    use mpi_mod
    implicit none
    character(len = 11), intent(in) :: flname

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit

    character(len = 80) :: section_name
    character(len = 80) :: variableName
    character(len = 20) :: formati='(2X, A32, I20.1)'
    character(len = 20) :: format2i='(2X, A32, 2I10.1)'
    character(len = 20) :: formatr='(2X, A32, F20.4)'
    character(len = 20) :: formate='(2X, A32, E20.4)'
    character(len = 20) :: format2r='(2X, A32, 2F10.2)'
    character(len = 20) :: formatir='(2X, A32, I10.1, F10.2)'
    character(len = 20) :: format2i2r='(2X, A32, 2I10.1, 2F10.2)'
    character(len = 20) :: formats='(2X, A32, A20)'
    integer :: slen

    integer :: icase
    integer  :: ibcy(1:5, 1:2) ! bc, type: 1-5: u,v,w,p,T;  1:2=yst, yen
    integer  :: ibcz(1:5, 1:2) ! bc, type: 1-5: u,v,w,p,T;  1:2=zst, zen
    real(WP) :: fbcy(1:5, 1:2)
    real(WP) :: fbcz(1:5, 1:2)
    real(WP) :: lzz, lyt, lyb
    integer  :: ncy, ncz
    integer  :: istret
    real(WP) :: rstret
    integer :: icoordinate

    if(nrank == 0) call Print_debug_start_msg("CHAPSim2.0 Starts ...")

!-------------------------------------------------------------------------------
! open file
!-------------------------------------------------------------------------------
    open ( newunit = inputUnit, &
           file    = flname, &
           status  = 'old', &
           action  = 'read', &
           iostat  = ioerr, &
           iomsg   = iotxt )
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', flname, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 1
    end if

    if(nrank == 0) call Print_debug_start_msg("Reading General Parameters from "//flname//" ...")
    
    do 
!-------------------------------------------------------------------------------
! reading headings/comments
!-------------------------------------------------------------------------------
      read(inputUnit, '(a)', iostat = ioerr) section_name
      slen = len_trim(section_name)
      if (ioerr /=0 ) exit
      if ( (section_name(1:1) == ';') .or. &
           (section_name(1:1) == '#') .or. &
           (section_name(1:1) == ' ') .or. &
           (slen == 0) ) then
        cycle
      end if
      if(nrank == 0) call Print_debug_start_msg("Reading "//section_name(1:slen))
!-------------------------------------------------------------------------------
! [decomposition]
!-------------------------------------------------------------------------------
      block_section: if ( section_name(1:slen) == '[decomposition]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ndomain
        read(inputUnit, *, iostat = ioerr) variableName, nrow
        read(inputUnit, *, iostat = ioerr) variableName, ncol
        if(nrank == 0) then
          write(*, formati) ' x-direction domain number                  :', ndomain
          write(*, formati) ' y-direction domain number (default Row)    :', nrow
          write(*, formati) ' z-direction domain number (default Column) :', ncol
        end if
!-------------------------------------------------------------------------------
! [flow]
!-------------------------------------------------------------------------------
      else if ( section_name(1:slen) == '[flow]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, icase
        read(inputUnit, *, iostat = ioerr) variableName, is_any_energyeq
        read(inputUnit, *, iostat = ioerr) variableName, ifluid
        read(inputUnit, *, iostat = ioerr) variableName, lenRef
        read(inputUnit, *, iostat = ioerr) variableName, t0Ref
        
        if(nrank == 0) then
          if(icase == ICASE_CHANNEL) write(*, formats) ' Case : ', "Channel flow" 
          if(icase == ICASE_PIPE)    write(*, formats) ' Case : ', "Pipe flow"
          if(icase == ICASE_ANNUAL)  write(*, formats) ' Case : ', "Annual flow"
          if(icase == ICASE_TGV2D)   write(*, formats) ' Case : ', "Taylor Green Vortex flow (2D)"
          if(icase == ICASE_TGV3D)   write(*, formats) ' Case : ', "Taylor Green Vortex flow (3D)"
        end if
!-------------------------------------------------------------------------------
! [boundary]
!-------------------------------------------------------------------------------
      else if ( section_name(1:slen) == '[boundary]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ibcy(1:5, 1), fbcy(1:5, 1)
        read(inputUnit, *, iostat = ioerr) variableName, ibcy(1:5, 2), fbcy(1:5, 2)
        read(inputUnit, *, iostat = ioerr) variableName, ibcz(1:5, 1), fbcz(1:5, 1)
        read(inputUnit, *, iostat = ioerr) variableName, ibcz(1:5, 2), fbcz(1:5, 2)
!-------------------------------------------------------------------------------
! [geometry]
!-------------------------------------------------------------------------------
      else if ( section_name(1:slen) == '[geometry]' )  then 

        read(inputUnit, *, iostat = ioerr) variableName, lyt
        read(inputUnit, *, iostat = ioerr) variableName, lyb
        read(inputUnit, *, iostat = ioerr) variableName, lzz
!-------------------------------------------------------------------------------
! [mesh]
!-------------------------------------------------------------------------------
      else if ( section_name(1:slen) == '[mesh]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ncy
        read(inputUnit, *, iostat = ioerr) variableName, ncz
        read(inputUnit, *, iostat = ioerr) variableName, istret
        read(inputUnit, *, iostat = ioerr) variableName, rstret
!-------------------------------------------------------------------------------
! [timestepping]
!-------------------------------------------------------------------------------
      else if ( section_name(1:slen) == '[timestepping]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, dt
        read(inputUnit, *, iostat = ioerr) variableName, iTimeScheme
        if(nrank == 0) write(*, formate) ' Physical Time Step : ', dt

      else
        exit
      end if block_section
    end do
!-------------------------------------------------------------------------------
! end of reading
!-------------------------------------------------------------------------------
    if(ioerr /= IOSTAT_END) &
    call Print_error_msg( 'Problem reading '//flname // &
    'in Subroutine: '// "Read_general_input")

    close(inputUnit)
!-------------------------------------------------------------------------------
! allocate udf 
!-------------------------------------------------------------------------------
    allocate( domain (ndomain) )
    allocate(   flow (ndomain) )
    allocate( thermo (ndomain) )
!-------------------------------------------------------------------------------
! assign
!-------------------------------------------------------------------------------
    if (icase == ICASE_CHANNEL) then
      icoordinate = ICARTESIAN
    else if (icase == ICASE_PIPE) then
      icoordinate = ICYLINDRICAL
      ifbcy(:, 1) = IBC_INTERIOR
    else if (icase == ICASE_ANNUAL) then
      icoordinate = ICYLINDRICAL
    else if (icase == ICASE_TGV2D) then
      icoordinate = ICARTESIAN
    else if (icase == ICASE_TGV3D) then
      icoordinate = ICARTESIAN
    else 
      icoordinate = ICARTESIAN
    end if

    call Set_periodic_bc ( ifbcy )
    call Set_periodic_bc ( ifbcz )

    if (icase == ICASE_CHANNEL) then
      if(istret /= ISTRET_2SIDES .and. nrank == 0) &
      call Print_warning_msg ("Grids are not two-side clustered.")
      lyb = - ONE
      lyt = ONE
    else if (icase == ICASE_PIPE) then
      if(istret /= ISTRET_TOP .and. nrank == 0) &
      call Print_warning_msg ("Grids are not near-wall clustered.")
      lyb = ZERO
      lyt = ONE
    else if (icase == ICASE_ANNUAL) then
      if(istret /= ISTRET_2SIDES .and. nrank == 0 ) &
      call Print_warning_msg ("Grids are not two-side clustered.")
      lyt = ONE
    else if (icase == ICASE_TGV2D) then
      if(istret /= ISTRET_NO .and. nrank == 0) &
      call Print_warning_msg ("Grids are clustered.")
      lxx = TWO * PI
      lzz = TWO * PI
      lyt =   PI
      lyb = - PI
    else if (icase == ICASE_TGV3D) then
      if(istret /= ISTRET_NO .and. nrank == 0) &
      call Print_warning_msg ("Grids are clustered.")
      lxx = TWO * PI
      lzz = TWO * PI
      lyt =   PI
      lyb = - PI
    else if (icase == ICASE_SINETEST) then
      if(istret /= ISTRET_NO .and. nrank == 0) &
      call Print_warning_msg ("Grids are clustered.")
      lxx = TWO * PI
      lzz = TWO * PI
      lyt =   PI
      lyb = - PI
    else 
      ! do nothing...
    end if

    domain(:)%ndom        = ndomain
    domain(:)%icase       = icase
    domain(:)%icoordinate = icoordinate

    domain(:)%ibcy(:, :)  = ibcy(:, :)
    domain(:)%ibcz(:, :)  = ibcz(:, :)
    domain(:)%fbcy(:, :)  = fbcy(:, :)
    domain(:)%fbcz(:, :)  = fbcz(:, :)

    domain(:)%lzz         = lzz
    domain(:)%lyt         = lyt
    domain(:)%lyb         = lyb
    domain(:)%nc(2)       = ncy
    domain(:)%nc(3)       = ncz

    domain(:)%istret      = istret
    domain(:)%rstret      = rstret

    flow(:)%dt            = dt
!-------------------------------------------------------------------------------
! time and scheme information is global
!-------------------------------------------------------------------------------
    call Set_timestepping_coefficients ()
!-------------------------------------------------------------------------------
! builidng up information
!-------------------------------------------------------------------------------
    domain%is_periodic(2) == .false.
    domain%is_periodic(3) == .false.
    if(domain%ibcy(1, 1) == IBC_PERIODIC) domain%is_periodic(2) = .true.
    if(domain%ibcz(1, 1) == IBC_PERIODIC) domain%is_periodic(3) = .true.

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine Read_general_input_common
!===============================================================================
!===============================================================================
!> \brief Reading the input parameters from the given file. The file name could
!> be changed in the above module.     
!>
!> This subroutine is called at beginning of solver.
!> [mpi] all ranks
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
  subroutine Read_general_input_subdomain (flname, fl, th, dm)
!===============================================================================
! Module files
!===============================================================================
    use iso_fortran_env,         only : ERROR_UNIT, IOSTAT_END
    use parameters_constant_mod, only : ZERO, ONE, TWO, PI
    use mpi_mod
    implicit none
    character(len = 11), intent(in) :: flname
    type(t_flow),   intent(inout)   :: fl
    type(t_thermo), intent(inout)   :: th
    type(t_domain), intent(inout)   :: dm

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit

    character(len = 80) :: section_name
    character(len = 80) :: variableName
    character(len = 20) :: formati='(2X, A32, I20.1)'
    character(len = 20) :: format2i='(2X, A32, 2I10.1)'
    character(len = 20) :: formatr='(2X, A32, F20.4)'
    character(len = 20) :: formate='(2X, A32, E20.4)'
    character(len = 20) :: format2r='(2X, A32, 2F10.2)'
    character(len = 20) :: formatir='(2X, A32, I10.1, F10.2)'
    character(len = 20) :: format2i2r='(2X, A32, 2I10.1, 2F10.2)'
    character(len = 20) :: formats='(2X, A32, A20)'
    integer :: slen

    integer  :: ibcx(1:5, 1:2) ! bc, type: 1-5: u,v,w,p,T;  1:2=xst, xen
    real(WP) :: fbcx(1:5, 1:2) ! bc, values 
    real(WP) :: lxx
    integer  :: ncx
    
    integer :: igravity
    
    real(WP) :: TiRef

!-------------------------------------------------------------------------------
! open file
!-------------------------------------------------------------------------------
    open ( newunit = inputUnit, &
            file    = flname, &
            status  = 'old', &
            action  = 'read', &
            iostat  = ioerr, &
            iomsg   = iotxt )
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', flname, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 1
    end if

    if(nrank == 0) call Print_debug_start_msg("Reading Parameters from "//flname//" ...")
!-------------------------------------------------------------------------------
! reading headings/comments
!-------------------------------------------------------------------------------
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
      if(nrank == 0) call Print_debug_start_msg("Reading "//section_name(1:slen))

!-------------------------------------------------------------------------------
! [subdomain]
!-------------------------------------------------------------------------------
      block_section: if ( section_name(1:slen) == '[subdomain]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, dm%idom

        if(nrank == 0) then
          write(*, formati) ' Current x-subdomain index :', dm%idom
        end if
!-------------------------------------------------------------------------------
! [flowtype]
!-------------------------------------------------------------------------------
      else if ( section_name(1:slen) == '[flowtype]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, dm%ithermo
        read(inputUnit, *, iostat = ioerr) variableName, dm%icht
        read(inputUnit, *, iostat = ioerr) variableName, fl%ren
        read(inputUnit, *, iostat = ioerr) variableName, fl%idriven
        read(inputUnit, *, iostat = ioerr) variableName, fl%drvf
        
        
        read(inputUnit, *, iostat = ioerr) variableName, th%igravity
        read(inputUnit, *, iostat = ioerr) variableName, th%tiRef

        if(nrank == 0) then
          if(dm%ithermo == 0) write(*, formats) ' Thermal field : ', 'No' 
          if(dm%ithermo == 1) write(*, formats) ' Thermal field : ', 'Yes' 
          if(dm%icht    == 0) write(*, formats) ' Conjugate Heat Transfer : ', 'No' 
          if(dm%icht    == 1) write(*, formats) ' Conjugate Heat Transfer : ', 'Yes'

          write(*, formatr) ' Reynolds number : ', fl%ren

          if(fl%idriven == IDRVF_MASSFLUX) then
            write(*, formats) ' Flow driven by : ', "Constant Mass Flux" 
          else if(fl%idriven == IDRVF_SKINFRIC) then
            write(*, formats) ' Flow driven by : ', "Provided Skin Friction" 
            write(*, formatr) ' Skin Friction : ', fl%drvf
          else if(fl%idriven == IDRVF_PRESLOSS) then
            write(*, formats) ' Flow driven by : ', "Pressure loss" 
            write(*, formatr) ' pressure loss per length : ', fl%drvf
          else 
            write(*, formats) ' No external driven force' 
          end if

          if(dm%ithermo /= 0) then
            th%ifluid   = ifluid
            th%lenRef   = lenRef
            th%T0Ref    = T0Ref
            if(nrank == 0)
              write(*, formati) ' Fluid type : ', th%ifluid
              write(*, formati) ' Gravity force direction : ', th%igravity
              write(*, formatr) ' Reference length for normalisation : ', th%lenRef
              write(*, formatr) ' Reference temperature for normalisation : ', th%t0Ref
              write(*, formatr) ' Initialisation temperature : ', th%tiRef
            end if
          end if

      else if ( section_name(1:slen) == '[geometry]' )  then 
!-------------------------------------------------------------------------------
! [geometry]
!-------------------------------------------------------------------------------
        read(inputUnit, *, iostat = ioerr) variableName, dm%lxx
        read(inputUnit, *, iostat = ioerr) variableName, dm%nc(1)

        if (dm%icase == ICASE_TGV2D .or. dm%icase == ICASE_TGV3D) then
          if(dm%istret /= ISTRET_NO .and. nrank == 0) &
          call Print_warning_msg ("Grids are clustered.")
          dm%lxx = TWO * PI
        end if
        
        if(nrank == 0) then
          write(*, formatr) ' length in x : ', dm%lxx
          write(*, formatr) ' length in z : ', dm%lzz
          write(*, formatr) ' length in y : ', dm%lyt - dm%lyb
          write(*, formatr) ' bottom in y : ', dm%lyb
          write(*, formatr) '    top in y : ', dm%lyt

          write(*, formati) ' Mesh Cell Number in x : ', dm%nc(1)
          write(*, formati) ' Mesh Cell Number in y : ', dm%nc(2)
          write(*, formati) ' Mesh Cell Number in z : ', dm%nc(3)
          if(dm%istret == ISTRET_NO) write(*, formats) ' Y mesh stretching : ', 'No' 
          if(dm%istret /= ISTRET_NO) then
            write(*, formats) ' Y mesh stretching : ', 'Yes' 
            write(*, formatr) ' Stretching factor beta : ', dm%rstret
          end if
        end if
        
      else if ( section_name(1:slen) == '[boundary]' ) then
!-------------------------------------------------------------------------------
! [boundary]
!-------------------------------------------------------------------------------
        read(inputUnit, *, iostat = ioerr) variableName, ibcx(1:5, 1), fbcx(1:5, 1)
        read(inputUnit, *, iostat = ioerr) variableName, ibcx(1:5, 2), fbcx(1:5, 2)

        call Set_periodic_bc ( ifbcx )
        dm%ibcx(:, :) = ibcx(:, :)
        dm%fbcx(:, :) = fbcx(:, :)

        if(nrank == 0) then
          write(*, formats) ' BC in x for variable: ', "u, v, w, p, T"
          do i = 1, 5
           write(*, format2i2r) ' BC for  in x-start side : ', dm%ibcx(i, 1), dm%fbcx(i, 1)
           write(*, format2i2r) ' BC for  in x-end   side : ', dm%ibcx(i, 2), dm%fbcx(i, 2)
          end do

          write(*, formats) ' BC in y for variable: ', "u, v, w, p, T"
          do i = 1, 5
           write(*, format2i2r) ' BC for  in y-start side : ', dm%ibcy(i, 1), dm%fbcy(i, 1)
           write(*, format2i2r) ' BC for  in y-end   side : ', dm%ibcy(i, 2), dm%fbcy(i, 2)
          end do

          write(*, formats) ' BC in z for variable: ', "u, v, w, p, T"
          do i = 1, 5
           write(*, format2i2r) ' BC for  in z-start side : ', dm%ibcz(i, 1), dm%fbcz(i, 1)
           write(*, format2i2r) ' BC for  in z-end   side : ', dm%ibcz(i, 2), dm%fbcz(i, 2)
          end do
        end if

      else if ( section_name(1:slen) == '[initialization]' ) then
!-------------------------------------------------------------------------------
! [initialization]
!-------------------------------------------------------------------------------
        read(inputUnit, *, iostat = ioerr) variableName, fl%irestart
        read(inputUnit, *, iostat = ioerr) variableName, fl%nrsttckpt
        read(inputUnit, *, iostat = ioerr) variableName, fl%renIni
        read(inputUnit, *, iostat = ioerr) variableName, fl%nIterIniRen
        read(inputUnit, *, iostat = ioerr) variableName, fl%initNoise

        if(nrank == 0) then
          fl%nrsttckpt = 0
          if(irestart == 0) then
            write(*, formats) ' Start from : ', 'Scratch'
            write(*, formatr) ' Initial Reynolds No : ', fl%renIni
            write(*, formati) ' Initial Re lasts until : ', fl%nIterIniRen
            write(*, formatr) ' Initial velocity perturbation : ', fl%initNoise
          else 
            write(*, formats) ' Start from : ', 'Restart' 
            write(*, formati) ' Restart iteration : ', fl%nrsttckpt
          end if
        end if

      else if ( section_name(1:slen) == '[schemes]' ) then
!-------------------------------------------------------------------------------
! [schemes]
!-------------------------------------------------------------------------------
        read(inputUnit, *, iostat = ioerr) variableName, dm%iAccuracy
        read(inputUnit, *, iostat = ioerr) variableName, dm%iviscous

        if(nrank == 0) then
          write(*, formati) ' Spatial Accuracy Order: ', dm%iAccuracy
          if(dm%iviscous == IVIS_EXPLICIT) write(*, formats) ' Viscous Term : ', 'Explicit Scheme'
          if(dm%iviscous == IVIS_SEMIMPLT) write(*, formats) ' Viscous Term : ', 'Semi-implicit Scheme'
        end if

      else if ( section_name(1:slen) == '[simcontrol]' ) then
!-------------------------------------------------------------------------------
! [simcontrol]
!-------------------------------------------------------------------------------
        read(inputUnit, *, iostat = ioerr) variableName, fl%nIterFlowStart
        read(inputUnit, *, iostat = ioerr) variableName, fl%nIterFlowEnd
        read(inputUnit, *, iostat = ioerr) variableName, th%nIterThermoStart
        read(inputUnit, *, iostat = ioerr) variableName, th%nIterThermoEnd

        if(nrank == 0) then
          write(*, format2i) ' Flow simulation lasts between : ', fl%nIterFlowStart, fl%nIterFlowEnd
          if(dm%ithermo == 1) write(*, format2i) ' Heat simulation lasts between : ', th%nIterThermoStart, th%nIterThermoEnd
        end if

      else if ( section_name(1:slen) == '[ioparams]' ) then
!-------------------------------------------------------------------------------
! [ioparams]
!-------------------------------------------------------------------------------
        read(inputUnit, *, iostat = ioerr) variableName, dm%nfreqckpt
        read(inputUnit, *, iostat = ioerr) variableName, dm%nvisu
        read(inputUnit, *, iostat = ioerr) variableName, dm%nIterStatsStart
        read(inputUnit, *, iostat = ioerr) variableName, dm%nfreqStats
        if(nrank == 0) then
          write(*, formati) ' Raw  data written at every : ', dm%nfreqckpt
          write(*, formati) ' Vis  data written at every : ', dm%nvisu
          write(*, formati) ' Stat data written from     : ', dm%nIterStatsStart
          write(*, formati) ' Stat data written at every : ', dm%nfreqStats
        end if                        
      else
        exit
      end if block_section
    end do

    if(ioerr /= IOSTAT_END) &
    call Print_error_msg( 'Problem reading '//flname // &
    'in Subroutine: '// "Read_general_input")

    close(inputUnit)

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine Read_general_input_subdomain
!===============================================================================
!===============================================================================
!> \brief Periodic B.C. configuration if one side of periodic bc is detected.     
!>
!> This subroutine is locally called once by \ref Read_general_input.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  bc            boundary condition index
!> \param[out]    flg           logical flag for periodic b.c.
!_______________________________________________________________________________
  subroutine Set_periodic_bc( bc )
    integer, intent(inout) :: bc(1:5,1:2)

    do i = 1, 5
      if ( (bc(i, 1) == IBC_PERIODIC) .or. (bc(i, 2) == IBC_PERIODIC) ) then
        bc(i, 1) = IBC_PERIODIC
        bc(i, 2) = IBC_PERIODIC
      end if
    end do

    return
  end subroutine Set_periodic_bc
!===============================================================================
!===============================================================================
!> \brief Define parameters for time stepping.     
!>
!> This subroutine is locally called once by \ref Read_general_input.
!> [mpi] all ranks
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     none          NA
!> \param[out]    none          NA
!_______________________________________________________________________________
  subroutine Set_timestepping_coefficients( )
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