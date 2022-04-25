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
  implicit none

  logical :: is_any_energyeq
  public  :: Read_input_parameters

contains
!===============================================================================
!> \brief Reading the input parameters from the given file.     
!> Scope:  mpi    called-freq    xdomain
!>         all    once           all
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[in]     none          NA
!> \param[out]    none          NA
!===============================================================================
  subroutine Read_input_parameters
    use wtformat_mod
    use mpi_mod
    use parameters_constant_mod
    use vars_df_mod
    use thermo_info_mod
    use boundary_conditions_mod
    implicit none

    character(len = 18) :: flname = 'input_chapsim.ini'
    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    integer  :: slen

    character(len = 80) :: secname
    character(len = 80) :: varname
    integer  :: itmp
    real(WP) :: rtmp
    real(WP), allocatable :: rtmpx(:)
    integer, allocatable  :: itmpx(:)
    integer :: i, j, k
    
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

    if(nrank == 0) &
    call Print_debug_start_msg("Reading General Parameters from "//flname//" ...")
!-------------------------------------------------------------------------------
! reading input
!-------------------------------------------------------------------------------
    do 
!-------------------------------------------------------------------------------
! reading headings/comments
!-------------------------------------------------------------------------------
      read(inputUnit, '(a)', iostat = ioerr) secname
      slen = len_trim(secname)
      if (ioerr /=0 ) exit
      if ( (secname(1:1) == ';') .or. &
           (secname(1:1) == '#') .or. &
           (secname(1:1) == ' ') .or. &
           (slen == 0) ) then
        cycle
      end if
      if(nrank == 0) call Print_debug_start_msg("Reading "//secname(1:slen))
!-------------------------------------------------------------------------------
! [decomposition]
!-------------------------------------------------------------------------------
      block_secname: if ( secname(1:slen) == '[decomposition]' ) then

        read(inputUnit, *, iostat = ioerr) varname, nxdomain
        read(inputUnit, *, iostat = ioerr) varname, nrow
        read(inputUnit, *, iostat = ioerr) varname, ncol
        allocate( domain (nxdomain) )
        allocate(   flow (nxdomain) )

        do i = 1, nxdomain
          domain(i)%idom = i
        end do

        if(nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  x-dir domain number             :', nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) '  y-dir domain number (mpi Row)   :', nrow
            write (OUTPUT_UNIT, wrtfmt1i) '  z-dir domain number (mpi Column):', ncol
          end do
        end if
!-------------------------------------------------------------------------------
! [flow]
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[domain]' ) then

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%icase
        domain(:)%icase = domain(1)%icase

        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%lxx

        read(inputUnit, *, iostat = ioerr) varname, rtmp
        domain(:)%lyt = rtmp

        read(inputUnit, *, iostat = ioerr) varname, rtmp
        domain(:)%lyb = rtmp

        read(inputUnit, *, iostat = ioerr) varname, rtmp
        domain(:)%lzz = rtmp

        if(nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  icase                  :', domain(i)%icase
          end do
        end if
!-------------------------------------------------------------------------------
! [boundary] 
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[boundary]' ) then

        do i = 1, nxdomain
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx(1:2, 1), domain(i)%fbcx(1:2, 1)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx(1:2, 2), domain(i)%fbcx(1:2, 2)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx(1:2, 3), domain(i)%fbcx(1:2, 3)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx(1:2, 4), domain(i)%fbcx(1:2, 4)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx(1:2, 5), domain(i)%fbcx(1:2, 5) ! dimensional
        end do
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy(1:2, 1), domain(1)%fbcy(1:2, 1)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy(1:2, 2), domain(1)%fbcy(1:2, 2)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy(1:2, 3), domain(1)%fbcy(1:2, 3)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy(1:2, 4), domain(1)%fbcy(1:2, 4)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy(1:2, 5), domain(1)%fbcy(1:2, 5) ! dimensional

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz(1:2, 1), domain(1)%fbcz(1:2, 1)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz(1:2, 2), domain(1)%fbcz(1:2, 2)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz(1:2, 3), domain(1)%fbcz(1:2, 3)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz(1:2, 4), domain(1)%fbcz(1:2, 4)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz(1:2, 5), domain(1)%fbcz(1:2, 5) ! dimensional

        do i = 1, nxdomain
          domain(i)%ibcy(:, :) = domain(1)%ibcy(:, :)
          domain(i)%fbcy(:, :) = domain(1)%fbcy(:, :)
          domain(i)%ibcz(:, :) = domain(1)%ibcz(:, :)
          domain(i)%fbcz(:, :) = domain(1)%fbcz(:, :)
        end do

        ! below for temporary fix
        do i = 1, nxdomain
          do j = 1, 2
            do k = 1, 5
              if(domain(i)%ibcx(j, k) > IBC_INTRPL) domain(i)%ibcx(j, k) = IBC_INTRPL
              if(domain(i)%ibcy(j, k) > IBC_INTRPL) domain(i)%ibcy(j, k) = IBC_INTRPL
              if(domain(i)%ibcz(j, k) > IBC_INTRPL) domain(i)%ibcz(j, k) = IBC_INTRPL
            end do
          end do
        end do

        if(nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt2i) '  u-x-bc-type  :', domain(i)%ibcx(1:2, 1)
            write (OUTPUT_UNIT, wrtfmt2r) '  u-x-bc-value :', domain(i)%fbcx(1:2, 1)
            write (OUTPUT_UNIT, wrtfmt2i) '  v-x-bc-type  :', domain(i)%ibcx(1:2, 2)
            write (OUTPUT_UNIT, wrtfmt2r) '  v-x-bc-value :', domain(i)%fbcx(1:2, 2)
            write (OUTPUT_UNIT, wrtfmt2i) '  w-x-bc-type  :', domain(i)%ibcx(1:2, 3)
            write (OUTPUT_UNIT, wrtfmt2r) '  w-x-bc-value :', domain(i)%fbcx(1:2, 3)
            write (OUTPUT_UNIT, wrtfmt2i) '  p-x-bc-type  :', domain(i)%ibcx(1:2, 4)
            write (OUTPUT_UNIT, wrtfmt2r) '  p-x-bc-value :', domain(i)%fbcx(1:2, 4)
            write (OUTPUT_UNIT, wrtfmt2i) '  T-x-bc-type  :', domain(i)%ibcx(1:2, 5)
            write (OUTPUT_UNIT, wrtfmt2r) '  T-x-bc-value :', domain(i)%fbcx(1:2, 5)
          end do
        end if
!-------------------------------------------------------------------------------
! [mesh] 
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[mesh]' ) then
        read(inputUnit, *, iostat = ioerr) varname, domain(1:nxdomain)%nc(1)
        read(inputUnit, *, iostat = ioerr) varname, itmp 
        domain(:)%nc(2) = itmp
        read(inputUnit, *, iostat = ioerr) varname, itmp 
        domain(:)%nc(3) = itmp
        read(inputUnit, *, iostat = ioerr) varname, itmp
        domain(:)%istret = itmp
        read(inputUnit, *, iostat = ioerr) varname, rtmp
        domain(:)%rstret = rtmp

        if(nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  mesh cell number - x     :', domain(i)%nc(1)
            write (OUTPUT_UNIT, wrtfmt1i) '  mesh cell number - y     :', domain(i)%nc(2)
            write (OUTPUT_UNIT, wrtfmt1i) '  mesh cell number - z     :', domain(i)%nc(3)
            write (OUTPUT_UNIT, wrtfmt1i) '  mesh y-stretching type   :', domain(i)%istret
            write (OUTPUT_UNIT, wrtfmt1r) '  mesh y-stretching factor :', domain(i)%rstret
          end do
        end if
!-------------------------------------------------------------------------------
! [timestepping]
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[timestepping]' ) then
        read(inputUnit, *, iostat = ioerr) varname, rtmp
        domain(:)%dt = rtmp
        read(inputUnit, *, iostat = ioerr) varname, itmp
        domain(:)%iTimeScheme = itmp

        if(nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1r) '  physical time step(dt) :', domain(i)%dt
            write (OUTPUT_UNIT, wrtfmt1i) '  time marching scheme   :', domain(i)%iTimeScheme
          end do
        end if
!-------------------------------------------------------------------------------
! [schemes]
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[schemes]' )  then
        read(inputUnit, *, iostat = ioerr) varname, itmp
        domain(:)%iAccuracy = itmp
        read(inputUnit, *, iostat = ioerr) varname, itmp
        domain(:)%iviscous = itmp

        if(nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  spatial accuracy scheme :', domain(i)%iAccuracy
            write (OUTPUT_UNIT, wrtfmt1i) '  viscous term treatment  :', domain(i)%iviscous
          end do
        end if
!-------------------------------------------------------------------------------
! [flow]
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[flow]' ) then
        read(inputUnit, *, iostat = ioerr) varname, is_any_energyeq
        if(is_any_energyeq) allocate( thermo(nxdomain) )
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%ren
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%idriven
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%drvfc      
        
        if(nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1r) '  flow Reynolds number      :', flow(i)%ren
            write (OUTPUT_UNIT, wrtfmt1i) '  flow driven force type    :', flow(i)%idriven
            write (OUTPUT_UNIT, wrtfmt1r) '  flow driven force         :', flow(i)%drvfc
          end do
        end if
!-------------------------------------------------------------------------------
! [thermo] 
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[thermo]' )  then 
        allocate ( itmpx(nxdomain) ); itmpx = 0
        allocate ( rtmpx(nxdomain) ); rtmpx = ZERO

        read(inputUnit, *, iostat = ioerr) varname, itmp
        if(is_any_energyeq) thermo(1 : nxdomain)%ifluid = itmp
        read(inputUnit, *, iostat = ioerr) varname, rtmp
        if(is_any_energyeq) thermo(1 : nxdomain)%lenRef = rtmp
        read(inputUnit, *, iostat = ioerr) varname, rtmp
        if(is_any_energyeq) thermo(1 : nxdomain)%T0Ref = rtmp
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%ithermo
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%icht
        read(inputUnit, *, iostat = ioerr) varname,   flow(1 : nxdomain)%igravity
        read(inputUnit, *, iostat = ioerr) varname, rtmpx(1: nxdomain)
        if(is_any_energyeq)  thermo(1 : nxdomain)%Tini0 = rtmpx(1: nxdomain)

        if(is_any_energyeq .and. nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  fluid medium               :', thermo(i)%ifluid
            write (OUTPUT_UNIT, wrtfmt1r) '  reference length (m)       :', thermo(i)%lenRef
            write (OUTPUT_UNIT, wrtfmt1r) '  reference temperature (K)  :', thermo(i)%T0Ref
            write (OUTPUT_UNIT, wrtfmt1i) '  is thermal field solved    ?', domain(i)%ithermo
            write (OUTPUT_UNIT, wrtfmt1i) '  is CHT solved              ?', domain(i)%icht
            write (OUTPUT_UNIT, wrtfmt1i) '  gravity direction          :', flow(i)%igravity
            write (OUTPUT_UNIT, wrtfmt1r) '  initial temperature (K)    :', thermo(i)%Tini0
          end do
        end if
!-------------------------------------------------------------------------------
! [initialization]
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[initialization]' ) then
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%irestart
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%nrsttckpt
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%nIterIniRen
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%renIni
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%initNoise
        do i = 1, nxdomain
          if(flow(i)%irestart /= 1) flow(i)%nrsttckpt = 0
        end do
        if( nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  flow restart                           :', flow(i)%irestart
            write (OUTPUT_UNIT, wrtfmt1i) '  restarting from iteration              :', flow(i)%nrsttckpt
            write (OUTPUT_UNIT, wrtfmt1i) '  Iteration for initial Reynolds No.     :', flow(i)%nIterIniRen
            write (OUTPUT_UNIT, wrtfmt1r) '  Initial Reynolds No.                   :', flow(i)%renIni
            write (OUTPUT_UNIT, wrtfmt1r) '  Initial velocity perturbation itensity :', flow(i)%initNoise
          end do
        end if
!-------------------------------------------------------------------------------
! [simcontrol]
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[simcontrol]' ) then
        read(inputUnit, *, iostat = ioerr) varname,   flow(1 : nxdomain)%nIterFlowStart
        read(inputUnit, *, iostat = ioerr) varname,   flow(1 : nxdomain)%nIterFlowEnd
        read(inputUnit, *, iostat = ioerr) varname,   itmpx(1:nxdomain)
        if(is_any_energyeq) thermo(1 : nxdomain)%nIterThermoStart = itmpx(1:nxdomain)
        read(inputUnit, *, iostat = ioerr) varname,   itmpx(1:nxdomain)
        if(is_any_energyeq) thermo(1 : nxdomain)%nIterThermoEnd = itmpx(1:nxdomain)

        if( nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  flow simulation starting from iteration    :', flow(i)%nIterFlowStart
            write (OUTPUT_UNIT, wrtfmt1i) '  flow simulation ending   at   iteration    :', flow(i)%nIterFlowEnd
            if(is_any_energyeq) then
            write (OUTPUT_UNIT, wrtfmt1i) '  thermal simulation starting from iteration :', thermo(i)%nIterThermoStart
            write (OUTPUT_UNIT, wrtfmt1i) '  thermal simulation ending   at   iteration :', thermo(i)%nIterThermoEnd
            end if
          end do
        end if
!-------------------------------------------------------------------------------
! [ioparams]
!-------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[ioparams]' ) then
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%nfreqckpt
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%nvisu
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%nIterStatsStart
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%nfreqStats

        if( nrank == 0) then
          do i = 1, nxdomain
            write (OUTPUT_UNIT, wrtfmt1i) 'For the domain-x  = ', i
            write (OUTPUT_UNIT, wrtfmt1i) '  raw  data written freqency  :', domain(i)%nfreqckpt
            write (OUTPUT_UNIT, wrtfmt1i) '  visu data written freqency  :', domain(i)%nvisu
            write (OUTPUT_UNIT, wrtfmt1i) '  statistics written from     :', domain(i)%nIterStatsStart
            write (OUTPUT_UNIT, wrtfmt1i) '  statistics written freqency :', domain(i)%nfreqStats
          end do
        end if
      else
        exit
      end if block_secname
    end do
!-------------------------------------------------------------------------------
! end of reading
!-------------------------------------------------------------------------------
    if(ioerr /= IOSTAT_END) &
    call Print_error_msg( 'Problem reading '//flname // &
    'in Subroutine: '// "Read_general_input")

    close(inputUnit)

    if(allocated(itmpx)) deallocate(itmpx)
    if(allocated(rtmpx)) deallocate(rtmpx)
!-------------------------------------------------------------------------------
! adjust input varnames
!-------------------------------------------------------------------------------
    do i = 1, nxdomain
!-------------------------------------------------------------------------------
!     coordinates
!-------------------------------------------------------------------------------
      if (domain(i)%icase == ICASE_CHANNEL) then
        domain(i)%icoordinate = ICARTESIAN
      else if (domain(i)%icase == ICASE_PIPE) then
        domain(i)%icoordinate = ICYLINDRICAL
        domain(i)%ibcy(:, 1) = IBC_INTERIOR
      else if (domain(i)%icase == ICASE_ANNUAL) then
        domain(i)%icoordinate = ICYLINDRICAL
      else if (domain(i)%icase == ICASE_TGV2D) then
        domain(i)%icoordinate = ICARTESIAN
      else if (domain(i)%icase == ICASE_TGV3D) then
        domain(i)%icoordinate = ICARTESIAN
      else if (domain(i)%icase == ICASE_BURGERS) then
        domain(i)%icoordinate = ICARTESIAN
      else 
        domain(i)%icoordinate = ICARTESIAN
      end if
!-------------------------------------------------------------------------------
!     stretching
!-------------------------------------------------------------------------------
      domain(i)%is_stretching(:) = .false.
      if(domain(i)%istret /= ISTRET_NO) domain(i)%is_stretching(2) = .true.
!-------------------------------------------------------------------------------
!     domain size
!-------------------------------------------------------------------------------
      if (domain(i)%icase == ICASE_CHANNEL) then
        if(domain(i)%istret /= ISTRET_2SIDES .and. nrank == 0) &
        call Print_warning_msg ("Grids are not two-side clustered.")
        domain(i)%lyb = - ONE
        domain(i)%lyt = ONE
      else if (domain(i)%icase == ICASE_PIPE) then
        if(domain(i)%istret /= ISTRET_TOP .and. nrank == 0) &
        call Print_warning_msg ("Grids are not near-wall clustered.")
        domain(i)%lyb = ZERO
        domain(i)%lyt = ONE
      else if (domain(i)%icase == ICASE_ANNUAL) then
        if(domain(i)%istret /= ISTRET_2SIDES .and. nrank == 0 ) &
        call Print_warning_msg ("Grids are not two-side clustered.")
        domain(i)%lyt = ONE
      else if (domain(i)%icase == ICASE_TGV2D) then
        if(domain(i)%istret /= ISTRET_NO .and. nrank == 0) &
        call Print_warning_msg ("Grids are clustered.")
        domain(i)%lxx = TWO * PI
        domain(i)%lzz = TWO * PI
        domain(i)%lyt =   PI
        domain(i)%lyb = - PI
      else if (domain(i)%icase == ICASE_TGV3D) then
        if(domain(i)%istret /= ISTRET_NO .and. nrank == 0) &
        call Print_warning_msg ("Grids are clustered.")
        domain(i)%lxx = TWO * PI
        domain(i)%lzz = TWO * PI
        domain(i)%lyt = TWO * PI
        domain(i)%lyb = ZERO
      else if (domain(i)%icase == ICASE_BURGERS) then
        if(domain(i)%istret /= ISTRET_NO .and. nrank == 0) &
        call Print_warning_msg ("Grids are clustered.")
        domain(i)%lxx = TWO
        domain(i)%lzz = TWO
        domain(i)%lyt = TWO
        domain(i)%lyb = ZERO
      else 
        ! do nothing...
      end if

      write (OUTPUT_UNIT, wrtfmt1r) '  scaled length in x-dir :', domain(i)%lxx
      write (OUTPUT_UNIT, wrtfmt1r) '  scaled length in y-dir :', domain(i)%lyt - domain(i)%lyb
      write (OUTPUT_UNIT, wrtfmt1r) '  scaled length in z-dir :', domain(i)%lzz
!-------------------------------------------------------------------------------
!     boundary
!-------------------------------------------------------------------------------
      if(is_any_energyeq) call Convert_thermo_BC_from_dim_to_undim(domain(i), thermo(i))

      domain(i)%is_periodic(:) = .false.
      do j = 1, 5
        if(domain(i)%ibcx(1, j) == IBC_PERIODIC .or. domain(i)%ibcx(2, j) == IBC_PERIODIC) then
          domain(i)%ibcx(1:2, j) = IBC_PERIODIC
          domain(i)%is_periodic(1) = .true.
        end if
        if(domain(i)%ibcy(1, j) == IBC_PERIODIC .or. domain(i)%ibcy(2, j) == IBC_PERIODIC) then
          domain(i)%ibcy(1:2, j) = IBC_PERIODIC
          domain(i)%is_periodic(2) = .true.
        end if
        if(domain(i)%ibcz(1, j) == IBC_PERIODIC .or. domain(i)%ibcz(2, j) == IBC_PERIODIC) then
          domain(i)%ibcz(1:2, j) = IBC_PERIODIC
          domain(i)%is_periodic(3) = .true.
        end if
      end do
!-------------------------------------------------------------------------------
! time and scheme 
!-------------------------------------------------------------------------------
    !option 1: to set up pressure treatment, for O(dt^2)
      domain(i)%sigma1p = ONE
      domain(i)%sigma2p = HALF
  
      !option 2: to set up pressure treatment, for O(dt)
      !sigma1p = ONE
      !sigma2p = ONE
  
      if(domain(i)%iTimeScheme == ITIME_RK3     .or. &
         domain(i)%iTimeScheme == ITIME_RK3_CN) then
        
        domain(i)%nsubitr = 3
        domain(i)%tGamma(0) = ONE
        domain(i)%tGamma(1) = EIGHT / FIFTEEN
        domain(i)%tGamma(2) = FIVE / TWELVE
        domain(i)%tGamma(3) = THREE / FOUR
  
        domain(i)%tZeta (0) = ZERO
        domain(i)%tZeta (1) = ZERO
        domain(i)%tZeta (2) = -SEVENTEEN / SIXTY
        domain(i)%tZeta (3) = -FIVE / TWELVE
  
      else if (domain(i)%iTimeScheme == ITIME_AB2) then
  
        domain(i)%nsubitr = 1
        domain(i)%tGamma(0) = ONE
        domain(i)%tGamma(1) = THREE / TWO
        domain(i)%tGamma(2) = ZERO
        domain(i)%tGamma(3) = ZERO
  
        domain(i)%tZeta (0) = ZERO
        domain(i)%tZeta (1) = -HALF
        domain(i)%tZeta (2) = ZERO
        domain(i)%tZeta (3) = ZERO
  
      else 
  
        domain(i)%nsubitr = 0
        domain(i)%tGamma(:) = ZERO
        domain(i)%tZeta (:) = ZERO
  
      end if 
      
      domain(i)%tAlpha(:) = domain(i)%tGamma(:) + domain(i)%tZeta(:)

    end do

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine Read_input_parameters

end module
