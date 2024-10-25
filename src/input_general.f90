!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================================================================================================
!> \file input_general.f90
!> \brief Reading the input parameters from the given file.
!> \author Wei Wang wei.wang@stfc.ac.uk
!> \date 11-05-2022, checked.
!==========================================================================================================
module input_general_mod
  use print_msg_mod
  implicit none
  public  :: Read_input_parameters

contains
!==========================================================================================================
!> \brief Reading the input parameters from the given file.     
!! Scope:  mpi    called-freq    xdomain
!!         all    once           all
!----------------------------------------------------------------------------------------------------------
! Arguments
!----------------------------------------------------------------------------------------------------------
!  mode           name          role                                           
!----------------------------------------------------------------------------------------------------------
!> \param[in]     none          NA
!> \param[out]    none          NA
!==========================================================================================================
  subroutine Read_input_parameters
    use wtformat_mod
    use mpi_mod
    use parameters_constant_mod
    use vars_df_mod
    use thermo_info_mod
    use boundary_conditions_mod
    use code_performance_mod
    implicit none
    character(len = 18) :: flinput = 'input_chapsim.ini'
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
    integer :: i, j, m, n
    logical :: is_any_energyeq
    
    if(nrank == 0) then
      call Print_debug_start_msg("CHAPSim2.0 Starts ...")
      write (*, wrtfmt1i) '  The precision is   :', WP
    end if
    is_any_energyeq = .false.

    !----------------------------------------------------------------------------------------------------------
    ! open file
    !----------------------------------------------------------------------------------------------------------
    open ( newunit = inputUnit, &
           file    = flinput, &
           status  = 'old', &
           action  = 'read', &
           iostat  = ioerr, &
           iomsg   = iotxt )
    if(ioerr /= 0) then
      write (*, *) 'Problem openning : ', flinput, ' for reading.'
      write (*, *) 'Message: ', trim (iotxt)
      stop 1
    end if

    if(nrank == 0) &
    call Print_debug_start_msg("Reading General Parameters from "//flinput//" ...")
    !----------------------------------------------------------------------------------------------------------
    ! reading input
    !----------------------------------------------------------------------------------------------------------
    do 
      !----------------------------------------------------------------------------------------------------------
      ! reading headings/comments
      !----------------------------------------------------------------------------------------------------------
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
      !----------------------------------------------------------------------------------------------------------
      ! [decomposition]
      !----------------------------------------------------------------------------------------------------------
      if ( secname(1:slen) == '[decomposition]' ) then

        read(inputUnit, *, iostat = ioerr) varname, nxdomain
        read(inputUnit, *, iostat = ioerr) varname, p_row
        read(inputUnit, *, iostat = ioerr) varname, p_col
        allocate( domain (nxdomain) )
        allocate(   flow (nxdomain) )

        do i = 1, nxdomain
          domain(i)%idom = i
        end do

        if(nrank == 0) then
          write (*, wrtfmt1i) '  x-dir domain number             :', nxdomain
          write (*, wrtfmt1i) '  y-dir domain number (mpi Row)   :', p_row
          write (*, wrtfmt1i) '  z-dir domain number (mpi Column):', p_col
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [domain]
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[domain]' ) then

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%icase
        domain(:)%icase = domain(1)%icase

        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%lxx

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%lyt
        domain(:)%lyt = domain(1)%lyt

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%lyb
        domain(:)%lyb = domain(1)%lyb

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%lzz
        domain(:)%lzz = domain(1)%lzz

        !----------------------------------------------------------------------------------------------------------
        !     restore domain size to default if not set properly
        !----------------------------------------------------------------------------------------------------------
        do i = 1, nxdomain
          if (domain(i)%icase == ICASE_CHANNEL) then
            domain(i)%lyb = - ONE
            domain(i)%lyt = ONE
          else if (domain(i)%icase == ICASE_PIPE) then
            domain(i)%lyb = ZERO
            domain(i)%lyt = ONE
          else if (domain(i)%icase == ICASE_ANNUAL) then
            domain(i)%lyt = ONE
          else if (domain(i)%icase == ICASE_TGV2D .or. domain(i)%icase == ICASE_TGV3D) then
            domain(i)%lxx = TWOPI
            domain(i)%lzz = TWOPI
            domain(i)%lyt =   PI
            domain(i)%lyb = - PI
          else if (domain(i)%icase == ICASE_BURGERS) then
            domain(i)%lxx = TWO
            domain(i)%lzz = TWO
            domain(i)%lyt = TWO
            domain(i)%lyb = ZERO
          else if (domain(i)%icase == ICASE_ALGTEST) then
            domain(i)%lxx = TWOPI
            domain(i)%lzz = TWOPI
            domain(i)%lyt = TWOPI
            domain(i)%lyb = ZERO
          else 
            ! do nothing...
          end if

          !----------------------------------------------------------------------------------------------------------
          ! coordinates type
          !----------------------------------------------------------------------------------------------------------
          if (domain(i)%icase == ICASE_PIPE) then
            domain(i)%icoordinate = ICYLINDRICAL
          else if (domain(i)%icase == ICASE_ANNUAL) then
            domain(i)%icoordinate = ICYLINDRICAL
          else 
            domain(i)%icoordinate = ICARTESIAN
          end if

        end do

        
        if(nrank == 0) then
          write (*, wrtfmt1s) '  icase options:'
          write (*, wrtfmt1s) '          0 = OTHERS'
          write (*, wrtfmt1s) '          1 = CHANNEL'
          write (*, wrtfmt1s) '          2 = PIPE'
          write (*, wrtfmt1s) '          3 = ANNUAL'
          write (*, wrtfmt1s) '          4 = TGV3D'
          write (*, wrtfmt1s) '          5 = TGV2D'
          write (*, wrtfmt1s) '          6 = BURGERS'
          write (*, wrtfmt1s) '          7 = ALGTEST'
          write (*, wrtfmt1s) '  coordinates system option :'
          write (*, wrtfmt1s) '          1 = Cartesian'
          write (*, wrtfmt1s) '          2 = Cylindrical'

          do i = 1, nxdomain
            write (*, wrtfmt1i) '------For the domain-x------ ', i
            write (*, wrtfmt1i) '  current icase id :', domain(i)%icase
            write (*, wrtfmt1i) '  current coordinates system :', domain(i)%icoordinate
            write (*, wrtfmt1r) '  scaled length in x-direction :', domain(i)%lxx
            write (*, wrtfmt1r) '  scaled length in y-direction :', domain(i)%lyt - domain(i)%lyb
            if((domain(i)%lyt - domain(i)%lyb) < ZERO) call Print_error_msg("Y length is smaller than zero.")
            write (*, wrtfmt1r) '  scaled length in z-direction :', domain(i)%lzz
          end do
        end if
        
      !----------------------------------------------------------------------------------------------------------
      ! [boundary] 
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[boundary]' ) then

        do i = 1, nxdomain
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx_nominal(1:2, 1), domain(i)%fbcx_const(1:2, 1)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx_nominal(1:2, 2), domain(i)%fbcx_const(1:2, 2)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx_nominal(1:2, 3), domain(i)%fbcx_const(1:2, 3)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx_nominal(1:2, 4), domain(i)%fbcx_const(1:2, 4)
          read(inputUnit, *, iostat = ioerr) varname, domain(i)%ibcx_nominal(1:2, 5), domain(i)%fbcx_const(1:2, 5) ! dimensional
        end do

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy_nominal(1:2, 1), domain(1)%fbcy_const(1:2, 1)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy_nominal(1:2, 2), domain(1)%fbcy_const(1:2, 2)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy_nominal(1:2, 3), domain(1)%fbcy_const(1:2, 3)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy_nominal(1:2, 4), domain(1)%fbcy_const(1:2, 4)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcy_nominal(1:2, 5), domain(1)%fbcy_const(1:2, 5) ! dimensional

        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz_nominal(1:2, 1), domain(1)%fbcz_const(1:2, 1)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz_nominal(1:2, 2), domain(1)%fbcz_const(1:2, 2)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz_nominal(1:2, 3), domain(1)%fbcz_const(1:2, 3)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz_nominal(1:2, 4), domain(1)%fbcz_const(1:2, 4)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ibcz_nominal(1:2, 5), domain(1)%fbcz_const(1:2, 5) ! dimensional


        do i = 2, nxdomain
          domain(i)%ibcy_nominal(:, :) = domain(1)%ibcy_nominal(:, :)
          domain(i)%ibcz_nominal(:, :) = domain(1)%ibcz_nominal(:, :)
          domain(i)%fbcy_const(:, :) = domain(1)%fbcy_const(:, :)
          domain(i)%fbcz_const(:, :) = domain(1)%fbcz_const(:, :)
        end do

        do i = 1, nxdomain
          domain(i)%is_periodic(:) = .false.
          do m = 1, NBC
            if(domain(i)%ibcx_nominal(1, m) == IBC_PERIODIC .or. &
               domain(i)%ibcx_nominal(2, m) == IBC_PERIODIC) then
               domain(i)%ibcx_nominal(1:2, m) = IBC_PERIODIC
               domain(i)%is_periodic(1) = .true.
            end if
            if(domain(i)%ibcy_nominal(1, m) == IBC_PERIODIC .or. &
               domain(i)%ibcy_nominal(2, m) == IBC_PERIODIC) then
               domain(i)%ibcy_nominal(1:2, m) = IBC_PERIODIC
               domain(i)%is_periodic(2) = .true.
            end if
            if(domain(i)%ibcz_nominal(1, m) == IBC_PERIODIC .or. &
               domain(i)%ibcz_nominal(2, m) == IBC_PERIODIC) then
               domain(i)%ibcz_nominal(1:2, m) = IBC_PERIODIC
               domain(i)%is_periodic(3) = .true.
            end if
          end do

          if (domain(i)%icase == ICASE_PIPE) then
            domain(i)%ibcy_nominal(1, :) = IBC_INTERIOR
            domain(i)%is_periodic(2) = .false.
          end if
          !----------------------------------------------------------------------------------------------------------
          ! to exclude non-resonable input
          !----------------------------------------------------------------------------------------------------------
          domain(i)%is_conv_outlet = .false.
          do m = 1, NBC
            if(domain(i)%ibcx_nominal(2, m) == IBC_PROFILE1D) call Print_error_msg(" This BC IBC_PROFILE1D is not supported.")
            do n = 1, 2
              if(domain(i)%ibcx_nominal(n, m) >  IBC_OTHERS   ) call Print_error_msg(" This xBC is not suported.")
              if(domain(i)%ibcy_nominal(n, m) >  IBC_OTHERS   ) call Print_error_msg(" This yBC is not suported.")
              if(domain(i)%ibcz_nominal(n, m) >  IBC_OTHERS   ) call Print_error_msg(" This zBC is not suported.")
              if(domain(i)%ibcy_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This yBC IBC_PROFILE1D is not supported.")
              if(domain(i)%ibcz_nominal(n, m) == IBC_PROFILE1D) call Print_error_msg(" This zBC IBC_PROFILE1D is not supported.")
            end do
            if(domain(i)%ibcx_nominal(2, m) == IBC_CONVECTIVE) domain(i)%is_conv_outlet = .true.
          end do 
        end do
      !----------------------------------------------------------------------------------------------------------
      ! [mesh] 
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[mesh]' ) then
        read(inputUnit, *, iostat = ioerr) varname, domain(1:nxdomain)%nc(1)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%nc(2)
        domain(:)%nc(2) = domain(1)%nc(2)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%nc(3)
        domain(:)%nc(3) = domain(1)%nc(3)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%istret
        domain(:)%istret = domain(1)%istret
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%rstret
        domain(:)%rstret = domain(1)%rstret

        do i = 1, nxdomain
          !----------------------------------------------------------------------------------------------------------
          !     stretching
          !----------------------------------------------------------------------------------------------------------
          domain(i)%is_stretching(:) = .false.
          if(domain(i)%istret /= ISTRET_NO) domain(i)%is_stretching(2) = .true.
          
          if (domain(i)%icase == ICASE_CHANNEL .and. &
              domain(i)%istret /= ISTRET_2SIDES .and. &
              domain(i)%istret /= ISTRET_NO ) then

            if(nrank == 0) call Print_warning_msg ("Grids are neither uniform nor two-side clustered.")
          
          else if (domain(i)%icase == ICASE_PIPE .and. &
                   domain(i)%istret /= ISTRET_TOP) then

            if(nrank == 0) call Print_warning_msg ("Grids are not near-wall clustered.")

          else if (domain(i)%icase == ICASE_ANNUAL .and. &
                   domain(i)%istret /= ISTRET_2SIDES .and. &
                   domain(i)%istret /= ISTRET_NO) then

            if(nrank == 0) call Print_warning_msg ("Grids are neither uniform nor two-side clustered.")

          else if (domain(i)%icase == ICASE_TGV2D .or. &
                   domain(i)%icase == ICASE_TGV3D .or. &
                   domain(i)%icase == ICASE_ALGTEST) then

            if(domain(i)%istret /= ISTRET_NO .and. nrank == 0) &
            call Print_warning_msg ("Grids are clustered.")

          else 
            ! do nothing...
          end if
        end do

        if(nrank == 0) then
          write (*, wrtfmt1s) '  mesh streching option '
          write (*, wrtfmt1s) '          0 = ISTRET_NO'
          write (*, wrtfmt1s) '          1 = ISTRET_CENTRE'
          write (*, wrtfmt1s) '          2 = ISTRET_2SIDES'
          write (*, wrtfmt1s) '          3 = ISTRET_BOTTOM'
          write (*, wrtfmt1s) '          4 = ISTRET_TOP'
          do i = 1, nxdomain
            write (*, wrtfmt1i) '------For the domain-x------ ', i
            write (*, wrtfmt1i) '  mesh cell number - x :', domain(i)%nc(1)
            write (*, wrtfmt1i) '  mesh cell number - y :', domain(i)%nc(2)
            write (*, wrtfmt1i) '  mesh cell number - z :', domain(i)%nc(3)
            write (*, wrtfmt3l) '  is mesh stretching in x, y, z :', domain(i)%is_stretching(1:3)
            write (*, wrtfmt1i) '  mesh y-stretching type   :', domain(i)%istret
            write (*, wrtfmt1r) '  mesh y-stretching factor :', domain(i)%rstret
          end do
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [timestepping]
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[timestepping]' ) then
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%dt
        domain(:)%dt = domain(1)%dt
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%iTimeScheme
        domain(:)%iTimeScheme = domain(1)%iTimeScheme

        if(nrank == 0) then
          do i = 1, nxdomain
            write (*, wrtfmt1i) '------For the domain-x------ ', i
            write (*, wrtfmt1e) '  physical time step(dt, unit = second) :', domain(i)%dt
            write (*, wrtfmt1i) '  time marching scheme   :', domain(i)%iTimeScheme
          end do
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [schemes]
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[schemes]' )  then
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%iAccuracy
        domain(:)%iAccuracy = domain(1)%iAccuracy
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%iviscous
        domain(:)%iviscous = domain(1)%iviscous

        if(domain(1)%iAccuracy == IACCU_CD2 .or. &
           domain(1)%iAccuracy ==IACCU_CD4) then

          domain(:)%is_compact_scheme = .false.

        else if (domain(1)%iAccuracy == IACCU_CP4 .or. &
                 domain(1)%iAccuracy == IACCU_CP6) then

          domain(:)%is_compact_scheme = .true.

        else
          call Print_error_msg("Input error for numerical schemes.")
        end if
        

        if(nrank == 0) then
          write (*, wrtfmt1s) '  spatial accuracy scheme options : '
          write (*, wrtfmt1s) '    1 = IACCU_CD2'
          write (*, wrtfmt1s) '    2 = IACCU_CD4'
          write (*, wrtfmt1s) '    3 = IACCU_CP4'
          write (*, wrtfmt1s) '    4 = IACCU_CP6'
          do i = 1, nxdomain
            write (*, wrtfmt1i) '  ------For the domain-x------ ', i
            write (*, wrtfmt1i) '  current spatial accuracy scheme :', domain(i)%iAccuracy
            write (*, wrtfmt1i) '  viscous term treatment  :', domain(i)%iviscous
          end do
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [flow]
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[flow]' ) then

        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%inittype
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%iterfrom
        read(inputUnit, *, iostat = ioerr) varname, flow(1)%init_velo3d(1:3)
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%noiselevel
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%reninit
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%initReTo
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%ren
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%idriven
        read(inputUnit, *, iostat = ioerr) varname, flow(1 : nxdomain)%drvfc   

        do i = 1, nxdomain

          if(flow(i)%inittype /= INIT_RESTART) flow(i)%iterfrom = 0
          
          flow(i)%init_velo3d(1:3) = flow(1)%init_velo3d(1:3)

          if(domain(i)%icase /= ICASE_CHANNEL .and. &
             domain(i)%icase /= ICASE_ANNUAL  .and. &
             domain(i)%icase /= ICASE_PIPE ) then
            flow(i)%idriven = IDRVF_NO
          end if
          
          if(domain(i)%ibcx_nominal(1, 1) /= IBC_PERIODIC .or. &
             domain(i)%ibcx_nominal(2, 1) /= IBC_PERIODIC) then 
            flow(i)%idriven = IDRVF_NO
          end if

        end do

        if( nrank == 0) then
          write (*, wrtfmt1s) '  Initialisation options:'
          write (*, wrtfmt1s) '          INIT_RESTART    = 0'
          write (*, wrtfmt1s) '          INIT_INTERPL    = 1'
          write (*, wrtfmt1s) '          INIT_RANDOM     = 2'
          write (*, wrtfmt1s) '          INIT_INLET      = 3'
          write (*, wrtfmt1s) '          INIT_GVCONST    = 4'
          write (*, wrtfmt1s) '          INIT_POISEUILLE = 5'
          write (*, wrtfmt1s) '          INIT_FUNCTION   = 6'
          do i = 1, nxdomain
            write (*, wrtfmt1i) '------For the domain-x------ ', i
            write (*, wrtfmt1i) '  flow initial type                  :', flow(i)%inittype
            write (*, wrtfmt1i) '  iteration starting from            :', flow(i)%iterfrom
            if(flow(i)%inittype == INIT_GVCONST) then
            write (*, wrtfmt3r) '  initial velocity u, v, w           :', flow(i)%init_velo3d(1:3)
            end if
            write (*, wrtfmt1r) '  Initial velocity influction level  :', flow(i)%noiselevel
            write (*, wrtfmt1r) '  Initial Reynolds No.               :', flow(i)%reninit
            write (*, wrtfmt1i) '  Iteration for initial Reynolds No. :', flow(i)%initReTo
            write (*, wrtfmt1r) '  flow Reynolds number               :', flow(i)%ren
            write (*, wrtfmt1i) '  flow driven force type             :', flow(i)%idriven
            write (*, wrtfmt1r) '  flow driven force(cf)              :', flow(i)%drvfc        
          end do
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [thermo] 
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[thermo]' )  then 
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%is_thermo
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%icht
        read(inputUnit, *, iostat = ioerr) varname,   flow(1 : nxdomain)%igravity

        if(ANY(domain(:)%is_thermo)) is_any_energyeq = .true.
        if(is_any_energyeq) allocate( thermo(nxdomain) )
        allocate ( itmpx(nxdomain) ); itmpx = 0
        allocate ( rtmpx(nxdomain) ); rtmpx = ZERO
        

        read(inputUnit, *, iostat = ioerr) varname, itmp
        if(is_any_energyeq) thermo(1 : nxdomain)%ifluid = itmp
        read(inputUnit, *, iostat = ioerr) varname, rtmp
        if(is_any_energyeq) thermo(1 : nxdomain)%ref_l0 = rtmp
        read(inputUnit, *, iostat = ioerr) varname, rtmp
        if(is_any_energyeq) thermo(1 : nxdomain)%ref_T0 = rtmp
        
        read(inputUnit, *, iostat = ioerr) varname, itmp
        if(is_any_energyeq) thermo(1 : nxdomain)%inittype  = itmp
        read(inputUnit, *, iostat = ioerr) varname, itmp
        if(is_any_energyeq) thermo(1 : nxdomain)%iterfrom = itmp
        read(inputUnit, *, iostat = ioerr) varname, rtmpx(1: nxdomain)
        if(is_any_energyeq) thermo(1 : nxdomain)%init_T0 = rtmpx(1: nxdomain)
        
        if(is_any_energyeq .and. nrank == 0) then
          do i = 1, nxdomain
            write (*, wrtfmt1i) '------For the domain-x------ ', i
            write (*, wrtfmt1l) '  is thermal field solved   ?', domain(i)%is_thermo
            write (*, wrtfmt1l) '  is CHT solved             ?', domain(i)%icht
            write (*, wrtfmt1i) '  gravity direction         :', flow(i)%igravity
            write (*, wrtfmt1i) '  fluid medium              :', thermo(i)%ifluid
            write (*, wrtfmt1r) '  reference length (m)      :', thermo(i)%ref_l0
            write (*, wrtfmt1r) '  reference temperature (K) :', thermo(i)%ref_T0
            write (*, wrtfmt1i) '  thermo field initial type :', thermo(i)%inittype
            write (*, wrtfmt1i) '  iteration starting from   :', thermo(i)%iterfrom
            write (*, wrtfmt1r) '  initial temperature (K)   :', thermo(i)%init_T0
          end do
        else if(nrank == 0) then
          write(*, *) ' Note: Thermal field is not considered. '
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [simcontrol]
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[simcontrol]' ) then
        read(inputUnit, *, iostat = ioerr) varname,   flow(1 : nxdomain)%nIterFlowStart
        read(inputUnit, *, iostat = ioerr) varname,   flow(1 : nxdomain)%nIterFlowEnd
        read(inputUnit, *, iostat = ioerr) varname,   itmpx(1:nxdomain)
        if(is_any_energyeq) thermo(1 : nxdomain)%nIterThermoStart = itmpx(1:nxdomain)
        read(inputUnit, *, iostat = ioerr) varname,   itmpx(1:nxdomain)
        if(is_any_energyeq) thermo(1 : nxdomain)%nIterThermoEnd = itmpx(1:nxdomain)

        if( nrank == 0) then
          do i = 1, nxdomain
            write (*, wrtfmt1i) '------For the domain-x------ ', i
            write (*, wrtfmt1i) '  flow simulation starting from iteration    :', flow(i)%nIterFlowStart
            write (*, wrtfmt1i) '  flow simulation ending   at   iteration    :', flow(i)%nIterFlowEnd
            if(is_any_energyeq) then
            write (*, wrtfmt1i) '  thermal simulation starting from iteration :', thermo(i)%nIterThermoStart
            write (*, wrtfmt1i) '  thermal simulation ending   at   iteration :', thermo(i)%nIterThermoEnd
            end if
          end do
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [ioparams]
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[ioparams]' ) then
        read(inputUnit, *, iostat = ioerr) varname, cpu_nfre
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%ckpt_nfre
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%visu_idim
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%visu_nfre
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%visu_nskip(1:3)
        read(inputUnit, *, iostat = ioerr) varname, domain(1 : nxdomain)%stat_istart
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%stat_nskip(1:3)
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%is_record_xoutlet, domain(1)%is_read_xinlet
        read(inputUnit, *, iostat = ioerr) varname, domain(1)%ndbfre, domain(1)%ndbend
        
        do i = 1, nxdomain
            domain(i)%visu_nskip(1:3) = domain(1)%visu_nskip(1:3)
            domain(i)%stat_nskip(1:3) = domain(1)%stat_nskip(1:3) 
           if(domain(i)%is_stretching(2)) domain(i)%visu_nskip(2) = 1
           if(domain(i)%is_stretching(2)) domain(i)%stat_nskip(2) = 1
        end do

        if( nrank == 0) then
          do i = 1, nxdomain
            write (*, wrtfmt1i) '------For the domain-x------ ', i
            write (*, wrtfmt1i) '  data check freqency  :', domain(i)%ckpt_nfre
            write (*, wrtfmt1i) '  visu data dimensions        :', domain(i)%visu_idim
            write (*, wrtfmt1i) '  visu data written freqency  :', domain(i)%visu_nfre
            write (*, wrtfmt3i) '  visu data skips in x, y, z  :', domain(i)%visu_nskip(1:3)
            write (*, wrtfmt1i) '  statistics written from     :', domain(i)%stat_istart
            write (*, wrtfmt3i) '  statistics skips in x, y, z :', domain(i)%stat_nskip(1:3)
          end do
        end if
      !----------------------------------------------------------------------------------------------------------
      ! [probe]
      !----------------------------------------------------------------------------------------------------------
      else if ( secname(1:slen) == '[probe]' ) then
        do i = 1, nxdomain
          read(inputUnit, *, iostat = ioerr) varname, itmp
          domain(i)%proben = itmp
          if(domain(i)%proben > 0) then
            allocate( domain(i)%probexyz(3, itmp))
            if( nrank == 0) write (*, wrtfmt1i) '------For the domain-x------ ', i
            do j = 1, domain(i)%proben
              read(inputUnit, *, iostat = ioerr) domain(i)%probexyz(1:3, j) 
              
              if( nrank == 0) write (*, wrtfmt3r) '  probed points x, y, z: ', domain(i)%probexyz(1:3, j) 
            end do 
          end if
        end do
      else
        exit
      end if
    end do
    !----------------------------------------------------------------------------------------------------------
    ! end of reading, clearing dummies
    !----------------------------------------------------------------------------------------------------------
    if(.not.IS_IOSTAT_END(ioerr)) &
    call Print_error_msg( 'Problem reading '//flinput // &
    'in Subroutine: '// "Read_general_input")

    close(inputUnit)

    if(allocated(itmpx)) deallocate(itmpx)
    if(allocated(rtmpx)) deallocate(rtmpx)
    !----------------------------------------------------------------------------------------------------------
    ! convert the input dimensional temperature/heat flux into undimensional
    !----------------------------------------------------------------------------------------------------------
    do i = 1, nxdomain
      call config_calc_basic_ibc(domain(i))
      call config_calc_eqs_ibc(domain(i))
    end do 
    !----------------------------------------------------------------------------------------------------------
    ! set up constant for time step marching 
    !----------------------------------------------------------------------------------------------------------
    do i = 1, nxdomain
      !option 1: to set up pressure treatment, for O(dt^2)
      !domain(i)%sigma1p = ONE
      !domain(i)%sigma2p = HALF
  
      !option 2: to set up pressure treatment, for O(dt)
      domain(i)%sigma1p = ONE
      domain(i)%sigma2p = ONE
  
      if(domain(i)%iTimeScheme == ITIME_RK3     .or. &
         domain(i)%iTimeScheme == ITIME_RK3_CN) then
        
        domain(i)%nsubitr = 3
        domain(i)%tGamma(0) = ONE
        domain(i)%tGamma(1) = EIGHT / FIFTEEN
        domain(i)%tGamma(2) = FIVE / TWELVE
        domain(i)%tGamma(3) = THREE * QUARTER
  
        domain(i)%tZeta (0) = ZERO
        domain(i)%tZeta (1) = ZERO
        domain(i)%tZeta (2) = - SEVENTEEN / SIXTY
        domain(i)%tZeta (3) = - FIVE / TWELVE
  
      else if (domain(i)%iTimeScheme == ITIME_AB2) then !Adams-Bashforth
  
        domain(i)%nsubitr = 1
        domain(i)%tGamma(0) = ONE
        domain(i)%tGamma(1) = ONEPFIVE
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
      
      domain(i)%tAlpha(0:3) = domain(i)%tGamma(0:3) + domain(i)%tZeta(0:3)

    end do

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine Read_input_parameters

end module
