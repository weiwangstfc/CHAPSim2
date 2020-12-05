!##############################################################################
module input_general_mod
  use precision_mod
  implicit none

  character(len = 9), parameter :: INPUT_FILE = 'input.ini'

  integer, parameter :: ICASE_CHANNEL = 1, &
                        ICASE_PIPE    = 2, &
                        ICASE_ANNUAL  = 3, &
                        ICASE_TGV     = 4

  integer, parameter :: ICARTESIAN   = 1, &
                        ICYLINDRICAL = 2
                        
  integer, parameter :: ISTRET_NO     = 0, &
                        ISTRET_SIDES  = 1, &
                        ISTRET_BOTTOM = 2, &
                        ISTRET_TOP    = 3

  integer, parameter :: ITIME_RK3 = 3, &
                        ITIME_AB1 = 1

  integer, parameter :: IBC_PERIODC     = 0, &
                        IBC_WALL_NOSLIP = 1, &
                        IBC_WALL_SLIP   = 2, &
                        IBC_INLET_MEAN  = 3, &
                        IBC_INLET_TG    = 4, &
                        IBC_INLET_MAP   = 5, &
                        IBC_INLET_DB    = 6, &
                        IBC_OUTLET_EXPO = 7, &
                        IBC_OUTLET_CONV = 8



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
  integer :: iterFlowFirst
  integer :: iterFlowLast
  ! boundary condition
  
  integer :: ifbcx(1:2)
  integer :: ifbcy(1:2)
  integer :: ifbcz(1:2)
  ! InOutParam
  integer :: irestart
  integer :: ncheckpoint
  integer :: nvisu
  integer :: iterStatsFirst
  integer :: nstats
  ! NumOption
  integer :: itimesteping
  integer :: iviscous
  integer :: ipressure
  ! initial fields
  real(WP) :: renIni
  integer :: iterRenIniEnd
  real(WP) :: initNoise
  ! PeriodicDrv
  integer :: idriven
  ! ThermoParam
  integer :: ifluid
  integer :: igravity
  real(WP) :: lenRef
  real(WP) :: t0Ref
  real(WP) :: tiRef
  integer :: itbcy(1:2)
  real(WP) :: tbcy(1:2)
  integer :: iterThermoFirst
  integer :: iterThermoLast

  ! derived parameters
  logical :: is_periodic(3)
  integer :: npx, npy, npz
  integer :: icoordinate

  integer :: ntInner
  real(WP) :: tGamma(0 : 3)
  real(WP) :: tZeta (0 : 3)
  real(WP) :: tAlpha(0 : 3)

  ! procedure
  private :: Set_periodic_bc
  private :: Set_timestepping_coefficients
  public :: Initialize_general_input

contains

  subroutine Set_timestepping_coefficients()
    use parameters_constant_mod

    if(itimesteping == ITIME_RK3) then
      
      ntInner = 3
      tGamma(0) = ONE
      tGamma(1) = EIGHT / FIFTEEN
      tGamma(2) = FIVE / TWELVE
      tGamma(3) = THREE / FOUR

      tZeta (0) = ZERO
      tZeta (1) = ZERO
      tZeta (2) = -SEVENTEEN / SIXTY
      tZeta (3) = -FIVE / TWELVE

    else if (itimesteping == ITIME_AB1) then

      ntInner = 1
      tGamma(0) = ONE
      tGamma(1) = ONEPFIVE
      tGamma(2) = ZERO
      tGamma(3) = ZERO

      tZeta (0) = ZERO
      tZeta (1) = -ZPFIVE
      tZeta (2) = ZERO
      tZeta (3) = ZERO

    else 

      ntInner = 0
      tGamma(:) = ZERO
      tZeta (:) = ZERO

    end if 
    
    tAlpha(:) = tGamma(:) + tZeta(:)

  end subroutine Set_timestepping_coefficients

  subroutine Set_periodic_bc( bc, flg )
    integer, intent(inout) :: bc(1:2)
    logical, intent(out) :: flg

    if ( (bc(1) == IBC_PERIODC) .or. (bc(2) == IBC_PERIODC) ) then
      bc(1) = IBC_PERIODC
      bc(2) = IBC_PERIODC
      flg = .true.
    else 
      flg = .false.
    end if

  end subroutine Set_periodic_bc

  subroutine Initialize_general_input ()
    use iso_fortran_env, only : ERROR_UNIT, IOSTAT_END
    implicit none

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit

    character(len = 80) :: section_name
    character(len = 80) :: variableName
    integer :: slen

    open ( newunit = inputUnit, file = INPUT_FILE, status = 'old', action  = 'read', &
          iostat = ioerr, iomsg = iotxt)
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

      block_section: if ( section_name(1:slen) == '[flowtype]' ) then

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
        read(inputUnit, *, iostat = ioerr) variableName, iterFlowFirst
        read(inputUnit, *, iostat = ioerr) variableName, iterFlowLast

      else if ( section_name(1:slen) == '[boundary]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ifbcx(1), ifbcx(2)
        read(inputUnit, *, iostat = ioerr) variableName, ifbcy(1), ifbcy(2)
        read(inputUnit, *, iostat = ioerr) variableName, ifbcz(1), ifbcz(2)

      else if ( section_name(1:slen) == '[ioparams]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, irestart
        read(inputUnit, *, iostat = ioerr) variableName, ncheckpoint
        read(inputUnit, *, iostat = ioerr) variableName, nvisu
        read(inputUnit, *, iostat = ioerr) variableName, iterStatsFirst
        read(inputUnit, *, iostat = ioerr) variableName, nstats

      else if ( section_name(1:slen) == '[schemes]' ) then
        
        read(inputUnit, *, iostat = ioerr) variableName, itimesteping
        read(inputUnit, *, iostat = ioerr) variableName, iviscous
        read(inputUnit, *, iostat = ioerr) variableName, ipressure

      else if ( section_name(1:slen) == '[initialization]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, renIni
        read(inputUnit, *, iostat = ioerr) variableName, iterRenIniEnd
        read(inputUnit, *, iostat = ioerr) variableName, initNoise

      else if ( section_name(1:slen) == '[periodicdriven]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, idriven

      else if ( section_name(1:slen) == '[thermohydraulics]' ) then

        read(inputUnit, *, iostat = ioerr) variableName, ifluid
        read(inputUnit, *, iostat = ioerr) variableName, igravity
        read(inputUnit, *, iostat = ioerr) variableName, lenRef
        read(inputUnit, *, iostat = ioerr) variableName, t0Ref
        read(inputUnit, *, iostat = ioerr) variableName, tiRef
        read(inputUnit, *, iostat = ioerr) variableName, itbcy(1), itbcy(2)
        read(inputUnit, *, iostat = ioerr) variableName, tbcy(1), tbcy(2)
        read(inputUnit, *, iostat = ioerr) variableName, iterThermoFirst
        read(inputUnit, *, iostat = ioerr) variableName, iterThermoLast

      else
        exit
      end if block_section
    end do

    if(ioerr /= IOSTAT_END) then
      write (ERROR_UNIT, *) 'Problem reading ', INPUT_FILE
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 3
    end if

    close(inputUnit)

    ! to set up periodic boundary conditions
    is_periodic(:) = .false.
    call Set_periodic_bc ( ifbcx, is_periodic(1) )
    call Set_periodic_bc ( ifbcy, is_periodic(2) )
    call Set_periodic_bc ( ifbcz, is_periodic(3) )

    ! to set up other variables derived from input variables
    npx = ncx + 1
    npy = ncy + 1
    npz = ncz + 1
    if ( is_periodic(1) ) npx = ncx
    if ( is_periodic(2) ) npy = ncy
    if ( is_periodic(3) ) npz = ncz

    ! to set up cooridnates
    if (icase == ICASE_CHANNEL) then
      icoordinate = ICARTESIAN
    else if (icase == ICASE_PIPE) then
      icoordinate = ICYLINDRICAL
    else if (icase == ICASE_ANNUAL) then
      icoordinate = ICYLINDRICAL
    else if (icase == ICASE_TGV) then
      icoordinate = ICARTESIAN
    else 
      icoordinate = ICARTESIAN
    end if

    call Set_timestepping_coefficients ( )

  end subroutine Initialize_general_input
end module 



