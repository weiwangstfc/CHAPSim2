!##############################################################################
module input_general_mod
  use precision_mod
  implicit none

  character(len = 9), parameter :: INPUT_FILE = 'input.ini'

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

  ! derive parameters
  logical :: is_x_periodic
  logical :: is_y_periodic
  logical :: is_z_periodic

  integer :: npx, npy, npz

  private :: Set_periodic_bc
  public :: Initialize_general_input

contains

  subroutine Set_periodic_bc( bc, flg )
    integer, intent(inout) :: bc(1:2)
    logical, intent(out) :: flg

    if ( (bc(1) == 0) .or. (bc(2) == 0) ) then
      bc(1) = 0
      bc(2) = 0
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
    call Set_periodic_bc ( ifbcx, is_x_periodic )
    call Set_periodic_bc ( ifbcy, is_y_periodic )
    call Set_periodic_bc ( ifbcz, is_z_periodic )

    ! to set up other variables derived from input variables
    npx = ncx + 1
    npy = ncy + 1
    npz = ncz + 1

  end subroutine Initialize_general_input
end module 





