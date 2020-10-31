!##############################################################################
module input_mod
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
  public :: Initialize_input

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

  subroutine Initialize_input ()
    use iso_fortran_env, only : ERROR_UNIT, IOSTAT_END
    implicit none

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, input_unit

    character(len = 80) :: section_name
    character(len = 80) :: variable_name
    integer :: slen

    open ( newunit = input_unit, file = INPUT_FILE, status = 'old', action  = 'read', &
          iostat = ioerr, iomsg = iotxt)
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', INPUT_FILE, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 1
    end if
    
    do 
      read(input_unit, '(a)', iostat = ioerr) section_name
      slen = len_trim(section_name)
      if (ioerr /=0 ) exit
      if ( (section_name(1:1) == ';') .or. &
          (section_name(1:1) == '#') .or. &
          (section_name(1:1) == ' ') .or. &
          (slen == 0) ) then
        cycle
      end if

      block_section: if ( section_name(1:slen) == '[flowtype]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, icase
        read(input_unit, *, iostat = ioerr) variable_name, ithermo
        read(input_unit, *, iostat = ioerr) variable_name, icht

      else if ( section_name(1:slen) == '[decomposition]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, p_row
        read(input_unit, *, iostat = ioerr) variable_name, p_col

      else if ( section_name(1:slen) == '[geometry]' )  then 

        read(input_unit, *, iostat = ioerr) variable_name, lxx
        read(input_unit, *, iostat = ioerr) variable_name, lyt
        read(input_unit, *, iostat = ioerr) variable_name, lyb
        read(input_unit, *, iostat = ioerr) variable_name, lzz

      else if ( section_name(1:slen) == '[mesh]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, ncx
        read(input_unit, *, iostat = ioerr) variable_name, ncy
        read(input_unit, *, iostat = ioerr) variable_name, ncz
        read(input_unit, *, iostat = ioerr) variable_name, istret
        read(input_unit, *, iostat = ioerr) variable_name, rstret

      else if ( section_name(1:slen) == '[flowparams]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, ren

      else if ( section_name(1:slen) == '[timestepping]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, dt
        read(input_unit, *, iostat = ioerr) variable_name, iterFlowFirst
        read(input_unit, *, iostat = ioerr) variable_name, iterFlowLast

      else if ( section_name(1:slen) == '[boundary]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, ifbcx(1), ifbcx(2)
        read(input_unit, *, iostat = ioerr) variable_name, ifbcy(1), ifbcy(2)
        read(input_unit, *, iostat = ioerr) variable_name, ifbcz(1), ifbcz(2)

      else if ( section_name(1:slen) == '[ioparams]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, irestart
        read(input_unit, *, iostat = ioerr) variable_name, ncheckpoint
        read(input_unit, *, iostat = ioerr) variable_name, nvisu
        read(input_unit, *, iostat = ioerr) variable_name, iterStatsFirst
        read(input_unit, *, iostat = ioerr) variable_name, nstats

      else if ( section_name(1:slen) == '[schemes]' ) then
        
        read(input_unit, *, iostat = ioerr) variable_name, itimesteping
        read(input_unit, *, iostat = ioerr) variable_name, iviscous
        read(input_unit, *, iostat = ioerr) variable_name, ipressure

      else if ( section_name(1:slen) == '[initialization]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, renIni
        read(input_unit, *, iostat = ioerr) variable_name, iterRenIniEnd
        read(input_unit, *, iostat = ioerr) variable_name, initNoise

      else if ( section_name(1:slen) == '[periodicdriven]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, idriven

      else if ( section_name(1:slen) == '[thermohydraulics]' ) then

        read(input_unit, *, iostat = ioerr) variable_name, ifluid
        read(input_unit, *, iostat = ioerr) variable_name, igravity
        read(input_unit, *, iostat = ioerr) variable_name, lenRef
        read(input_unit, *, iostat = ioerr) variable_name, t0Ref
        read(input_unit, *, iostat = ioerr) variable_name, tiRef
        read(input_unit, *, iostat = ioerr) variable_name, itbcy(1), itbcy(2)
        read(input_unit, *, iostat = ioerr) variable_name, tbcy(1), tbcy(2)
        read(input_unit, *, iostat = ioerr) variable_name, iterThermoFirst
        read(input_unit, *, iostat = ioerr) variable_name, iterThermoLast

      else
        exit
      end if block_section
    end do

    if(ioerr /= IOSTAT_END) then
      write (ERROR_UNIT, *) 'Problem reading ', INPUT_FILE
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 3
    end if

    close(input_unit)

    ! to set up periodic boundary conditions
    call Set_periodic_bc ( ifbcx, is_x_periodic )
    call Set_periodic_bc ( ifbcy, is_y_periodic )
    call Set_periodic_bc ( ifbcz, is_z_periodic )

    ! to set up other variables derived from input variables
    npx = ncx + 1
    npy = ncy + 1
    npz = ncz + 1

  end subroutine Initialize_input
end module 





