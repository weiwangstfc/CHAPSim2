module property_table_mod
  use input_thermo_mod
  implicit none

  integer :: nlist
  real(WP), save, allocatable, dimension(:) :: listH
  real(WP), save, allocatable, dimension(:) :: listT
  real(WP), save, allocatable, dimension(:) :: listD
  real(WP), save, allocatable, dimension(:) :: listM
  real(WP), save, allocatable, dimension(:) :: listK
  real(WP), save, allocatable, dimension(:) :: listB
  real(WP), save, allocatable, dimension(:) :: listCp
  real(WP), save, allocatable, dimension(:) :: listDH

  real(WP) :: d0ref
  real(WP) :: h0ref
  real(WP) :: m0ref
  real(WP) :: k0ref
  real(WP) :: b0ref
  real(WP) :: Cp0ref

  real(WP) :: diref
  real(WP) :: hiref
  real(WP) :: miref
  real(WP) :: kiref
  real(WP) :: biref
  real(WP) :: Cpiref

  public :: Building_property_relations

contains

  subroutine Building_property_relations ( )
    use iso_fortran_env, only : ERROR_UNIT, IOSTAT_END
    use table_index_locating_mod

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    character(len = 80) :: str
    real(WP) :: rtmp
    real(WP) :: bufT, bufH, bufD, bufM, bufK, bufB, bufCp
    integer :: i, k
    integer :: i1, i2
    real(WP) :: weight1, weight2
 
    open ( newunit = inputUnit, file = inputProperty, status = 'old', action  = 'read', &
          iostat = ioerr, iomsg = iotxt)
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', INPUT_FILE, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 4
    end if

    read(inputUnit, *, iostat = ioerr) str
    read(inputUnit, *, iostat = ioerr) nlist

    allocate ( listH (nlist) )
    allocate ( listT (nlist) )
    allocate ( listD (nlist) )
    allocate ( listM (nlist) )
    allocate ( listK (nlist) )
    allocate ( listB (nlist) )
    allocate ( listCp(nlist) )
    allocate ( listDH(nlist) )

    block_tablereading: do i = 1, nlist
      read(inputUnit, *, iostat = ioerr) rtmp, listH(i), listT(i), listD(i), &
        listM(i), listK(i), listCp(i), listB(i)
    end do block_tablereading

    block_Tsorting: do i =  1, nlist
      k = minloc ( listT(i : nlist) , dim = 1) + i - 1
      bufT  = listT(i)
      bufD  = listD(i)
      bufH  = listH(i)
      bufM  = listM(i)
      bufK  = listK(i)
      bufB  = listB(i)
      bufCp = listCp(i)

      listT(i)  = listT(k)
      listD(i)  = listD(k)
      listH(i)  = listH(k)
      listM(i)  = listM(k)
      listK(i)  = listK(k)
      listB(i)  = listB(k)
      listCp(i) = listCp(k)

      listT(k)  = bufT
      listD(k)  = bufD
      listH(k)  = bufH
      listM(k)  = bufM
      listK(k)  = bufK
      listB(k)  = bufB
      listCp(k) = bufCp
    end do block_Tsorting

    d0ref  = Map_variables_from_list( t0Ref, listT(:), listD(:) )
    h0ref  = Map_variables_from_list( t0Ref, listT(:), listH(:) )
    m0ref  = Map_variables_from_list( t0Ref, listT(:), listM(:) )
    k0ref  = Map_variables_from_list( t0Ref, listT(:), listK(:) )
    b0ref  = Map_variables_from_list( t0Ref, listT(:), listB(:) )
    cp0ref = Map_variables_from_list( t0Ref, listT(:), listCp(:) )

    listH(:)  = (listH(:) - h0ref)/ t0ref / cp0ref
    listT(:)  = listT(:) / t0ref
    listD(:)  = listD(:) / d0ref
    listM(:)  = listM(:) / m0ref
    listK(:)  = listK(:) / k0ref
    listB(:)  = listB(:) / b0ref
    listCp(:) = listCp(:) / cp0ref
    listDH(:) = listD(:) * listH(:)

  end subroutine Building_property_relations
end module