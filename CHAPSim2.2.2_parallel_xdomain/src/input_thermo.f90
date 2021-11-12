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
!> \file input_thermo.f90
!>
!> \brief Reading the input parameters from the given file and building up the
!> relationships between properties.
!>
!===============================================================================
module input_thermo_mod
  use parameters_constant_mod
  implicit none

  type thermoProperty_t
    real(WP) :: t  ! temperature
    real(WP) :: d  ! density
    real(WP) :: m  ! dynviscosity
    real(WP) :: k  ! thermconductivity
    real(WP) :: h  ! enthalpy
    real(WP) :: dh ! mass enthalpy
    real(WP) :: cp ! specific heat capacity 
    real(WP) :: b  ! thermal expansion
  contains
    private
    procedure, public :: Get_initialized_thermal_properties
    procedure, public :: is_T_in_scope
    procedure, public :: Refresh_thermal_properties_from_T_dimensional
    procedure, public :: Refresh_thermal_properties_from_T_undim
    procedure, public :: Refresh_thermal_properties_from_H
    procedure, public :: Refresh_thermal_properties_from_DH
    procedure :: Print_debug
    generic :: Print => Print_debug
    generic :: write(formatted) => Print_debug
  end type thermoProperty_t

  type(thermoProperty_t), save, allocatable, dimension(:) :: listTP
  type(thermoProperty_t) :: tpRef0 ! dim, reference state

  private
  character(len = 64) :: inputProperty
  integer :: ipropertyState
  real(WP) :: u0dim
  real(WP) :: TM0
  real(WP) :: TB0
  real(WP) :: HM0
  real(WP) :: CoD(0:1)
  real(WP) :: CoK(0:2)
  real(WP) :: CoB
  real(WP) :: CoCp(-2:2)
  real(WP) :: CoH(-1:3)
  real(WP) :: CoM(-1:1)

  integer :: nlist

  private :: Buildup_property_relations_from_table
  private :: Buildup_property_relations_from_function
  private :: Check_monotonicity_DH_of_HT_list
  public  :: Buildup_thermo_mapping_relations
  private :: Initialize_thermo_parameters
  private :: Sort_listTP_Tsmall2big
  private :: Write_thermo_property

contains
!===============================================================================
!===============================================================================
!> \brief Defination of a procedure in the type thermoProperty_t.
!>  to initialize the default thermal properties.     
!>
!> This subroutine is called as required to initialize the default
!> thermal properties. It is used only when the \ref ithermo is 1.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  this          a cell element with udf property
!_______________________________________________________________________________
  subroutine Get_initialized_thermal_properties ( this )
    implicit none

    class(thermoProperty_t), intent(inout) :: this
    
    this%t  = ONE
    this%d  = ONE
    this%m  = ONE
    this%k  = ONE
    this%cp = ONE
    this%b  = ONE
    this%h  = ZERO
    this%dh = ZERO

    return
  end subroutine Get_initialized_thermal_properties
!===============================================================================
!===============================================================================
!> \brief Defination of a procedure in the type thermoProperty_t.
!>  to check the temperature limitations.     
!>
!> This subroutine is called as required to check the temperature
!> of given element is within the given limits as a single phase
!> flow.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  this          a cell element with udf property
!_______________________________________________________________________________
  subroutine is_T_in_scope ( this )
    implicit none

    class( thermoProperty_t ), intent( inout ) :: this
    
    if(ipropertyState == IPROPERTY_TABLE) then
      if ( ( this%t < listTP(1)%t     )  .OR. &
           ( this%t > listTP(nlist)%t ) ) then
        print*, this%t, listTP(1)%t, listTP(nlist)%t
        stop 'temperature exceeds specified range.'
      end if
    end if

    if(ipropertyState == IPROPERTY_FUNCS) then
      if ( ( this%t < ( TM0 / tpRef0%t ) ) .OR. &
           ( this%t > ( TB0 / tpRef0%t ) ) ) then 
        print*, this%t, TM0 / tpRef0%t, TB0 / tpRef0%t
        stop 'temperature exceeds specified range.'
      end if
    end if

    return
  end subroutine is_T_in_scope
!===============================================================================
!===============================================================================
!> \brief Defination of a procedure in the type thermoProperty_t.
!>  to update the thermal properties based on the known temperature.     
!>
!> This subroutine is called as required to update all thermal properties from
!> the known temperature (dimensional or dimensionless).
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  this          a cell element with udf property
!> \param[in]     dim           an optional indicator of dimensional
!>                              /dimensionless T. Exsiting of \ref dim indicates
!>                              the known/given T is dimensional. Otherwise, not.
!_______________________________________________________________________________
  subroutine Refresh_thermal_properties_from_T_dimensional ( this )
    use parameters_constant_mod, only : MINP, ONE, ZERO
    use input_general_mod, only : ifluid
    implicit none
    class(thermoProperty_t), intent(inout) :: this
    
    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2
    real(WP) :: t1, dummy

    if(ipropertyState == IPROPERTY_TABLE) then 
      i1 = 1
      i2 = nlist
      do while ( (i2 - i1) > 1)
        im = i1 + (i2 - i1) / 2
        d1 = listTP(i1)%t - this%t
        dm = listTP(im)%t - this%t
        if ( (d1 * dm) > MINP ) then
          i1 = im
        else
          i2 = im
        end if
      end do

      w1 = (listTP(i2)%t - this%t) / (listTP(i2)%t - listTP(i1)%t) 
      w2 = ONE - w1

      this%d = w1 * listTP(i1)%d + w2 * listTP(i2)%d
      this%m = w1 * listTP(i1)%m + w2 * listTP(i2)%m
      this%k = w1 * listTP(i1)%k + w2 * listTP(i2)%k
      this%h = w1 * listTP(i1)%h + w2 * listTP(i2)%h
      this%b = w1 * listTP(i1)%b + w2 * listTP(i2)%b
      this%cp = w1 * listTP(i1)%cp + w2 * listTP(i2)%cp
      this%dh = this%d * this%h

    else if(ipropertyState == IPROPERTY_FUNCS) then 
      
      t1 = this%t
  
      ! D = density = f(T)
      this%d = CoD(0) + CoD(1) * t1
  
      ! K = thermal conductivity = f(T)
      this%k = CoK(0) + CoK(1) * t1 + CoK(2) * t1**2
  
      ! Cp = f(T)
      this%cp = CoCp(-2) * t1**(-2) + CoCp(-1) * t1**(-1) + CoCp(0) + CoCp(1) * t1 + CoCp(2) * t1**2
  
      ! H = entropy = f(T)
      this%h = Hm0 + &
        CoH(-1) * (ONE / t1 - ONE / Tm0) + &
        CoH(0) + &
        CoH(1) * (t1 - Tm0) + &
        CoH(2) * (t1**2 - Tm0**2) + &
        CoH(3) * (t1**3 - Tm0**3)
  
      ! B = f(T)
      this%b = ONE / (CoB - t1)
  
      ! dynamic viscosity = f(T)
      select case (ifluid)
        ! unit: T(Kelvin), M(Pa S)
        case (ILIQUID_SODIUM)
          dummy = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
        case (ILIQUID_LEAD)
          dummy = CoM_Pb(0) * EXP (CoM_Pb(-1) / t1)
        case (ILIQUID_BISMUTH)
          dummy = CoM_Bi(0) * EXP (CoM_Bi(-1) / t1)
        case (ILIQUID_LBE)
          dummy = CoM_LBE(0) * EXP (CoM_LBE(-1) / t1)
        case default
          dummy = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
      end select
      this%m = dummy
      this%dh = this%d * this%h
    else
      this%t  = ONE
      this%d  = ONE
      this%m  = ONE
      this%k  = ONE
      this%cp = ONE
      this%b  = ONE
      this%h  = ZERO
      this%dh = ZERO
    end if
    return
  end subroutine Refresh_thermal_properties_from_T_dimensional
!===============================================================================
!===============================================================================
  subroutine Refresh_thermal_properties_from_T_undim ( this)
    use parameters_constant_mod, only : MINP, ONE, ZERO
    use input_general_mod, only : ifluid
    implicit none
    class(thermoProperty_t), intent(inout) :: this
    
    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2
    real(WP) :: t1, dummy

    if(ipropertyState == IPROPERTY_TABLE) then 
      i1 = 1
      i2 = nlist
      do while ( (i2 - i1) > 1)
        im = i1 + (i2 - i1) / 2
        d1 = listTP(i1)%t - this%t
        dm = listTP(im)%t - this%t
        if ( (d1 * dm) > MINP ) then
          i1 = im
        else
          i2 = im
        end if
      end do

      w1 = (listTP(i2)%t - this%t) / (listTP(i2)%t - listTP(i1)%t) 
      w2 = ONE - w1

      this%d = w1 * listTP(i1)%d + w2 * listTP(i2)%d
      this%m = w1 * listTP(i1)%m + w2 * listTP(i2)%m
      this%k = w1 * listTP(i1)%k + w2 * listTP(i2)%k
      this%h = w1 * listTP(i1)%h + w2 * listTP(i2)%h
      this%b = w1 * listTP(i1)%b + w2 * listTP(i2)%b
      this%cp = w1 * listTP(i1)%cp + w2 * listTP(i2)%cp
      this%dh = this%d * this%h

    else if(ipropertyState == IPROPERTY_FUNCS) then 
      
      ! convert undim to dim 
      t1 = this%t * tpRef0%t

      ! D = density = f(T)
      dummy = CoD(0) + CoD(1) * t1
      this%d = dummy / tpRef0%d
  
      ! K = thermal conductivity = f(T)
      dummy = CoK(0) + CoK(1) * t1 + CoK(2) * t1**2
      this%k = dummy / tpRef0%k
  
      ! Cp = f(T)
      dummy = CoCp(-2) * t1**(-2) + CoCp(-1) * t1**(-1) + CoCp(0) + CoCp(1) * t1 + CoCp(2) * t1**2
      this%cp = dummy / tpRef0%cp
  
      ! H = entropy = f(T)
      dummy = Hm0 + &
        CoH(-1) * (ONE / t1 - ONE / Tm0) + &
        CoH(0) + &
        CoH(1) * (t1 - Tm0) + &
        CoH(2) * (t1**2 - Tm0**2) + &
        CoH(3) * (t1**3 - Tm0**3)
      this%h = (dummy - tpRef0%h) / (tpRef0%cp * tpRef0%t)
  
      ! B = f(T)
      dummy = ONE / (CoB - t1)
      this%b = dummy / tpRef0%b
  
      ! dynamic viscosity = f(T)
      select case (ifluid)
        ! unit: T(Kelvin), M(Pa S)
        case (ILIQUID_SODIUM)
          dummy = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
        case (ILIQUID_LEAD)
          dummy = CoM_Pb(0) * EXP (CoM_Pb(-1) / t1)
        case (ILIQUID_BISMUTH)
          dummy = CoM_Bi(0) * EXP (CoM_Bi(-1) / t1)
        case (ILIQUID_LBE)
          dummy = CoM_LBE(0) * EXP (CoM_LBE(-1) / t1)
        case default
          dummy = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
      end select
      this%m = dummy / tpRef0%m
      this%dh = this%d * this%h
    else
      this%t  = ONE
      this%d  = ONE
      this%m  = ONE
      this%k  = ONE
      this%cp = ONE
      this%b  = ONE
      this%h  = ZERO
      this%dh = ZERO
    end if
    return
  end subroutine Refresh_thermal_properties_from_T_undim
!===============================================================================
!===============================================================================
!> \brief Defination of a procedure in the type thermoProperty_t.
!>  to update the thermal properties based on the known enthalpy.     
!>
!> This subroutine is called as required to update all thermal properties from
!> the known enthalpy (dimensionless only).
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  this          a cell element with udf property
!_______________________________________________________________________________
  subroutine Refresh_thermal_properties_from_H(this)
    use parameters_constant_mod, only : MINP, ONE
    class(thermoProperty_t), intent(inout) :: this

    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2

    i1 = 1
    i2 = nlist

    do while ( (i2 - i1) > 1)
      im = i1 + (i2 - i1) / 2
      d1 = listTP(i1)%h - this%h
      dm = listTP(im)%h - this%h
      if ( (d1 * dm) > MINP ) then
        i1 = im
      else
        i2 = im
      end if
    end do

    w1 = (listTP(i2)%h - this%h) / (listTP(i2)%h - listTP(i1)%h) 
    w2 = ONE - w1

    this%d = w1 * listTP(i1)%d + w2 * listTP(i2)%d
    this%m = w1 * listTP(i1)%m + w2 * listTP(i2)%m
    this%k = w1 * listTP(i1)%k + w2 * listTP(i2)%k
    this%t = w1 * listTP(i1)%t + w2 * listTP(i2)%t
    this%b = w1 * listTP(i1)%b + w2 * listTP(i2)%b
    this%cp = w1 * listTP(i1)%cp + w2 * listTP(i2)%cp
    this%dh = this%d * this%h
    return
  end subroutine Refresh_thermal_properties_from_H
!===============================================================================
!===============================================================================
!> \brief Defination of a procedure in the type thermoProperty_t.
!>  to update the thermal properties based on the known enthalpy per unit mass.     
!>
!> This subroutine is called as required to update all thermal properties from
!> the known enthalpy per unit mass (dimensionless only).
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  this          a cell element with udf property
!_______________________________________________________________________________
  subroutine Refresh_thermal_properties_from_DH(this)
    use parameters_constant_mod, only : MINP, ONE
    class(thermoProperty_t), intent(inout) :: this

    integer :: i1, i2, im
    real(WP) :: d1, dm
    real(WP) :: w1, w2
    logical :: is_dim

    i1 = 1
    i2 = nlist

    do while ( (i2 - i1) > 1)
      im = i1 + (i2 - i1) / 2
      d1 = listTP(i1)%dh - this%dh
      dm = listTP(im)%dh - this%dh
      if ( (d1 * dm) > MINP ) then
        i1 = im
      else
        i2 = im
      end if
    end do

    w1 = (listTP(i2)%dh - this%dh) / (listTP(i2)%dh - listTP(i1)%dh) 
    w2 = ONE - w1
    
    if(ipropertyState == IPROPERTY_TABLE) then 
      this%h = w1 * listTP(i1)%h + w2 * listTP(i2)%h
      call this%Refresh_thermal_properties_from_H
    else if (ipropertyState == IPROPERTY_FUNCS) then 
      this%t = w1 * listTP(i1)%t + w2 * listTP(i2)%t
      call this%Refresh_thermal_properties_from_T_undim
    else  
      STOP 'Error. No such option of ipropertyState.'
    end if
    return
  end subroutine Refresh_thermal_properties_from_DH
!===============================================================================
!===============================================================================
!> \brief Defination of a procedure in the type thermoProperty_t.
!>  to print out the thermal properties at the given element.     
!>
!> This subroutine is called as required to print out thermal properties 
!> for degbugging.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     this          a cell element with udf property
!> \param[in]     unit          
!> \param[in]     iotype        
!> \param[in]     v_list        
!> \param[out]    iostat        
!> \param[inout]  iomsg         
!_______________________________________________________________________________
  subroutine Print_debug(this, unit, iotype, v_list, iostat, iomsg)
    use iso_fortran_env, only : error_unit
    class(thermoProperty_t), intent(in) :: this
    integer, intent(in)                 :: unit
    character(len = *), intent(in)      :: iotype
    integer, intent(in)                 :: v_list(:)
    integer, intent(out)                :: iostat
    character(len = *), intent(inout)   :: iomsg

    integer                             :: i_pass

    iostat = 0
    iomsg = ""
    
    this_block: do i_pass = 1, 1
      !write(unit, *, iostat = iostat, iomsg = iomsg) 'thermalProperty'
      !if(iostat /= 0) exit this_block
      if(iotype(1:2) == 'DT' .and. len(iotype) > 2) &
        write(unit, *, iostat = iostat, iomsg = iomsg) iotype(3:)
      if(iostat /= 0) exit this_block

      write(unit, *, iostat = iostat, iomsg = iomsg) &
      this%h, this%t, this%d, this%m, this%k, this%cp, this%b, this%dh
      
      if(iostat /= 0) exit this_block
    end do this_block

    if(iostat /= 0) then
      write (error_unit, "(A)") "print error : " // trim(iomsg)
      write (error_unit, "(A, I0)") "  iostat : ", iostat
    end if
    return
  end subroutine Print_debug
!===============================================================================
!> \brief Sort out the user given thermal property table based on the temperature
!>  from small to big.
!>
!> This subroutine is called locally once reading in the given thermal table.
!>
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]  list         the thermal table element array
!===============================================================================
  subroutine Sort_listTP_Tsmall2big(list)
    type(thermoProperty_t),intent(inout) :: list(:)
    integer :: i, n, k
    real(WP) :: buf

    n = size( list )

    do i = 1, n
      k = minloc( list(i:n)%t, dim = 1) + i - 1

      buf = list(i)%t
      list(i)%t = list(k)%t
      list(k)%t = buf

      buf = list(i)%d
      list(i)%d = list(k)%d
      list(k)%d = buf

      buf = list(i)%m
      list(i)%m = list(k)%m
      list(k)%m = buf

      buf = list(i)%k
      list(i)%k = list(k)%k
      list(k)%k = buf

      buf = list(i)%b
      list(i)%b = list(k)%b
      list(k)%b = buf

      buf = list(i)%cp
      list(i)%cp = list(k)%cp
      list(k)%cp = buf

      buf = list(i)%h
      list(i)%h = list(k)%h
      list(k)%h = buf

      buf = list(i)%dh
      list(i)%dh = list(k)%dh
      list(k)%dh = buf

    end do
    return
  end subroutine Sort_listTP_Tsmall2big
!===============================================================================
!===============================================================================
!> \brief Check the monotonicity of the $\rho h$ along $h$ and $T$
!>
!> This subroutine is called locally once building up the thermal property
!> relationships. Non-monotonicity could happen in fluids at supercritical
!> pressure when inproper reference temperature is given.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  none          NA
!_______________________________________________________________________________
  subroutine Check_monotonicity_DH_of_HT_list
    use parameters_constant_mod, only : MINP
    integer :: i
    real(WP) :: ddh1, dt1, dh1
    real(WP) :: ddh2, dt2, dh2
    real(WP) :: ddt, ddh

    do i = 2, nlist - 1
        ddh1 = listTP(i)%dh - listTP(i - 1)%dh
        dt1 = listTP(i)%t - listTP(i - 1)%t
        dh1 = listTP(i)%h - listTP(i - 1)%h

        ddh2 = listTP(i + 1)%dh - listTP(i)%dh
        dt2 = listTP(i + 1)%t - listTP(i)%t
        dh2 = listTP(i + 1)%h - listTP(i)%h

        ddt = ddh1 / dt1 * ddh2 / dt2 
        ddh = ddh1 / dh1 * ddh2 / dh2

        if (ddt < MINP) STOP 'Error. The relation (rho * h) = FUNCTION (T) is not monotonicity.'
        if (ddh < MINP) STOP 'Error. The relation (rho * h) = FUNCTION (H) is not monotonicity.'

    end do
    return
  end subroutine Check_monotonicity_DH_of_HT_list
!===============================================================================
!> \brief Building up the thermal property relations from the given table.
!>
!> This subroutine is called once after reading the table.
!> [mpi] all ranks
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]  none          NA
!===============================================================================
  subroutine Buildup_property_relations_from_table(t0Ref)
    use mpi_mod
    use iso_fortran_env, only : ERROR_UNIT, IOSTAT_END
    use parameters_constant_mod, only : ZERO
    real(WP), intent(in) :: t0Ref

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    character(len = 80) :: str
    real(WP) :: rtmp
    integer :: i
    logical :: is_dim
!-------------------------------------------------------------------------------
! to read given table of thermal properties
!-------------------------------------------------------------------------------
    open ( newunit = inputUnit,     &
           file    = inputProperty, &
           status  = 'old',         &
           action  = 'read',        &
           iostat  = ioerr,         &
           iomsg   = iotxt)
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', inputProperty, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 4
    end if

    nlist = 0
    read(inputUnit, *, iostat = ioerr) str
    do
      read(inputUnit, *, iostat = ioerr) rtmp, rtmp, rtmp, rtmp, &
      rtmp, rtmp, rtmp, rtmp
      if(ioerr /= 0) exit
      nlist = nlist + 1
    end do
    rewind(inputUnit)
!-------------------------------------------------------------------------------
! to read given table of thermal properties
!-------------------------------------------------------------------------------
    allocate ( listTP (nlist) )

    read(inputUnit, *, iostat = ioerr) str
    block_tablereading: do i = 1, nlist
      call listTP(i)%Get_initialized_thermal_properties()
      read(inputUnit, *, iostat = ioerr) rtmp, listTP(i)%h, listTP(i)%t, listTP(i)%d, &
      listTP(i)%m, listTP(i)%k, listTP(i)%cp, listTP(i)%b
      listTP(i)%dh = listTP(i)%d * listTP(i)%h
    end do block_tablereading
    close(inputUnit)
!-------------------------------------------------------------------------------
! to sort input date based on Temperature (small to big)
!-------------------------------------------------------------------------------
    call Sort_listTP_Tsmall2big ( listTP(:) )

    ! to update reference of thermal properties
    tpRef0%t = t0Ref
    call tpRef0%Refresh_thermal_properties_from_T_dimensional

    ! to unify/undimensionalize the table of thermal property
    listTP(:)%t = listTP(:)%t / tpRef0%t
    listTP(:)%d = listTP(:)%d / tpRef0%d
    listTP(:)%m = listTP(:)%m / tpRef0%m
    listTP(:)%k = listTP(:)%k / tpRef0%k
    listTP(:)%b = listTP(:)%b / tpRef0%b
    listTP(:)%cp = listTP(:)%cp / tpRef0%cp
    listTP(:)%h = (listTP(:)%h - tpRef0%h) / tpRef0%t / tpRef0%cp
    listTP(:)%dh = listTP(:)%d * listTP(:)%h
    return
  end subroutine Buildup_property_relations_from_table
!===============================================================================

!> \brief Building up the thermal property relations from defined relations.
!>
!> This subroutine is called once after defining the relations.
!>
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]  none          NA
!===============================================================================
  subroutine Buildup_property_relations_from_function (T0Ref)
    implicit none
    real(WP), intent(in) :: T0Ref
    integer :: i

    ! to update reference of thermal properties
    tpRef0%t = T0Ref
    call tpRef0%Refresh_thermal_properties_from_T_dimensional

    
    nlist = 1024
    allocate ( listTP (nlist) )
    do i = 1, nlist
      call listTP(i)%Get_initialized_thermal_properties()
      listTP(i)%t = ( Tm0 + (Tb0 - Tm0) * real(i, WP) / real(nlist, WP) ) / tpRef0%t
      call listTP(i)%Refresh_thermal_properties_from_T_undim
    end do
    return
  end subroutine Buildup_property_relations_from_function
!===============================================================================
!===============================================================================
!> \brief Write out the rebuilt thermal property relations.
!>
!> This subroutine is called for testing.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  none          NA
!_______________________________________________________________________________
  subroutine Write_thermo_property
    !use mpi_mod ! for test
    use parameters_constant_mod, only : ZERO, TRUNCERR
    implicit none
    type(thermoProperty_t) :: tp
    integer :: n, i
    real(WP) :: dhmax1, dhmin1
    real(WP) :: dhmax, dhmin
    integer :: tp_unit
    logical :: is_dim

    !if (nrank /= 0) return

    n = 128
    call tp%Get_initialized_thermal_properties

    dhmin = ZERO
    dhmax = ZERO
    if(ipropertyState == IPROPERTY_TABLE) then 
      dhmin = listTP(1)%dh + TRUNCERR
      dhmax = listTP(nlist)%dh - TRUNCERR
    end if

    if(ipropertyState == IPROPERTY_FUNCS) then 
      is_dim = .false.
      
      tp%t  = TB0 / tpRef0%t
      call tp%Refresh_thermal_properties_from_T_undim
      dhmin1 = tp%dh

      tp%t  = TM0 / tpRef0%t
      call tp%Refresh_thermal_properties_from_T_undim
      dhmax1 = tp%dh
      
      dhmin = dmin1( dhmin1, dhmax1) + TRUNCERR
      dhmax = dmax1( dhmin1, dhmax1) - TRUNCERR
    end if

    open (newunit = tp_unit, file = 'check_tp_from_dh.dat')
    do i = 1, n
      tp%dh = dhmin + (dhmax - dhmin) * real(i - 1, WP) / real(n - 1, WP)
      call tp%Refresh_thermal_properties_from_DH()
      call tp%is_T_in_scope()
      write(tp_unit, '(dt)') tp
    end do
    close (tp_unit)
    return
  end subroutine Write_thermo_property
!===============================================================================
!> \brief Identify table or equations for thermal properties based on input
!>  fluid material.
!>
!> This subroutine is called once in setting up thermal relations.
!> [mpi] all ranks
!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!  mode           name          role                                           
!-------------------------------------------------------------------------------
!> \param[inout]  none          NA
!===============================================================================
  subroutine Initialize_thermo_parameters(ifluid)
    use mpi_mod
    integer, intent(in) :: ifluid
    implicit none

    if(nrank == 0) call Print_debug_start_msg("Initializing thermal parameters ...")
    select case (ifluid)
    case (ISCP_WATER)
      ipropertyState = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_WATER)

    case (ISCP_CO2)
      ipropertyState = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_CO2)

    case (ILIQUID_SODIUM)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_Na
      TB0 = TB0_Na
      HM0 = HM0_Na
      CoD(0:1) = CoD_Na(0:1)
      CoK(0:2) = CoK_Na(0:2)
      CoB = CoB_Na
      CoCp(-2:2) = CoCp_Na(-2:2)
      CoH(-1:3) = CoH_Na(-1:3)
      CoM(-1:1) = CoM_Na(-1:1)

    case (ILIQUID_LEAD)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_Pb
      TB0 = TB0_Pb
      HM0 = HM0_Pb
      CoD(0:1) = CoD_Pb(0:1)
      CoK(0:2) = CoK_Pb(0:2)
      CoB = CoB_Pb
      CoCp(-2:2) = CoCp_Pb(-2:2)
      CoH(-1:3) = CoH_Pb(-1:3)
      CoM(-1:1) = CoM_Pb(-1:1)

    case (ILIQUID_BISMUTH)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_BI
      TB0 = TB0_BI
      HM0 = HM0_BI
      CoD(0:1) = CoD_BI(0:1)
      CoK(0:2) = CoK_BI(0:2)
      CoB = CoB_BI
      CoCp(-2:2) = CoCp_BI(-2:2)
      CoH(-1:3) = CoH_BI(-1:3)
      CoM(-1:1) = CoM_BI(-1:1)

    case (ILIQUID_LBE)
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_LBE
      TB0 = TB0_LBE
      HM0 = HM0_LBE
      CoD(0:1) = CoD_LBE(0:1)
      CoK(0:2) = CoK_LBE(0:2)
      CoB = CoB_LBE
      CoCp(-2:2) = CoCp_LBE(-2:2)
      CoH(-1:3) = CoH_LBE(-1:3)
      CoM(-1:1) = CoM_LBE(-1:1)
    case default
      ipropertyState = IPROPERTY_FUNCS
      TM0 = TM0_Na
      TB0 = TB0_Na
      HM0 = HM0_Na
      CoD(0:1) = CoD_Na(0:1)
      CoK(0:2) = CoK_Na(0:2)
      CoB = CoB_Na
      CoCp(-2:2) = CoCp_Na(-2:2)
      CoH(-1:3) = CoH_Na(-1:3)
      CoM(-1:1) = CoM_Na(-1:1)

    end select


    if(nrank == 0) call Print_debug_end_msg
    return
  end subroutine Initialize_thermo_parameters

!===============================================================================
!===============================================================================
!> \brief The main code for thermal property initialisation.
!>
!> This subroutine is called once in \ref Initialize_chapsim.
!> [mpi] all ranks
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  none          NA
!_______________________________________________________________________________
  subroutine Buildup_thermo_mapping_relations
    use input_general_mod, only : is_any_energyeq
    use var_dft_mod, only : thermo
    implicit none
    
    if ( .not. is_any_energyeq) return

    call Initialize_thermo_parameters(thermo(1)%ifluid)
    if (ipropertyState == IPROPERTY_TABLE) call Buildup_property_relations_from_table(thermo(1)%T0Ref)
    if (ipropertyState == IPROPERTY_FUNCS) call Buildup_property_relations_from_function(thermo(1)%T0Ref)
    call Check_monotonicity_DH_of_HT_list
    call Write_thermo_property ! for test
    return
  end subroutine Buildup_thermo_mapping_relations
  
end module input_thermo_mod



