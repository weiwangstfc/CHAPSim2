module input_thermo_mod
  use precision_mod
  use input_general_mod, only: t0ref, tiref
  implicit none

  integer, parameter :: ISCP_WATER      = 1, &
                        ISCP_CO2        = 2, &
                        ILIQUID_SODIUM  = 3, &
                        ILIQUID_LEAD    = 4, &
                        ILIQUID_BISMUTH = 5, &
                        ILIQUID_LBE     = 6

  integer, parameter :: IPROPERTY_TABLE = 1, &
                        IPROPERTY_FUNCS = 2

  character(len = 64), parameter :: INPUT_SCP_WATER = 'NIST_WATER_23.5MP.DAT'
  character(len = 64), parameter :: INPUT_SCP_CO2   = 'NIST_CO2_8MP.DAT'
  character(len = 64) :: inputProperty

  integer :: iproperty

  real(WP), parameter :: TM0_Na = 371.0 ! unit: K, melting temperature at 1 atm for Na
  real(WP), parameter :: TM0_Pb = 600.6 ! unit: K, melting temperature at 1 atm for Lead
  real(WP), parameter :: TM0_BI = 544.6 ! unit: K, melting temperature at 1 atm for Bismuth
  real(WP), parameter :: TM0_LBE = 398.0 ! unit: K, melting temperature at 1 atm for LBE
  real(WP) :: TM0

  real(WP), parameter :: TB0_Na = 1155.0 ! unit: K, boling temperature at 1 atm for Na
  real(WP), parameter :: TB0_Pb = 2021.0 ! unit: K, boling temperature at 1 atm for Lead
  real(WP), parameter :: TB0_BI = 1831.0 ! unit: K, boling temperature at 1 atm for Bismuth
  real(WP), parameter :: TB0_LBE = 1927.0 ! unit: K, boling temperature at 1 atm for LBE
  real(WP) :: TB0

  real(WP), parameter :: HM0_Na = 113.0e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP), parameter :: HM0_Pb = 23.07e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP), parameter :: HM0_BI = 53.3e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP), parameter :: HM0_LBE = 38.6e3 ! unit: J / Kg, latent melting heat, enthalpy
  real(WP) :: HM0
  ! D = CoD(0) + CoD(1) * T
  real(WP), parameter :: CoD_Na(0:1) = (/1014.0, -0.235/)
  real(WP), parameter :: CoD_Pb(0:1) = (/11441.0, -1.2795/)
  real(WP), parameter :: CoD_Bi(0:1) = (/10725.0, -1.22 /)
  real(WP), parameter :: CoD_LBE(0:1) = (/11065.0, 1.293 /)
  real(WP) :: CoD(0:1)
  ! K = CoK(0) + CoK(1) * T + CoK(2) * T^2
  real(WP), parameter :: CoK_Na(0:2) = (/104.0, -0.047, 0.0/)
  real(WP), parameter :: CoK_Pb(0:2) = (/9.2, 0.011, 0.0/)
  real(WP), parameter :: CoK_Bi(0:2) = (/7.34, 9.5E-3, 0.0/)
  real(WP), parameter :: CoK_LBE(0:2) = (/ 3.284, 1.617E-2, -2.305E-6/)
  real(WP) :: CoK(0:2)
  ! B = 1 / (CoB - T)
  real(WP), parameter :: CoB_Na = 4316.0
  real(WP), parameter :: CoB_Pb = 8942.0
  real(WP), parameter :: CoB_BI = 8791.0
  real(WP), parameter :: CoB_LBE = 8558.0
  real(WP) :: CoB
  ! Cp = CoCp(-2) * T^(-2) + CoCp(-1) * T^(-1) + CoCp(0) + CoCp(1) * T + CoCp(2) * T^2
  real(WP), parameter :: CoCp_Na(-2:2) = (/- 3.001e6, 0.0, 1658.0, -0.8479, 4.454E-4/)
  real(WP), parameter :: CoCp_Pb(-2:2) = (/- 1.524e6, 0.0, 176.2, -4.923E-2, 1.544E-5/)
  real(WP), parameter :: CoCp_Bi(-2:2) = (/7.183e6, 0.0, 118.2, 5.934E-3, 0.0/)
  real(WP), parameter :: CoCp_LBE(-2:2) = (/-4.56e5, 0.0, 164.8, - 3.94E-2, 1.25E-5/)
  real(WP) :: CoCp(-2:2)
  ! H = HM0 + CoH(-1) * (1 / T - 1 / Tm0) + CoH(0) + CoH(1) * (T - Tm0) +  CoH(2) * (T^2 - Tm0^2) +  CoH(3) * (T^3- Tm0^3)
  real(WP), parameter :: CoH_Na(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)
  real(WP), parameter :: CoH_Pb(-1:3) = (/1.524e6, 0.0, 176.2, -2.4615E-2, 5.147E-6/)
  real(WP), parameter :: CoH_Bi(-1:3) = (/-7.183e6, 0.0, 118.2, 2.967E-3, 0.0/)
  real(WP), parameter :: CoH_LBE(-1:3) = (/4.56e5, 0.0, 164.8, -1.97E-2, 4.167E-4/)! check, WRong from literature.
  real(WP) :: CoH(-1:3)
  ! M = vARies
  real(WP), parameter :: CoM_Na(-1:1) = (/556.835, -6.4406, -0.3958/) ! M = exp ( CoM(-1) / T + CoM(0) + CoM(1) * ln(T) )
  real(WP), parameter :: CoM_Pb(-1:1) = (/1069.0, 4.55E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_Bi(-1:1) = (/780.0, 4.456E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP), parameter :: CoM_LBE(-1:1) = (/754.1, 4.94E-4, 0.0/) ! M = CoM(0) * exp (CoM(-1) / T)
  real(WP) :: CoM(-1:1)

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

  private :: Initialize_thermo_parameters
  private :: Building_property_relations_from_table
  private :: Building_property_relations_from_equation
  private :: Prepare_function_DH_of_T

  public  :: Initialize_thermo_input
  public  :: Calculate_thermal_variables_from_temperature
  


contains
  !!--------------------------------
  subroutine Initialize_thermo_input
    use input_general_mod, only : ithermo
    if (ithermo /= 1) return 
    call Initialize_thermo_parameters
    if (iproperty == IPROPERTY_TABLE) call Building_property_relations_from_table
    if (iproperty == IPROPERTY_FUNCS) call Building_property_relations_from_equation

  end subroutine Initialize_thermo_input

  !!--------------------------------
  subroutine Initialize_thermo_parameters
    use input_general_mod, only: ifluid

    select case (ifluid)
    case (ISCP_WATER)
      iproperty = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_WATER)

    case (ISCP_CO2)
      iproperty = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_CO2)

    case (ILIQUID_SODIUM)
      iproperty = IPROPERTY_FUNCS
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
      iproperty = IPROPERTY_FUNCS
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
      iproperty = IPROPERTY_FUNCS
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
      iproperty = IPROPERTY_FUNCS
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
      iproperty = IPROPERTY_FUNCS
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


  end subroutine Initialize_thermo_parameters
  !!--------------------------------
  subroutine Building_property_relations_from_table ( )
    use iso_fortran_env, only : ERROR_UNIT, IOSTAT_END
    use parameters_constant_mod, only : ZERO
    use table_index_locating_mod

    integer, parameter :: IOMSG_LEN = 200
    character(len = IOMSG_LEN) :: iotxt
    integer :: ioerr, inputUnit
    character(len = 80) :: str
    real(WP) :: rtmp
    real(WP) :: bufT, bufH, bufD, bufM, bufK, bufB, bufCp
    integer :: i, k
 
    open ( newunit = inputUnit, file = inputProperty, status = 'old', action  = 'read', &
          iostat = ioerr, iomsg = iotxt)
    if(ioerr /= 0) then
      write (ERROR_UNIT, *) 'Problem openning : ', inputProperty, ' for reading.'
      write (ERROR_UNIT, *) 'Message: ', trim (iotxt)
      stop 4
    end if

    read(inputUnit, *, iostat = ioerr) str
    read(inputUnit, *, iostat = ioerr) nlist

    allocate ( listH (nlist) ); listH(:) = ZERO
    allocate ( listT (nlist) ); listT(:) = ZERO
    allocate ( listD (nlist) ); listD(:) = ZERO
    allocate ( listM (nlist) ); listM(:) = ZERO
    allocate ( listK (nlist) ); listK(:) = ZERO
    allocate ( listB (nlist) ); listB(:) = ZERO
    allocate ( listCp(nlist) ); listCp(:) = ZERO
    allocate ( listDH(nlist) ); listDH(:) = ZERO

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

  end subroutine Building_property_relations_from_table

  !!--------------------------------
  subroutine Building_property_relations_from_equation ( )
    
    d0ref  = Calculate_thermal_variables_from_temperature (t0ref, 'TD', 1)
    m0ref  = Calculate_thermal_variables_from_temperature (t0ref, 'TM', 1)
    k0ref  = Calculate_thermal_variables_from_temperature (t0ref, 'TK', 1)
    cp0ref = Calculate_thermal_variables_from_temperature (t0ref, 'TCp', 1)
    b0ref  = Calculate_thermal_variables_from_temperature (t0ref, 'TB', 1)
    h0ref  = Calculate_thermal_variables_from_temperature (t0ref, 'TH', 1)

    call Prepare_function_DH_of_T()

  end subroutine Building_property_relations_from_equation

  !!--------------------------------  
  pure function Calculate_thermal_variables_from_temperature(t, ft, dim) result(d)
    use parameters_constant_mod, only: ONE
    use input_general_mod, only: ifluid
    real(WP), intent(in) :: t        
    character(*), intent(in) :: ft
    integer, intent(in), optional :: dim ! without = undim, with = dim.
    real(WP) :: d1, d
    real(WP) :: t1

    if (present(dim)) then 
      t1 = t
    else 
      ! convert undim to dim 
      t1 = t * t0ref
    end if

    d = ONE
    ! entropy = f(T)
    block_TH: if (trim(ft) == 'TH') then
      d1 = Hm0 + &
        CoH(-1) * (ONE / t1 - ONE / Tm0) + &
        CoH(0) + &
        CoH(1) * (t1 - Tm0) + &
        CoH(2) * (t1**2 - Tm0**2) + &
        CoH(3) * (t1**3 - Tm0**3)
      if (present(dim)) then 
        d = d1
      else 
        d = (d1 - h0ref) / (cp0ref * t0ref)
      end if
    end if block_TH

    ! density = f(T)
    block_TD: if (trim(ft) == 'TD') then
      d1 = CoD(0) + CoD(1) * t1
      if (present(dim)) then 
        d = d1
      else 
        d = d1 / d0ref
      end if
    end if block_TD

    ! thermal conductivity = f(T)
    block_TK: if (trim(ft) == 'TK') then
      d1 = CoK(0) + CoK(1) * t1 + CoK(2) * t1**2
      if (present(dim)) then 
        d = d1
      else 
        d = d1 / k0ref
      end if

    end if block_TK

    ! dynamic viscosity = f(T)
    block_TM: if (trim(ft) == 'TM') then
      select case (ifluid)
        ! unit: T(Kelvin), M(Pa S)
        case (ILIQUID_SODIUM)
          d1 = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
        case (ILIQUID_LEAD)
          d1 = CoM_Pb(0) * EXP (CoM_Pb(-1) / t1)
        case (ILIQUID_BISMUTH)
          d1 = CoM_Bi(0) * EXP (CoM_Bi(-1) / t1)
        case (ILIQUID_LBE)
          d1 = CoM_LBE(0) * EXP (CoM_LBE(-1) / t1)
        case default
          d1 = EXP( CoM_Na(-1) / t1 + CoM_Na(0) + CoM_Na(1) * LOG(t1) )
      end select
      if (present(dim)) then 
        d = d1
      else 
        d = d1 / m0ref
      end if
    
    end if block_TM

    ! Cp = f(T)
    block_TCp: if (trim(ft) == 'TCp') then
      d1 = CoCp(-2) * t1**(-2) + CoCp(-1) * t1**(-1) + CoCp(0) + CoCp(1) * t1 + CoCp(2) * t1**2
      if (present(dim)) then 
        d = d1
      else 
        d = d1 / cp0ref
      end if

    end if block_TCp
    
    ! B = f(T)
    block_TB: if (trim(ft) == 'TB') then
      d1 = ONE / (CoB - t1)
      if (present(dim)) then 
        d = d1
      else 
        d = d1 / b0ref
     end if
    end if block_TB

  end function Calculate_thermal_variables_from_temperature

  subroutine Prepare_function_DH_of_T
    use parameters_constant_mod, only: ZERO
    use sort_array_mod
    integer :: i
    real(WP) :: t, d, h

    nlist = 1024
    allocate ( listH (nlist) ); listH(:) = ZERO
    allocate ( listT (nlist) ); listT(:) = ZERO
    allocate ( listDH (nlist) ); listDH(:) = ZERO

    do i = 1, nlist
      t = ( Tm0 + (Tb0 - Tm0) * real(i, WP) / real(nlist, WP) ) / t0ref
      d = Calculate_thermal_variables_from_temperature(t, 'TD')
      h = Calculate_thermal_variables_from_temperature(t, 'TH')

      listT(i) = t
      listH(i) = h
      listDH(i) = d * h 
    end do

    call Sort_array_small2big(listT(:), listH(:), listDH(:))
    
  end if


  end subroutine Prepare_function_DH_of_T
  

end module input_thermo_mod

