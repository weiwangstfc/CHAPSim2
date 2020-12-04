module thermo_variables_mod
    use precision_mod
    use input_general_mod, only: ifluid, t0ref, tiref
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
  
    integer :: ipropertyState
  
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
  
    type thermoProperty_t
        real(WP) :: t  !temperature
        real(WP) :: d  !density
        real(WP) :: m  !dynviscosity
        real(WP) :: k  !thermconductivity
        real(WP) :: h  !enthalpy
        real(WP) :: dh ! mass enthalpy
        real(WP) :: cp ! specific heat capacity 
        real(WP) :: b  ! thermal expansion
    contains
        private
        procedure :: is_T_in_scope
        procedure :: Get_initilized_tp
        procedure :: GetTP_from_T
        procedure :: GetTP_from_H
        procedure :: GetTP_from_DH
        procedure :: Print_debug
        generic :: Print => Print_debug
        generic :: write(formatted) => Print_debug
    end type thermoProperty_t
  
    type(thermoProperty_t), save, allocatable, dimension(:) :: listTP
    type(thermoProperty_t) :: tpRef0 ! dim
    type(thermoProperty_t) :: tpIni0 ! dim
contains
    subroutine is_T_in_scope ( tp )
        class(thermoProperty_t) :: tp

        if(ipropertyState == IPROPERTY_TABLE) then
        if ( ( tp%t < listTP(1)%t     )  .OR. &
            ( tp%t > listTP(nlist)%t ) ) then
            print*, tp%t, listTP(1)%t, listTP(nlist)%t
            stop 'temperature exceeds specified range.'
        end if
        end if

        if(ipropertyState == IPROPERTY_FUNCS) then
        if ( ( tp%t < ( TM0 / tpRef0%t ) ) .OR. &
            ( tp%t > ( TB0 / tpRef0%t ) ) ) then 
            print*, tp%t, TM0 / tpRef0%t, TB0 / tpRef0%t
            stop 'temperature exceeds specified range.'
        end if
        end if

    end subroutine is_T_in_scope



end module thermo_variables_mod

