module input_thermo_mod
  use input_mod, only : ifluid, igravity, lenRef, t0Ref, tiRef
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

  



contains

  subroutine Initialize_thermo_parameters

    ! property
    select case (ifluid)
    case (ifluid == ISCP_WATER)

      iproperty = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_WATER)

    case (ifluid == ISCP_CO2)

      iproperty = IPROPERTY_TABLE
      inputProperty = TRIM(INPUT_SCP_CO2)

    case (ifluid == ILIQUID_SODIUM)

      iproperty = IPROPERTY_FUNCS

    case (ifluid == ILIQUID_LEAD)

      iproperty = IPROPERTY_FUNCS

    case (ifluid == ILIQUID_BISMUTH)

      iproperty = IPROPERTY_FUNCS

    case (ifluid == ILIQUID_LBE)

      iproperty = IPROPERTY_FUNCS
      
    case default

    end select

  end subroutine Initialize_thermo_parameters
end module input_thermo_mod