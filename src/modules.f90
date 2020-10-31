!##############################################################################
module mpi_info_mod
  include "mpif.h"
  integer :: ierror
  integer :: nrank
  integer :: nproc
end module mpi_info_mod

!##############################################################################
module precision_mod
  use mpi_info_mod
  integer, parameter :: I4 = selected_int_kind( 4 )
  integer, parameter :: I8 = selected_int_kind( 8 )
  integer, parameter :: I15 = selected_int_kind( 15 )
  integer, parameter :: SP = selected_real_kind( p = 6, r = 37 )
  integer, parameter :: DP = selected_real_kind( p = 15, r = 307 )
  integer, parameter :: QP = selected_real_kind( p = 33, r = 4931 )

  integer, parameter :: WP = DP

end module precision_mod
!##############################################################################
module parameters_input_mod
  use precision_mod

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
  integer :: npx, npy, npz
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

end module parameters_input_mod

!##############################################################################
module parameters_constant_mod
  use precision_mod

  real(WP), parameter :: ZPONE    = 0.1_WP
  real(WP), parameter :: ZPTWO    = 0.2_WP
  real(WP), parameter :: ZPTHREE  = 0.3_WP
  real(WP), parameter :: ZPFOUR   = 0.4_WP
  real(WP), parameter :: ZPFIVE   = 0.5_WP
  real(WP), parameter :: ZPSIX    = 0.6_WP
  real(WP), parameter :: ZPSEVEN  = 0.7_WP
  real(WP), parameter :: ZPEIGHT  = 0.8_WP
  real(WP), parameter :: ZPNINE   = 0.9_WP

  real(WP), parameter :: HALF     = 0.5_WP
  real(WP), parameter :: ZERO     = 0.0_WP
  real(WP), parameter :: ONE      = 1.0_WP
  real(WP), parameter :: ONEPFIVE = 1.5_WP
  real(WP), parameter :: TWO      = 2.0_WP
  real(WP), parameter :: THREE    = 3.0_WP
  real(WP), parameter :: FOUR     = 4.0_WP
  real(WP), parameter :: FIVE     = 5.0_WP
  real(WP), parameter :: SIX      = 6.0_WP
  real(WP), parameter :: SEVEN    = 7.0_WP
  real(WP), parameter :: EIGHT    = 8.0_WP
  real(WP), parameter :: NINE     = 9.0_WP

  real(WP),parameter :: PI = dacos( -ONE )
  real(WP),parameter :: TWOPI = TWO * dacos( -ONE )


  integer, parameter :: ICASE_CHANNEL = 1, &
                        ICASE_PIPE    = 2, &
                        ICASE_ANNUAL  = 3, &
                        ICASE_TGV     = 4
                        
  integer, parameter :: ISTRET_NO     = 0, &
                        ISTRET_CENTRE = 1, &
                        ISTRET_SIDES  = 2, &
                        ISTRET_BOTTOM = 3, &
                        ISTRET_TOP    = 4

  integer, parameter :: ISCP_WATER      = 1, &
                        ISCP_CO2        = 2, &
                        ILIQUID_SODIUM  = 3, &
                        ILIQUID_LEAD    = 4, &
                        ILIQUID_BISMUTH = 5, &
                        ILIQUID_LBE     = 6

end module parameters_constant_mod

!##############################################################################
module parameters_properties_mod
  use precision_mod

end module parameters_properties_mod


