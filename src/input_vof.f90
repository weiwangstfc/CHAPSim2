!----------------------------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------------------------
!==========================================================================================================
!> \file input_vof.f90
!>
!> \brief initialise and update properties within the interface
!>
!==========================================================================================================
module vof_info_mod
  use parameters_constant_mod
  use udf_type_mod
  implicit none

  public  :: Initialize_vof_properties
  public  :: Update_vof_properties

contains

!==========================================================================================================
!==========================================================================================================
  subroutine Initialize_vof_properties (fl, vf)
    use parameters_constant_mod
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_vof),    intent(inout) :: vf

    if(nrank == 0) call Print_debug_mid_msg("Initialize vof interface properties ...")

    !fl%dDens(:, :, :) = (ONE-vf%phi(:, :, :))*vf%rho1 + vf%phi(:, :, :)*vf%rho2
    !fl%mVisc(:, :, :) = (ONE-vf%phi(:, :, :))*vf%mu1 + vf%phi(:, :, :)*vf%mu2

    !fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    !fl%dDensm2(:, :, :) = fl%dDensm1(:, :, :)

    if(vf%rho1.le.vf%rho2) then
      vf%rho0 = vf%rho1
      vf%mu0  = vf%mu1
    else
      vf%rho0 = vf%rho2
      vf%mu0  = vf%mu2
    end if

    vf%rhostar(:, :, :) = (ONE-vf%phi(:, :, :))*vf%rho1 + vf%phi(:, :, :)*vf%rho2
    vf%mustar(:, :, :)  = (ONE-vf%phi(:, :, :))*vf%mu1 + vf%phi(:, :, :)*vf%mu2
    vf%rhostar(:, :, :) = vf%rhostar(:, :, :)/vf%rho0
    vf%mustar(:, :, :)  = vf%mustar(:, :, :)/vf%mu0

    if(nrank == 0) call Print_debug_end_msg

    return
  end subroutine Initialize_vof_properties

!==========================================================================================================
!==========================================================================================================
  subroutine Update_vof_properties(dm, fl, vf)
    use udf_type_mod
    implicit none
    type(t_domain), intent(in) :: dm
    type(t_flow),   intent(inout) :: fl
    type(t_vof),    intent(inout) :: vf

!----------------------------------------------------------------------------------------------------------
!   x-pencil
!----------------------------------------------------------------------------------------------------------
    !fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    !fl%dDensm2(:, :, :) = fl%dDensm1(:, :, :)

    !fl%dDens(:, :, :) = (ONE-vf%phi(:, :, :))*vf%rho1 + vf%phi(:, :, :)*vf%rho2
    !fl%mVisc(:, :, :) = (ONE-vf%phi(:, :, :))*vf%mu1 + vf%phi(:, :, :)*vf%mu2

    ! putting this here will ignore drho/dt in the continuity equation
    !fl%dDensm1(:, :, :) = fl%dDens(:, :, :)
    !fl%dDensm2(:, :, :) = fl%dDensm1(:, :, :)

    vf%rhostar(:, :, :) = (ONE-vf%phi(:, :, :))*vf%rho1 + vf%phi(:, :, :)*vf%rho2
    vf%mustar(:, :, :)  = (ONE-vf%phi(:, :, :))*vf%mu1 + vf%phi(:, :, :)*vf%mu2
    vf%rhostar(:, :, :) = vf%rhostar(:, :, :)/vf%rho0
    vf%mustar(:, :, :)  = vf%mustar(:, :, :)/vf%mu0

  return
  end subroutine Update_vof_properties


end module vof_info_mod
