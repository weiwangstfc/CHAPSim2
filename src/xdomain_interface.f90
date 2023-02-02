module xdomain_interface


contains
  subroutine buildup_turbgen_flow_inlet_xbc(dm0, dm1, fl0, fl1, )
    use parameters_constant_mod
    !----------------------------------------------------------------------------------------------------------
    ! default, x-pencil
    !----------------------------------------------------------------------------------------------------------
    if(.not. dm0%is_turbgen) return

    if(dm1%ibcx(1, 1) == IBC_TURBGEN) then
      fl1%qx(1, :, :)= fl0%qx(1, :, :)
      fl1%
    end if



  end subroutine


end module