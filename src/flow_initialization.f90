module poiseuille_flow_initialization
  implicit none
  private 

  private :: Generate_poiseuille_flow_profile
  public  :: Initialize_poiseuille_flow

contains
  subroutine Generate_poiseuille_flow_profile(u_laminar, domain_dummy, cell_dummy)
    use input_general_mod, only: lyt, lyb, icase
    use chapsim_abort_mod
    real(WP), intent(out) :: u_laminar(:)
    type(domain_t), intent(in) :: domain_dummy
    type(cell_t), intent(in) :: cell_dummy(:, :, :)
    
    real(WP) :: a, b, c, yy

    u_laminar (:) = ZERO

    if (icase == ICASE_CHANNEL) then
      a = (lyt - lyb) / TWO
      b = ZERO
      c = ONEPFIVE
    else if (icase == ICASE_PIPE) then
      a = lyt - lyb
      b = ZERO
      c = TWO
    else if (icase == ICASE_ANNUAL) then
      a = (lyt - lyb) / TWO
      b = (lyt + lyb) / TWO
      c = TWO
    else 
      a = MAXP
      b = ZERO
      c = ONE
    end if

    do j = 1, domain_dummy%nc(2)
      yy = cell_dummy(1, j, 1)%y
      u_laminar(j) = ( ONE - ( (yy - b)**2 ) / a / a ) * c
    end do

    ! check the mean velocity
    umean = ZERO
    do j = 1, domain_dummy%nc(2)
      umean = umean + u_laminar(j) * cell_dummy(1, j, 1)%dy
    end do
    umean = umean / (lyt - lyb)
    if ( abs_wp(umean - ONE) > TRUNCERR) &
    call Print_error_msg(2, "Error in poiseuille_flow_profile.")

  end subroutine Generate_poiseuille_flow_profile

  subroutine Initialize_poiseuille_flow(flow_dummy, domain_dummy, node_dummy, cell_dummy)
    use random_number_generation_mod
    type(flow_t), intent(out) :: flow_dummy(:, :, :)
    type(domain_t), intent(in) :: domain_dummy
    type(cell_t), intent(in) :: cell_dummy(:, :, :)
    type(node_t), intent(in) :: node_dummy(:, :, :)
    
    real(WP), allocatable, dimension(:) :: u_laminar
    integer :: i, j, k
    integer :: seed
    real(WP) :: rd(3)

    allocate (u_laminar (domain_dummy%nc) )
    u_laminar(:) = ZERO
    call Generate_poiseuille_flow_profile ( u_laminar, domain_dummy, cell_dummy(1, :, 1) )

    flow_dummy(:, :, :)%pre = ZERO
    flow_dummy(:, :, :)%phi = ZERO

    flow_dummy(:, :, :)%u(:) = ZERO
    flow_dummy(:, :, :)%g(:) = ZERO

    seed = 0 ! real random
    do k = domain_dummy%krange(1), domain_dummy%krange(2)
      do j = domain_dummy%jrange(1), domain_dummy%jrange(2)
        do i = domain_dummy%irange(1), domain_dummy%irange(2)
          seed = seed + k + j + i ! repeatable random
          call Initialize_random_number ( seed )
          call Generate_rvec_random ( -ONE, ONE, 3, rd)

          flow_dummy(i, j, k)%u(1) = initNoise * rd(1) + u_laminar (j)
          flow_dummy(i, j, k)%u(2) = initNoise * rd(2) / node_dummy(i, j, k)%ri
          flow_dummy(i, j, k)%u(3) = initNoise * rd(3) / cell_dummy(i, j, k)%ri
        end do
      end do
    end do
    deallocate (u_laminar)
  end subroutine  Initialize_poiseuille_flow

end module  poiseuille_flow_initialization



module vortexgreen_flow_initialization
  implicit none
  private 
  public  :: Initialize_vortexgreen_flow

contains

  subroutine  Initialize_poiseuille_flow(flow_dummy, domain_dummy, node_dummy, cell_dummy)
    type(flow_t), intent(out) :: flow_dummy(:, :, :)
    type(domain_t), intent(in) :: domain_dummy
    type(cell_t), intent(in) :: cell_dummy(:, :, :)
    type(node_t), intent(in) :: node_dummy(:, :, :)

    flow_dummy(:, :, :)%pre = ZERO
    flow_dummy(:, :, :)%phi = ZERO

    flow_dummy(:, :, :)%u(:) = ZERO
    flow_dummy(:, :, :)%g(:) = ZERO

    do k = domain_dummy%krange(1), domain_dummy%krange(2)
      do j = domain_dummy%jrange(1), domain_dummy%jrange(2)
        do i = domain_dummy%irange(1), domain_dummy%irange(2)
          flow_dummy(i, j, k)%u(1) =  sin_wp ( node_dummy(i, j, k)%x ) * &
                                      cos_wp ( cell_dummy(i, j, k)%y ) * &
                                      cos_wp ( cell_dummy(i, j, k)%z )
          flow_dummy(i, j, k)%u(2) = -cos_wp ( cell_dummy(i, j, k)%x ) * &
                                      sin_wp ( node_dummy(i, j, k)%y ) * &
                                      cos_wp ( cell_dummy(i, j, k)%z )
          flow_dummy(i, j, k)%u(3) = ZERO
          flow_dummy(i, j, k)%pre = ( cos_wp( TWO * cell_dummy(i, j, k)%x       ) + &
                                      cos_wp( TWO * cell_dummy(i, j, k)%y       ) ) * &
                                    ( cos_wp( TWO * cell_dummy(i, j, k)%z + TWO ) ) / SIXTEEN
        end do
      end do
    end do


  end subroutine Initialize_vortexgreen_flow

end module