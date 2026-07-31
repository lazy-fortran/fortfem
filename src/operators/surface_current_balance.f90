module fortfem_surface_current_balance
    !! Discrete junction ledger for integrated surface currents.
    !!
    !! The integer junction incidence is topology metadata.  A column with
    !! -1 at a start boundary and +1 at an end boundary distributes the
    !! integrated current of one manifold to its boundary ledger.  Closed
    !! columns are zero, so their global balance is identically zero.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_surface_current_junction_balance
    public :: assemble_surface_current_junction_balance_jvp
    public :: assemble_surface_current_junction_balance_vjp

contains

    subroutine assemble_surface_current_junction_balance( &
            junction_incidence, manifold_current, junction_balance, &
            global_balance, status)
        integer, intent(in) :: junction_incidence(:, :)
        real(dp), intent(in) :: manifold_current(:, :)
        real(dp), intent(out) :: junction_balance(:, :), global_balance(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: junction, manifold

        junction_balance = 0.0_dp
        global_balance = 0.0_dp
        call validate_balance_inputs( &
            junction_incidence, manifold_current, junction_balance, &
            global_balance, status)
        if (status%code /= FORTSPARSE_OK) return

        do manifold = 1, size(junction_incidence, 2)
            do junction = 1, size(junction_incidence, 1)
                junction_balance(:, junction) = junction_balance(:, junction) &
                    + real(junction_incidence(junction, manifold), dp) &
                    * manifold_current(:, manifold)
            end do
        end do
        global_balance = sum(junction_balance, dim=2)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_junction_balance

    subroutine assemble_surface_current_junction_balance_jvp( &
            junction_incidence, manifold_current, manifold_current_dot, &
            junction_balance_dot, global_balance_dot, status)
        integer, intent(in) :: junction_incidence(:, :)
        real(dp), intent(in) :: manifold_current(:, :), manifold_current_dot(:, :)
        real(dp), intent(out) :: junction_balance_dot(:, :), global_balance_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: junction, manifold

        junction_balance_dot = 0.0_dp
        global_balance_dot = 0.0_dp
        call validate_balance_inputs( &
            junction_incidence, manifold_current, junction_balance_dot, &
            global_balance_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(manifold_current_dot, 1) /= 3 .or. &
            size(manifold_current_dot, 2) /= size(junction_incidence, 2) .or. &
            any(.not. ieee_is_finite(manifold_current_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current balance JVP has incompatible increments")
            return
        end if

        do manifold = 1, size(junction_incidence, 2)
            do junction = 1, size(junction_incidence, 1)
                junction_balance_dot(:, junction) = &
                    junction_balance_dot(:, junction) + &
                    real(junction_incidence(junction, manifold), dp) &
                    * manifold_current_dot(:, manifold)
            end do
        end do
        global_balance_dot = sum(junction_balance_dot, dim=2)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_junction_balance_jvp

    subroutine assemble_surface_current_junction_balance_vjp( &
            junction_incidence, manifold_current, junction_balance_bar, &
            global_balance_bar, manifold_current_bar, status)
        integer, intent(in) :: junction_incidence(:, :)
        real(dp), intent(in) :: manifold_current(:, :)
        real(dp), intent(in) :: junction_balance_bar(:, :), global_balance_bar(:)
        real(dp), intent(out) :: manifold_current_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: junction, manifold

        manifold_current_bar = 0.0_dp
        call validate_balance_inputs( &
            junction_incidence, manifold_current, junction_balance_bar, &
            global_balance_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(manifold_current_bar, 1) /= 3 .or. &
            size(manifold_current_bar, 2) /= size(junction_incidence, 2) .or. &
            any(.not. ieee_is_finite(junction_balance_bar)) .or. &
            any(.not. ieee_is_finite(global_balance_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current balance VJP has incompatible cotangents")
            return
        end if

        do manifold = 1, size(junction_incidence, 2)
            do junction = 1, size(junction_incidence, 1)
                manifold_current_bar(:, manifold) = &
                    manifold_current_bar(:, manifold) + &
                    real(junction_incidence(junction, manifold), dp) * &
                    (junction_balance_bar(:, junction) + global_balance_bar)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_junction_balance_vjp

    subroutine validate_balance_inputs( &
            junction_incidence, manifold_current, junction_balance, &
            global_balance, status)
        integer, intent(in) :: junction_incidence(:, :)
        real(dp), intent(in) :: manifold_current(:, :)
        real(dp), intent(in) :: junction_balance(:, :), global_balance(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: manifold

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface-current balance has incompatible arrays")
        if (size(manifold_current, 1) /= 3 .or. &
            size(manifold_current, 2) /= size(junction_incidence, 2) .or. &
            size(junction_balance, 1) /= 3 .or. &
            size(junction_balance, 2) /= size(junction_incidence, 1) .or. &
            size(global_balance) /= 3) return
        if (any(.not. ieee_is_finite(manifold_current))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current balance received non-finite current")
            return
        end if
        do manifold = 1, size(junction_incidence, 2)
            if (sum(junction_incidence(:, manifold)) /= 0) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_balance_inputs

end module fortfem_surface_current_balance
