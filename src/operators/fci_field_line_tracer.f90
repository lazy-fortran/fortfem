module fortfem_fci_field_line_tracer
    !! Fixed-step field-line tracing for FCI geometry services.
    !!
    !! The tracer integrates a user-supplied right-hand side with toroidal
    !! coordinate (or any other monotonically advancing parameter) `phi`.
    !! Mesh lookup, interpolation stencils, stopping conditions, and support
    !! volumes deliberately remain outside this small numerical service.
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    abstract interface
        pure subroutine fci_field_line_rhs(phi, point, derivative)
            import dp
            real(dp), intent(in) :: phi
            real(dp), intent(in) :: point(:)
            real(dp), intent(out) :: derivative(:)
        end subroutine fci_field_line_rhs
    end interface

    public :: trace_fci_field_line_rk4

contains

    subroutine trace_fci_field_line_rk4( &
            initial_point, toroidal_span, step_count, rhs, endpoint, status)
        !! Trace `d point / d phi = rhs(phi, point)` using classical RK4.
        real(dp), intent(in) :: initial_point(:)
        real(dp), intent(in) :: toroidal_span
        integer, intent(in) :: step_count
        procedure(fci_field_line_rhs) :: rhs
        real(dp), intent(out) :: endpoint(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: point(:), trial(:)
        real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:)
        real(dp) :: phi, step_size
        integer :: step

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI field-line tracing received incompatible arrays")
        if (size(initial_point) < 1) return
        if (size(endpoint) /= size(initial_point)) return
        if (step_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI field-line tracing requires a positive step count")
            return
        end if

        allocate( &
            point(size(initial_point)), trial(size(initial_point)), &
            k1(size(initial_point)), k2(size(initial_point)), &
            k3(size(initial_point)), k4(size(initial_point)))
        point = initial_point
        phi = 0.0_dp
        step_size = toroidal_span/real(step_count, dp)
        do step = 1, step_count
            call rhs(phi, point, k1)
            trial = point + 0.5_dp*step_size*k1
            call rhs(phi + 0.5_dp*step_size, trial, k2)
            trial = point + 0.5_dp*step_size*k2
            call rhs(phi + 0.5_dp*step_size, trial, k3)
            trial = point + step_size*k3
            call rhs(phi + step_size, trial, k4)
            point = point + step_size*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
            phi = phi + step_size
        end do
        endpoint = point
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine trace_fci_field_line_rk4

end module fortfem_fci_field_line_tracer
