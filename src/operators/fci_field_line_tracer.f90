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

    abstract interface
        pure subroutine fci_field_line_rhs_jvp( &
                phi, point, point_dot, derivative, derivative_dot)
            import dp
            real(dp), intent(in) :: phi
            real(dp), intent(in) :: point(:), point_dot(:)
            real(dp), intent(out) :: derivative(:), derivative_dot(:)
        end subroutine fci_field_line_rhs_jvp
    end interface

    public :: trace_fci_field_line_rk4
    public :: trace_fci_field_line_rk4_jvp

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

    subroutine trace_fci_field_line_rk4_jvp( &
            initial_point, initial_point_dot, toroidal_span, step_count, rhs, &
            endpoint, endpoint_dot, status)
        !! Trace a field line and its tangent with classical RK4.
        !!
        !! The callback supplies both `d point / d phi` and its directional
        !! derivative for a supplied `point_dot`.  This is a fixed-step,
        !! matrix-free sensitivity path; topology or stopping events remain
        !! outside this service.
        real(dp), intent(in) :: initial_point(:), initial_point_dot(:)
        real(dp), intent(in) :: toroidal_span
        integer, intent(in) :: step_count
        procedure(fci_field_line_rhs_jvp) :: rhs
        real(dp), intent(out) :: endpoint(:), endpoint_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: point(:), point_dot(:), trial(:), trial_dot(:)
        real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:)
        real(dp), allocatable :: k1_dot(:), k2_dot(:), k3_dot(:), k4_dot(:)
        real(dp) :: phi, step_size
        integer :: step

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI field-line tracer JVP received incompatible arrays")
        endpoint = 0.0_dp
        endpoint_dot = 0.0_dp
        if (size(initial_point) < 1) return
        if (size(initial_point_dot) /= size(initial_point)) return
        if (size(endpoint) /= size(initial_point)) return
        if (size(endpoint_dot) /= size(initial_point)) return
        if (step_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI field-line tracer JVP requires a positive step count")
            return
        end if

        allocate( &
            point(size(initial_point)), point_dot(size(initial_point)), &
            trial(size(initial_point)), trial_dot(size(initial_point)), &
            k1(size(initial_point)), k2(size(initial_point)), &
            k3(size(initial_point)), k4(size(initial_point)), &
            k1_dot(size(initial_point)), k2_dot(size(initial_point)), &
            k3_dot(size(initial_point)), k4_dot(size(initial_point)))
        point = initial_point
        point_dot = initial_point_dot
        phi = 0.0_dp
        step_size = toroidal_span/real(step_count, dp)
        do step = 1, step_count
            call rhs(phi, point, point_dot, k1, k1_dot)
            trial = point + 0.5_dp*step_size*k1
            trial_dot = point_dot + 0.5_dp*step_size*k1_dot
            call rhs(phi + 0.5_dp*step_size, trial, trial_dot, k2, k2_dot)
            trial = point + 0.5_dp*step_size*k2
            trial_dot = point_dot + 0.5_dp*step_size*k2_dot
            call rhs(phi + 0.5_dp*step_size, trial, trial_dot, k3, k3_dot)
            trial = point + step_size*k3
            trial_dot = point_dot + step_size*k3_dot
            call rhs(phi + step_size, trial, trial_dot, k4, k4_dot)
            point = point + step_size*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
            point_dot = point_dot + step_size*(k1_dot + 2.0_dp*k2_dot + &
                2.0_dp*k3_dot + k4_dot)/6.0_dp
            phi = phi + step_size
        end do
        endpoint = point
        endpoint_dot = point_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine trace_fci_field_line_rk4_jvp

end module fortfem_fci_field_line_tracer
