module fortfem_surface_delta_load
    !! Explicit trace-basis assembly for a surface delta source.
    !!
    !! For quadrature traces T(q,i), weights w(q), and surface source g(q),
    !! the weak load is l_i = sum_q T(q,i) w(q) g(q).  Keeping this term
    !! separate from volume assembly makes a fitted current sheet or
    !! distributional source explicit instead of approximating a delta in a
    !! volume cell.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_surface_delta_load

contains

    subroutine assemble_surface_delta_load( &
            trace_basis, surface_weights, surface_source, load, status)
        !! Assemble a scalar surface-delta weak load from trace values.
        real(dp), intent(in) :: trace_basis(:, :)
        real(dp), intent(in) :: surface_weights(:), surface_source(:)
        real(dp), intent(out) :: load(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, dof_count, quadrature, dof
        real(dp) :: weighted_source

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface delta load received incompatible arrays")
        load = 0.0_dp
        quadrature_count = size(trace_basis, 1)
        dof_count = size(trace_basis, 2)
        if (quadrature_count < 1 .or. dof_count < 1) return
        if (size(surface_weights) /= quadrature_count .or. &
            size(surface_source) /= quadrature_count .or. &
            size(load) /= dof_count) return
        if (any(.not. ieee_is_finite(trace_basis)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(surface_source))) return
        if (any(surface_weights <= 0.0_dp)) return

        do quadrature = 1, quadrature_count
            weighted_source = surface_weights(quadrature)* &
                surface_source(quadrature)
            do dof = 1, dof_count
                load(dof) = load(dof) + trace_basis(quadrature, dof)* &
                    weighted_source
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_delta_load

end module fortfem_surface_delta_load
