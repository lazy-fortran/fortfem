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
    public :: assemble_surface_delta_load_jvp
    public :: assemble_surface_delta_load_vjp
    public :: assemble_surface_vector_delta_load
    public :: assemble_surface_vector_delta_load_jvp
    public :: assemble_surface_vector_delta_load_vjp

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

        load = 0.0_dp
        call validate_scalar_delta_inputs( &
            trace_basis, surface_weights, surface_source, &
            quadrature_count, dof_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(load) /= dof_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface delta load output has incompatible size")
            return
        end if

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

    subroutine assemble_surface_delta_load_jvp( &
            trace_basis, surface_weights, surface_source, trace_basis_dot, &
            surface_weights_dot, surface_source_dot, load_dot, status)
        !! Apply the product-rule directional derivative of a scalar load.
        real(dp), intent(in) :: trace_basis(:, :)
        real(dp), intent(in) :: surface_weights(:), surface_source(:)
        real(dp), intent(in) :: trace_basis_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:), surface_source_dot(:)
        real(dp), intent(out) :: load_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, dof_count, quadrature, dof

        load_dot = 0.0_dp
        call validate_scalar_delta_inputs( &
            trace_basis, surface_weights, surface_source, &
            quadrature_count, dof_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(load_dot) /= dof_count .or. &
            size(trace_basis_dot, 1) /= quadrature_count .or. &
            size(trace_basis_dot, 2) /= dof_count .or. &
            size(surface_weights_dot) /= quadrature_count .or. &
            size(surface_source_dot) /= quadrature_count .or. &
            any(.not. ieee_is_finite(trace_basis_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot)) .or. &
            any(.not. ieee_is_finite(surface_source_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface delta load JVP received incompatible increments")
            return
        end if

        do quadrature = 1, quadrature_count
            do dof = 1, dof_count
                load_dot(dof) = load_dot(dof) + &
                    trace_basis_dot(quadrature, dof)*surface_weights(quadrature)* &
                    surface_source(quadrature) + &
                    trace_basis(quadrature, dof)*surface_weights_dot(quadrature)* &
                    surface_source(quadrature) + &
                    trace_basis(quadrature, dof)*surface_weights(quadrature)* &
                    surface_source_dot(quadrature)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_delta_load_jvp

    subroutine assemble_surface_delta_load_vjp( &
            trace_basis, surface_weights, surface_source, load_bar, &
            trace_basis_bar, surface_weights_bar, surface_source_bar, status)
        !! Apply the real transpose of the scalar surface-load map.
        real(dp), intent(in) :: trace_basis(:, :)
        real(dp), intent(in) :: surface_weights(:), surface_source(:)
        real(dp), intent(in) :: load_bar(:)
        real(dp), intent(out) :: trace_basis_bar(:, :)
        real(dp), intent(out) :: surface_weights_bar(:), surface_source_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, dof_count, quadrature, dof

        trace_basis_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        surface_source_bar = 0.0_dp
        call validate_scalar_delta_inputs( &
            trace_basis, surface_weights, surface_source, &
            quadrature_count, dof_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(load_bar) /= dof_count .or. &
            size(trace_basis_bar, 1) /= quadrature_count .or. &
            size(trace_basis_bar, 2) /= dof_count .or. &
            size(surface_weights_bar) /= quadrature_count .or. &
            size(surface_source_bar) /= quadrature_count .or. &
            any(.not. ieee_is_finite(load_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface delta load VJP received incompatible cotangents")
            return
        end if

        do quadrature = 1, quadrature_count
            do dof = 1, dof_count
                trace_basis_bar(quadrature, dof) = &
                    trace_basis_bar(quadrature, dof) + &
                    surface_weights(quadrature)*surface_source(quadrature)* &
                    load_bar(dof)
                surface_weights_bar(quadrature) = &
                    surface_weights_bar(quadrature) + &
                    trace_basis(quadrature, dof)*surface_source(quadrature)* &
                    load_bar(dof)
                surface_source_bar(quadrature) = &
                    surface_source_bar(quadrature) + &
                    trace_basis(quadrature, dof)*surface_weights(quadrature)* &
                    load_bar(dof)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_delta_load_vjp

    subroutine assemble_surface_vector_delta_load( &
            tangential_trace_basis, surface_weights, surface_current, load, &
            status)
        !! Assemble a tangential trace pairing with a surface current.
        real(dp), intent(in) :: tangential_trace_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:)
        real(dp), intent(in) :: surface_current(:, :)
        real(dp), intent(out) :: load(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, dof_count, quadrature, dof
        real(dp) :: weighted_current

        load = 0.0_dp
        call validate_vector_delta_inputs( &
            tangential_trace_basis, surface_weights, surface_current, &
            quadrature_count, dof_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(load) /= dof_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface vector delta load output has incompatible size")
            return
        end if

        do quadrature = 1, quadrature_count
            do dof = 1, dof_count
                weighted_current = dot_product( &
                    tangential_trace_basis(quadrature, dof, :), &
                    surface_current(quadrature, :))*surface_weights(quadrature)
                load(dof) = load(dof) + weighted_current
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_vector_delta_load

    subroutine assemble_surface_vector_delta_load_jvp( &
            tangential_trace_basis, surface_weights, surface_current, &
            tangential_trace_basis_dot, surface_weights_dot, &
            surface_current_dot, load_dot, status)
        !! Apply the product-rule directional derivative of a vector load.
        real(dp), intent(in) :: tangential_trace_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:), surface_current(:, :)
        real(dp), intent(in) :: tangential_trace_basis_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:), surface_current_dot(:, :)
        real(dp), intent(out) :: load_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, dof_count, quadrature, dof

        load_dot = 0.0_dp
        call validate_vector_delta_inputs( &
            tangential_trace_basis, surface_weights, surface_current, &
            quadrature_count, dof_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(load_dot) /= dof_count .or. &
            size(tangential_trace_basis_dot, 1) /= quadrature_count .or. &
            size(tangential_trace_basis_dot, 2) /= dof_count .or. &
            size(tangential_trace_basis_dot, 3) /= 3 .or. &
            size(surface_weights_dot) /= quadrature_count .or. &
            size(surface_current_dot, 1) /= quadrature_count .or. &
            size(surface_current_dot, 2) /= 3 .or. &
            any(.not. ieee_is_finite(tangential_trace_basis_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot)) .or. &
            any(.not. ieee_is_finite(surface_current_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface vector delta load JVP received incompatible increments")
            return
        end if

        do quadrature = 1, quadrature_count
            do dof = 1, dof_count
                load_dot(dof) = load_dot(dof) + &
                    surface_weights_dot(quadrature)*dot_product( &
                    tangential_trace_basis(quadrature, dof, :), &
                    surface_current(quadrature, :)) + &
                    surface_weights(quadrature)*dot_product( &
                    tangential_trace_basis_dot(quadrature, dof, :), &
                    surface_current(quadrature, :)) + &
                    surface_weights(quadrature)*dot_product( &
                    tangential_trace_basis(quadrature, dof, :), &
                    surface_current_dot(quadrature, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_vector_delta_load_jvp

    subroutine assemble_surface_vector_delta_load_vjp( &
            tangential_trace_basis, surface_weights, surface_current, load_bar, &
            tangential_trace_basis_bar, surface_weights_bar, &
            surface_current_bar, status)
        !! Apply the real transpose of the vector surface-load map.
        real(dp), intent(in) :: tangential_trace_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:), surface_current(:, :)
        real(dp), intent(in) :: load_bar(:)
        real(dp), intent(out) :: tangential_trace_basis_bar(:, :, :)
        real(dp), intent(out) :: surface_weights_bar(:), surface_current_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, dof_count, quadrature, dof

        tangential_trace_basis_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        surface_current_bar = 0.0_dp
        call validate_vector_delta_inputs( &
            tangential_trace_basis, surface_weights, surface_current, &
            quadrature_count, dof_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(load_bar) /= dof_count .or. &
            size(tangential_trace_basis_bar, 1) /= quadrature_count .or. &
            size(tangential_trace_basis_bar, 2) /= dof_count .or. &
            size(tangential_trace_basis_bar, 3) /= 3 .or. &
            size(surface_weights_bar) /= quadrature_count .or. &
            size(surface_current_bar, 1) /= quadrature_count .or. &
            size(surface_current_bar, 2) /= 3 .or. &
            any(.not. ieee_is_finite(load_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface vector delta load VJP received incompatible cotangents")
            return
        end if

        do quadrature = 1, quadrature_count
            do dof = 1, dof_count
                tangential_trace_basis_bar(quadrature, dof, :) = &
                    tangential_trace_basis_bar(quadrature, dof, :) + &
                    surface_weights(quadrature)*surface_current(quadrature, :)* &
                    load_bar(dof)
                surface_current_bar(quadrature, :) = &
                    surface_current_bar(quadrature, :) + &
                    surface_weights(quadrature)* &
                    tangential_trace_basis(quadrature, dof, :)*load_bar(dof)
                surface_weights_bar(quadrature) = &
                    surface_weights_bar(quadrature) + load_bar(dof)*dot_product( &
                    tangential_trace_basis(quadrature, dof, :), &
                    surface_current(quadrature, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_vector_delta_load_vjp

    subroutine validate_scalar_delta_inputs( &
            trace_basis, surface_weights, surface_source, &
            quadrature_count, dof_count, status)
        real(dp), intent(in) :: trace_basis(:, :), surface_weights(:)
        real(dp), intent(in) :: surface_source(:)
        integer, intent(out) :: quadrature_count, dof_count
        type(fortsparse_status_t), intent(out) :: status

        quadrature_count = size(trace_basis, 1)
        dof_count = size(trace_basis, 2)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface delta load received incompatible arrays")
        if (quadrature_count < 1 .or. dof_count < 1) return
        if (size(surface_weights) /= quadrature_count .or. &
            size(surface_source) /= quadrature_count) return
        if (any(.not. ieee_is_finite(trace_basis)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(surface_source)) .or. &
            any(surface_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_scalar_delta_inputs

    subroutine validate_vector_delta_inputs( &
            tangential_trace_basis, surface_weights, surface_current, &
            quadrature_count, dof_count, status)
        real(dp), intent(in) :: tangential_trace_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:), surface_current(:, :)
        integer, intent(out) :: quadrature_count, dof_count
        type(fortsparse_status_t), intent(out) :: status

        quadrature_count = size(tangential_trace_basis, 1)
        dof_count = size(tangential_trace_basis, 2)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface vector delta load received incompatible arrays")
        if (quadrature_count < 1 .or. dof_count < 1) return
        if (size(tangential_trace_basis, 3) /= 3 .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(surface_current, 1) /= quadrature_count .or. &
            size(surface_current, 2) /= 3) return
        if (any(.not. ieee_is_finite(tangential_trace_basis)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(surface_current)) .or. &
            any(surface_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vector_delta_inputs

end module fortfem_surface_delta_load
