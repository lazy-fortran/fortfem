module fortfem_surface_vector_trace
    !! Geometry-neutral normal/tangential decomposition on a surface.
    !!
    !! FEM, BEM, DtN, PML, and visualization clients often receive a
    !! three-dimensional field sampled at common surface points.  This
    !! contract turns those caller-owned samples and normals into the scalar
    !! normal trace and the vector tangential trace without selecting an
    !! equation, basis, or boundary convention.  Non-unit normals are
    !! normalized locally, so meshes and exact-curved surfaces share the same
    !! value/JVP/VJP behavior.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_surface_vector_trace
    public :: evaluate_surface_vector_trace_jvp
    public :: evaluate_surface_vector_trace_vjp

contains

    subroutine evaluate_surface_vector_trace( &
            field, normals, normal_component, tangential_field, status)
        !! Decompose each supplied field sample into normal and tangent parts.
        real(dp), intent(in) :: field(:, :), normals(:, :)
        real(dp), intent(out) :: normal_component(:), tangential_field(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: normal_length, unit_normal(3)

        normal_component = 0.0_dp
        tangential_field = 0.0_dp
        if (.not. valid_base_inputs(field, normals, normal_component, &
                tangential_field)) then
            call invalid_status(status, &
                "surface vector trace has incompatible or non-finite samples")
            return
        end if
        do sample = 1, size(field, 2)
            normal_length = sqrt(dot_product(normals(:, sample), normals(:, sample)))
            unit_normal = normals(:, sample)/normal_length
            normal_component(sample) = dot_product(unit_normal, field(:, sample))
            tangential_field(:, sample) = field(:, sample) - &
                unit_normal*normal_component(sample)
        end do
        if (.not. all(ieee_is_finite(normal_component)) .or. &
            .not. all(ieee_is_finite(tangential_field))) then
            normal_component = 0.0_dp
            tangential_field = 0.0_dp
            call invalid_status(status, "surface vector trace is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_vector_trace

    subroutine evaluate_surface_vector_trace_jvp( &
            field, normals, field_dot, normals_dot, normal_component_dot, &
            tangential_field_dot, status)
        !! Fixed-topology directional derivative of the trace decomposition.
        real(dp), intent(in) :: field(:, :), normals(:, :)
        real(dp), intent(in) :: field_dot(:, :), normals_dot(:, :)
        real(dp), intent(out) :: normal_component_dot(:), tangential_field_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: normal_length, unit_normal(3), unit_normal_dot(3)
        real(dp) :: normal_component, normal_component_dot_local

        normal_component_dot = 0.0_dp
        tangential_field_dot = 0.0_dp
        if (.not. valid_base_inputs(field, normals, normal_component_dot, &
                tangential_field_dot) .or. .not. valid_direction_inputs( &
                field, normals, field_dot, normals_dot)) then
            call invalid_status(status, &
                "surface vector trace JVP has incompatible samples")
            return
        end if
        do sample = 1, size(field, 2)
            normal_length = sqrt(dot_product(normals(:, sample), normals(:, sample)))
            unit_normal = normals(:, sample)/normal_length
            unit_normal_dot = (normals_dot(:, sample) - unit_normal* &
                dot_product(unit_normal, normals_dot(:, sample)))/normal_length
            normal_component = dot_product(unit_normal, field(:, sample))
            normal_component_dot_local = dot_product(unit_normal_dot, &
                field(:, sample)) + &
                dot_product(unit_normal, field_dot(:, sample))
            normal_component_dot(sample) = normal_component_dot_local
            tangential_field_dot(:, sample) = field_dot(:, sample) - &
                unit_normal_dot*normal_component - &
                unit_normal*normal_component_dot_local
        end do
        if (.not. all(ieee_is_finite(normal_component_dot)) .or. &
            .not. all(ieee_is_finite(tangential_field_dot))) then
            normal_component_dot = 0.0_dp
            tangential_field_dot = 0.0_dp
            call invalid_status(status, "surface vector trace JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_vector_trace_jvp

    subroutine evaluate_surface_vector_trace_vjp( &
            field, normals, normal_component_bar, tangential_field_bar, &
            field_bar, normals_bar, status)
        !! Real Euclidean adjoint of the normal/tangential trace map.
        real(dp), intent(in) :: field(:, :), normals(:, :)
        real(dp), intent(in) :: normal_component_bar(:), tangential_field_bar(:, :)
        real(dp), intent(out) :: field_bar(:, :), normals_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: normal_length, unit_normal(3), projector(3, 3)
        real(dp) :: normal_component, unit_normal_bar(3), tangent_bar(3)

        field_bar = 0.0_dp
        normals_bar = 0.0_dp
        if (.not. valid_base_inputs(field, normals, normal_component_bar, &
                tangential_field_bar) .or. size(field_bar, 1) /= 3 .or. &
            size(field_bar, 2) /= size(field, 2) .or. &
            size(normals_bar, 1) /= 3 .or. size(normals_bar, 2) /= size(field, 2) .or. &
            .not. all(ieee_is_finite(normal_component_bar)) .or. &
            .not. all(ieee_is_finite(tangential_field_bar))) then
            call invalid_status(status, "surface vector trace VJP has invalid samples")
            return
        end if
        do sample = 1, size(field, 2)
            normal_length = sqrt(dot_product(normals(:, sample), normals(:, sample)))
            unit_normal = normals(:, sample)/normal_length
            projector = 0.0_dp
            projector(1, 1) = 1.0_dp
            projector(2, 2) = 1.0_dp
            projector(3, 3) = 1.0_dp
            projector = projector - outer_product(unit_normal, unit_normal)
            normal_component = dot_product(unit_normal, field(:, sample))
            tangent_bar = tangential_field_bar(:, sample)
            field_bar(:, sample) = normal_component_bar(sample)*unit_normal + &
                matmul(projector, tangent_bar)
            unit_normal_bar = normal_component_bar(sample)*field(:, sample) - &
                normal_component*tangent_bar - &
                dot_product(tangent_bar, unit_normal)*field(:, sample)
            normals_bar(:, sample) = matmul(projector, unit_normal_bar)/normal_length
        end do
        if (.not. all(ieee_is_finite(field_bar)) .or. &
            .not. all(ieee_is_finite(normals_bar))) then
            field_bar = 0.0_dp
            normals_bar = 0.0_dp
            call invalid_status(status, "surface vector trace VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_vector_trace_vjp

    logical function valid_base_inputs( &
            field, normals, scalar_output, vector_output) result(valid)
        real(dp), intent(in) :: field(:, :), normals(:, :)
        real(dp), intent(in) :: scalar_output(:), vector_output(:, :)
        real(dp) :: normal_length
        integer :: sample

        valid = size(field, 1) == 3 .and. size(field, 2) > 0 .and. &
            size(normals, 1) == 3 .and. size(normals, 2) == size(field, 2) .and. &
            size(scalar_output) == size(field, 2) .and. &
            size(vector_output, 1) == 3 .and. &
            size(vector_output, 2) == size(field, 2) .and. &
            all(ieee_is_finite(field)) .and. all(ieee_is_finite(normals))
        if (.not. valid) return
        do sample = 1, size(field, 2)
            normal_length = dot_product(normals(:, sample), normals(:, sample))
            if (.not. ieee_is_finite(normal_length) .or. &
                normal_length <= tiny(1.0_dp)) then
                valid = .false.
                return
            end if
        end do
    end function valid_base_inputs

    logical function valid_direction_inputs( &
            field, normals, field_dot, normals_dot) result(valid)
        real(dp), intent(in) :: field(:, :), normals(:, :)
        real(dp), intent(in) :: field_dot(:, :), normals_dot(:, :)

        valid = size(field_dot, 1) == size(field, 1) .and. &
            size(field_dot, 2) == size(field, 2) .and. &
            size(normals_dot, 1) == size(normals, 1) .and. &
            size(normals_dot, 2) == size(normals, 2) .and. &
            all(ieee_is_finite(field_dot)) .and. all(ieee_is_finite(normals_dot))
    end function valid_direction_inputs

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3, 3)
        integer :: row, column

        do row = 1, 3
            do column = 1, 3
                product(row, column) = left(row)*right(column)
            end do
        end do
    end function outer_product

    subroutine invalid_status(status, message)
        type(fortsparse_status_t), intent(out) :: status
        character(len=*), intent(in) :: message

        call status_set(status, FORTSPARSE_INVALID_MATRIX, message)
    end subroutine invalid_status

end module fortfem_surface_vector_trace
