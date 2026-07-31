module fortfem_fci_parallel_operator
    !! Sparse support-operator building blocks for field-coordinate-independent
    !! (FCI) parallel derivatives.
    !!
    !! The map arrays describe interpolation from the lower and upper
    !! poloidal planes of each toroidal segment to its staggered flux boxes.
    !! The gradient is assembled as
    !!
    !!   Q = (Q_plus - Q_minus) / ell,
    !!
    !! while the support divergence is the negative volume-weighted adjoint
    !!
    !!   P = -W_c^{-1} Q^T W_s.
    !!
    !! This is the conservative pairing used by PARALLAX-style FCI support
    !! operators.  Field-line tracing and interpolation coefficient
    !! construction remain separate geometry services; this module only owns
    !! the algebraic contract and its sparse representation.
    use fortfem_kinds, only: dp
    use fortfem_generated_fci_parallel_gradient, only: &
        generated_fci_parallel_gradient
    use fortfem_generated_fci_parallel_gradient_jvp, only: &
        generated_fci_parallel_gradient_jvp
    use fortfem_generated_fci_parallel_gradient_vjp, only: &
        generated_fci_parallel_gradient_vjp
    use fortsparse, only: csc_from_triplet, csc_is_valid, csc_t, &
        fortsparse_status_t, status_set, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_OK
    implicit none
    private

    public :: assemble_fci_parallel_gradient_csc
    public :: assemble_fci_parallel_support_divergence_csc
    public :: apply_fci_parallel_gradient
    public :: apply_fci_parallel_diffusion
    public :: apply_fci_parallel_diffusion_field_vjp
    public :: apply_fci_parallel_gradient_jvp
    public :: apply_fci_parallel_gradient_vjp

contains

    subroutine assemble_fci_parallel_gradient_csc( &
            forward_map, backward_map, line_lengths, gradient, status)
        !! Assemble Q from mapped upper/lower plane interpolation matrices.
        !!
        !! `forward_map(:, :, segment)` maps the upper plane of a segment and
        !! `backward_map(:, :, segment)` maps its lower plane.  The canonical
        !! unknown vector is ordered plane-by-plane, and the staggered result
        !! is ordered segment-by-segment.  Each map may be sparse in practice;
        !! zero entries are omitted from the resulting CSC matrix.
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        type(csc_t), intent(out) :: gradient
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered
        integer :: segment, sample, plane_node, entry_count
        integer :: lower_column, upper_column, row_count, column_count
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        real(dp) :: coefficient

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI gradient assembly received incompatible maps")
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        if (n_staggered < 1 .or. n_plane < 1 .or. n_segment < 1) return
        if (any(shape(backward_map) /= shape(forward_map))) return
        if (size(line_lengths, 1) /= n_staggered) return
        if (size(line_lengths, 2) /= n_segment) return
        if (any(line_lengths <= 0.0_dp)) return

        row_count = n_staggered*n_segment
        column_count = n_plane*(n_segment + 1)
        allocate( &
            rows(2*n_staggered*n_plane*n_segment), &
            columns(2*n_staggered*n_plane*n_segment), &
            values(2*n_staggered*n_plane*n_segment))
        entry_count = 0

        do segment = 1, n_segment
            do sample = 1, n_staggered
                lower_column = (segment - 1)*n_plane
                upper_column = segment*n_plane
                do plane_node = 1, n_plane
                    coefficient = -backward_map(sample, plane_node, segment) / &
                        line_lengths(sample, segment)
                    if (coefficient /= 0.0_dp) then
                        entry_count = entry_count + 1
                        rows(entry_count) = &
                            sample + (segment - 1)*n_staggered
                        columns(entry_count) = lower_column + plane_node
                        values(entry_count) = coefficient
                    end if
                    coefficient = forward_map(sample, plane_node, segment) / &
                        line_lengths(sample, segment)
                    if (coefficient /= 0.0_dp) then
                        entry_count = entry_count + 1
                        rows(entry_count) = &
                            sample + (segment - 1)*n_staggered
                        columns(entry_count) = upper_column + plane_node
                        values(entry_count) = coefficient
                    end if
                end do
            end do
        end do

        call csc_from_triplet( &
            row_count, column_count, rows(:entry_count), columns(:entry_count), &
            values(:entry_count), gradient, status)
        if (status%code == FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI gradient CSC construction failed")
    end subroutine assemble_fci_parallel_gradient_csc

    subroutine assemble_fci_parallel_support_divergence_csc( &
            gradient, canonical_volumes, staggered_volumes, divergence, status)
        !! Assemble P = -W_c^{-1} Q^T W_s from a FCI gradient Q.
        !!
        !! The returned operator satisfies, up to floating-point roundoff,
        !!
        !!   u^T W_c P f = -(Q u)^T W_s f,
        !!
        !! which is the support-operator conservation/adjointness contract.
        type(csc_t), intent(in) :: gradient
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: staggered_volumes(:)
        type(csc_t), intent(out) :: divergence
        type(fortsparse_status_t), intent(out) :: status

        integer :: column, entry, row
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI support divergence received an invalid gradient")
        if (.not. csc_is_valid(gradient)) return
        if (size(canonical_volumes) /= gradient%ncol) return
        if (size(staggered_volumes) /= gradient%nrow) return
        if (any(canonical_volumes <= 0.0_dp)) return
        if (any(staggered_volumes <= 0.0_dp)) return

        allocate(rows(gradient%nnz), columns(gradient%nnz))
        allocate(values(gradient%nnz))
        entry = 0
        do column = 1, gradient%ncol
            do row = gradient%col_ptr(column), gradient%col_ptr(column + 1) - 1
                entry = entry + 1
                rows(entry) = column
                columns(entry) = gradient%row_idx(row)
                values(entry) = -staggered_volumes(gradient%row_idx(row)) / &
                    canonical_volumes(column)*gradient%val(row)
            end do
        end do

        call csc_from_triplet( &
            gradient%ncol, gradient%nrow, rows, columns, values, divergence, &
            status)
        if (status%code == FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI support divergence CSC construction failed")
    end subroutine assemble_fci_parallel_support_divergence_csc

    subroutine apply_fci_parallel_gradient( &
            forward_map, backward_map, line_lengths, field, gradient_field, &
            status)
        !! Apply the mapped FCI gradient without assembling a global matrix.
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        real(dp), intent(in) :: field(:)
        real(dp), intent(out) :: gradient_field(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered
        integer :: segment, sample, plane_node, row
        integer :: lower_column, upper_column
        real(dp) :: contribution

        call validate_fci_action_shapes( &
            forward_map, backward_map, line_lengths, field, gradient_field, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        gradient_field = 0.0_dp
        do segment = 1, n_segment
            lower_column = (segment - 1)*n_plane
            upper_column = segment*n_plane
            do sample = 1, n_staggered
                row = sample + (segment - 1)*n_staggered
                do plane_node = 1, n_plane
                    call generated_fci_parallel_gradient( &
                        forward_map(sample, plane_node, segment), &
                        field(upper_column + plane_node), &
                        backward_map(sample, plane_node, segment), &
                        field(lower_column + plane_node), &
                        line_lengths(sample, segment), contribution)
                    gradient_field(row) = gradient_field(row) + contribution
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fci_parallel_gradient

    subroutine apply_fci_parallel_diffusion( &
            forward_map, backward_map, line_lengths, parallel_coefficient, &
            canonical_volumes, staggered_volumes, field, diffusion_field, &
            status)
        !! Apply the matrix-free FCI parallel diffusion action
        !!
        !!   L_parallel = -W_c^{-1} Q^T W_s K_parallel Q.
        !!
        !! Positive coefficients and volumes therefore give the weighted
        !! negative-energy identity
        !!
        !!   u^T W_c L_parallel u = -(Q u)^T W_s K_parallel Q u.
        !!
        !! This is the support-operator ingredient; perpendicular terms,
        !! boundary conditions, and preconditioning remain higher-level
        !! assembly concerns.
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        real(dp), intent(in) :: parallel_coefficient(:)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: staggered_volumes(:)
        real(dp), intent(in) :: field(:)
        real(dp), intent(out) :: diffusion_field(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered
        integer :: segment, sample, plane_node, row
        integer :: lower_column, upper_column
        real(dp), allocatable :: gradient_field(:)
        real(dp) :: weighted_flux, coefficient

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI parallel diffusion received incompatible arrays")
        diffusion_field = 0.0_dp
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        if (n_staggered < 1 .or. n_plane < 1 .or. n_segment < 1) return
        if (any(shape(backward_map) /= shape(forward_map))) return
        if (size(line_lengths, 1) /= n_staggered) return
        if (size(line_lengths, 2) /= n_segment) return
        if (size(parallel_coefficient) /= n_staggered*n_segment) return
        if (size(staggered_volumes) /= n_staggered*n_segment) return
        if (size(canonical_volumes) /= n_plane*(n_segment + 1)) return
        if (size(field) /= n_plane*(n_segment + 1)) return
        if (size(diffusion_field) /= n_plane*(n_segment + 1)) return
        if (any(line_lengths <= 0.0_dp)) return
        if (any(parallel_coefficient <= 0.0_dp)) return
        if (any(canonical_volumes <= 0.0_dp)) return
        if (any(staggered_volumes <= 0.0_dp)) return

        allocate(gradient_field(n_staggered*n_segment))
        call apply_fci_parallel_gradient( &
            forward_map, backward_map, line_lengths, field, gradient_field, &
            status)
        if (status%code /= FORTSPARSE_OK) return

        do segment = 1, n_segment
            lower_column = (segment - 1)*n_plane
            upper_column = segment*n_plane
            do sample = 1, n_staggered
                row = sample + (segment - 1)*n_staggered
                weighted_flux = staggered_volumes(row)* &
                    parallel_coefficient(row)*gradient_field(row)
                do plane_node = 1, n_plane
                    coefficient = backward_map(sample, plane_node, segment) / &
                        line_lengths(sample, segment)
                    diffusion_field(lower_column + plane_node) = &
                        diffusion_field(lower_column + plane_node) + &
                        coefficient*weighted_flux/ &
                        canonical_volumes(lower_column + plane_node)
                    coefficient = forward_map(sample, plane_node, segment) / &
                        line_lengths(sample, segment)
                    diffusion_field(upper_column + plane_node) = &
                        diffusion_field(upper_column + plane_node) - &
                        coefficient*weighted_flux/ &
                        canonical_volumes(upper_column + plane_node)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fci_parallel_diffusion

    subroutine apply_fci_parallel_diffusion_field_vjp( &
            forward_map, backward_map, line_lengths, parallel_coefficient, &
            canonical_volumes, staggered_volumes, diffusion_field_bar, &
            field_bar, status)
        !! Apply the VJP with respect to the input field of the FCI action.
        !!
        !! For fixed maps, coefficients, and volumes, the action is
        !!
        !!   L = -W_c^{-1} Q^T W_s K_parallel Q,
        !!
        !! so its Euclidean transpose is evaluated by first applying Q to
        !! W_c^{-1} times the output cotangent and then using the generated
        !! gradient VJP.  Keeping this composition explicit preserves the
        !! same FortSym-generated scalar product as the primal gradient.
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        real(dp), intent(in) :: parallel_coefficient(:)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: staggered_volumes(:)
        real(dp), intent(in) :: diffusion_field_bar(:)
        real(dp), intent(out) :: field_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered
        real(dp), allocatable :: scaled_bar(:), gradient_bar(:)
        real(dp), allocatable :: zero_field(:)
        real(dp), allocatable :: forward_map_bar(:, :, :)
        real(dp), allocatable :: backward_map_bar(:, :, :)
        real(dp), allocatable :: line_lengths_bar(:, :)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI diffusion VJP received incompatible arrays")
        field_bar = 0.0_dp
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        if (n_staggered < 1 .or. n_plane < 1 .or. n_segment < 1) return
        if (any(shape(backward_map) /= shape(forward_map))) return
        if (size(line_lengths, 1) /= n_staggered) return
        if (size(line_lengths, 2) /= n_segment) return
        if (size(parallel_coefficient) /= n_staggered*n_segment) return
        if (size(staggered_volumes) /= n_staggered*n_segment) return
        if (size(canonical_volumes) /= n_plane*(n_segment + 1)) return
        if (size(diffusion_field_bar) /= n_plane*(n_segment + 1)) return
        if (size(field_bar) /= n_plane*(n_segment + 1)) return
        if (any(line_lengths <= 0.0_dp)) return
        if (any(parallel_coefficient <= 0.0_dp)) return
        if (any(canonical_volumes <= 0.0_dp)) return
        if (any(staggered_volumes <= 0.0_dp)) return

        allocate(scaled_bar(size(diffusion_field_bar)))
        allocate(gradient_bar(n_staggered*n_segment))
        allocate(zero_field(size(diffusion_field_bar)))
        allocate(forward_map_bar(n_staggered, n_plane, n_segment))
        allocate(backward_map_bar(n_staggered, n_plane, n_segment))
        allocate(line_lengths_bar(n_staggered, n_segment))
        scaled_bar = diffusion_field_bar/canonical_volumes
        call apply_fci_parallel_gradient( &
            forward_map, backward_map, line_lengths, scaled_bar, gradient_bar, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        gradient_bar = -staggered_volumes*parallel_coefficient*gradient_bar
        zero_field = 0.0_dp
        call apply_fci_parallel_gradient_vjp( &
            forward_map, backward_map, line_lengths, zero_field, gradient_bar, &
            forward_map_bar, backward_map_bar, line_lengths_bar, field_bar, &
            status)
    end subroutine apply_fci_parallel_diffusion_field_vjp

    subroutine apply_fci_parallel_gradient_jvp( &
            forward_map, backward_map, line_lengths, field, forward_map_dot, &
            backward_map_dot, line_lengths_dot, field_dot, gradient_field_dot, &
            status)
        !! Apply the analytical JVP of the mapped FCI gradient action.
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        real(dp), intent(in) :: field(:)
        real(dp), intent(in) :: forward_map_dot(:, :, :)
        real(dp), intent(in) :: backward_map_dot(:, :, :)
        real(dp), intent(in) :: line_lengths_dot(:, :)
        real(dp), intent(in) :: field_dot(:)
        real(dp), intent(out) :: gradient_field_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered
        integer :: segment, sample, plane_node, row
        integer :: lower_column, upper_column
        real(dp) :: contribution_dot

        call validate_fci_action_shapes( &
            forward_map, backward_map, line_lengths, field, gradient_field_dot, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(forward_map_dot) /= shape(forward_map))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient JVP received an incompatible forward map")
            return
        end if
        if (any(shape(backward_map_dot) /= shape(backward_map))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient JVP received an incompatible backward map")
            return
        end if
        if (any(shape(line_lengths_dot) /= shape(line_lengths))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient JVP received incompatible line lengths")
            return
        end if
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        if (size(field_dot) /= n_plane*(n_segment + 1)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient JVP received an incompatible field tangent")
            return
        end if
        gradient_field_dot = 0.0_dp
        do segment = 1, n_segment
            lower_column = (segment - 1)*n_plane
            upper_column = segment*n_plane
            do sample = 1, n_staggered
                row = sample + (segment - 1)*n_staggered
                do plane_node = 1, n_plane
                    call generated_fci_parallel_gradient_jvp( &
                        forward_map(sample, plane_node, segment), &
                        field(upper_column + plane_node), &
                        backward_map(sample, plane_node, segment), &
                        field(lower_column + plane_node), &
                        line_lengths(sample, segment), &
                        forward_map_dot(sample, plane_node, segment), &
                        field_dot(upper_column + plane_node), &
                        backward_map_dot(sample, plane_node, segment), &
                        field_dot(lower_column + plane_node), &
                        line_lengths_dot(sample, segment), contribution_dot)
                    gradient_field_dot(row) = &
                        gradient_field_dot(row) + contribution_dot
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fci_parallel_gradient_jvp

    subroutine apply_fci_parallel_gradient_vjp( &
            forward_map, backward_map, line_lengths, field, gradient_field_bar, &
            forward_map_bar, backward_map_bar, line_lengths_bar, field_bar, &
            status)
        !! Apply the real VJP of the mapped FCI gradient action.
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        real(dp), intent(in) :: field(:)
        real(dp), intent(in) :: gradient_field_bar(:)
        real(dp), intent(out) :: forward_map_bar(:, :, :)
        real(dp), intent(out) :: backward_map_bar(:, :, :)
        real(dp), intent(out) :: line_lengths_bar(:, :)
        real(dp), intent(out) :: field_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered
        integer :: segment, sample, plane_node, row
        integer :: lower_column, upper_column
        real(dp) :: forward_bar, upper_bar, backward_bar, lower_bar
        real(dp) :: length_bar

        call validate_fci_action_shapes( &
            forward_map, backward_map, line_lengths, field, gradient_field_bar, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(forward_map_bar) /= shape(forward_map))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient VJP received an incompatible forward cotangent")
            return
        end if
        if (any(shape(backward_map_bar) /= shape(backward_map))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient VJP received an incompatible backward cotangent")
            return
        end if
        if (any(shape(line_lengths_bar) /= shape(line_lengths))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient VJP received incompatible line-length cotangents")
            return
        end if
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        if (size(field_bar) /= n_plane*(n_segment + 1)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI gradient VJP received an incompatible field cotangent")
            return
        end if
        forward_map_bar = 0.0_dp
        backward_map_bar = 0.0_dp
        line_lengths_bar = 0.0_dp
        field_bar = 0.0_dp
        do segment = 1, n_segment
            lower_column = (segment - 1)*n_plane
            upper_column = segment*n_plane
            do sample = 1, n_staggered
                row = sample + (segment - 1)*n_staggered
                do plane_node = 1, n_plane
                    call generated_fci_parallel_gradient_vjp( &
                        forward_map(sample, plane_node, segment), &
                        field(upper_column + plane_node), &
                        backward_map(sample, plane_node, segment), &
                        field(lower_column + plane_node), &
                        line_lengths(sample, segment), gradient_field_bar(row), &
                        forward_bar, upper_bar, backward_bar, lower_bar, &
                        length_bar)
                    forward_map_bar(sample, plane_node, segment) = &
                        forward_map_bar(sample, plane_node, segment) + forward_bar
                    backward_map_bar(sample, plane_node, segment) = &
                        backward_map_bar(sample, plane_node, segment) + backward_bar
                    line_lengths_bar(sample, segment) = &
                        line_lengths_bar(sample, segment) + length_bar
                    field_bar(upper_column + plane_node) = &
                        field_bar(upper_column + plane_node) + upper_bar
                    field_bar(lower_column + plane_node) = &
                        field_bar(lower_column + plane_node) + lower_bar
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fci_parallel_gradient_vjp

    subroutine validate_fci_action_shapes( &
            forward_map, backward_map, line_lengths, field, output, status)
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        real(dp), intent(in) :: line_lengths(:, :)
        real(dp), intent(in) :: field(:)
        real(dp), intent(in) :: output(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_plane, n_segment, n_staggered

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI gradient action received incompatible arrays")
        n_staggered = size(forward_map, 1)
        n_plane = size(forward_map, 2)
        n_segment = size(forward_map, 3)
        if (n_staggered < 1 .or. n_plane < 1 .or. n_segment < 1) return
        if (any(shape(backward_map) /= shape(forward_map))) return
        if (size(line_lengths, 1) /= n_staggered) return
        if (size(line_lengths, 2) /= n_segment) return
        if (any(line_lengths <= 0.0_dp)) return
        if (size(field) /= n_plane*(n_segment + 1)) return
        if (size(output) /= n_staggered*n_segment) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_fci_action_shapes

end module fortfem_fci_parallel_operator
