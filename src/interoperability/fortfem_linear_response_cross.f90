module fortfem_linear_response_cross
    !! Neutral cross-composition oracle for modal linear-response clients.
    !!
    !! The fixed-topology residual is the concatenation
    !!
    !!   r_e = K u - lambda M u,
    !!   r_r = A y - G u - f,
    !!   r_s = H y - h.
    !!
    !! Here the first block is a generalized eigen residual, the second is a
    !! retained complex response block, and the third is an explicitly
    !! supplied shielding/trace block.  The metadata records whether the
    !! latter is interpreted as an ideal-shielding trace; it does not select
    !! a plasma closure or a response-matrix convention.  All matrices,
    !! source terms, and topology are caller-owned.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    type, public :: linear_response_cross_metadata_t
        character(len=32) :: schema_version = "fortfem-linear-cross-1"
        character(len=64) :: normalization_label = ""
        character(len=128) :: provenance = ""
        real(dp) :: normalization_scale = 1.0_dp
        logical :: ideal_shielding = .false.
        integer :: retained_response_count = 0
        integer :: shielding_trace_count = 0
        integer, allocatable :: retained_indices(:)
        character(len=64), allocatable :: response_labels(:)
    contains
        procedure, private :: assign_linear_response_cross_metadata
        generic :: assignment(=) => assign_linear_response_cross_metadata
    end type linear_response_cross_metadata_t

    public :: initialize_linear_response_cross_metadata
    public :: validate_linear_response_cross_metadata
    public :: assemble_linear_response_eigen_cross_residual
    public :: assemble_linear_response_eigen_cross_residual_jvp
    public :: assemble_linear_response_eigen_cross_residual_vjp

    interface finite_complex
        module procedure finite_complex_scalar
        module procedure finite_complex_vector
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine initialize_linear_response_cross_metadata( &
            metadata, normalization_label, normalization_scale, provenance, &
            retained_indices, response_labels, shielding_trace_count, &
            ideal_shielding, status)
        type(linear_response_cross_metadata_t), intent(out) :: metadata
        character(len=*), intent(in) :: normalization_label, provenance
        real(dp), intent(in) :: normalization_scale
        integer, intent(in) :: retained_indices(:)
        character(len=*), intent(in) :: response_labels(:)
        integer, intent(in) :: shielding_trace_count
        logical, intent(in) :: ideal_shielding
        integer, intent(out) :: status

        metadata = linear_response_cross_metadata_t()
        status = 1
        if (len_trim(normalization_label) < 1 .or. &
            len_trim(normalization_label) > len(metadata%normalization_label)) return
        if (len_trim(provenance) < 1 .or. &
            len_trim(provenance) > len(metadata%provenance)) return
        if (.not. ieee_is_finite(normalization_scale) .or. &
            normalization_scale <= 0.0_dp) return
        if (size(retained_indices) < 1 .or. &
            size(response_labels) /= size(retained_indices) .or. &
            shielding_trace_count < 0) return
        if (any(retained_indices < 1) .or. &
            .not. unique_indices(retained_indices)) return
        if (.not. unique_labels(response_labels)) return

        metadata%normalization_label = trim(normalization_label)
        metadata%normalization_scale = normalization_scale
        metadata%provenance = trim(provenance)
        metadata%ideal_shielding = ideal_shielding
        metadata%retained_response_count = size(retained_indices)
        metadata%shielding_trace_count = shielding_trace_count
        allocate(metadata%retained_indices, source=retained_indices)
        allocate(metadata%response_labels(size(response_labels)))
        metadata%response_labels = response_labels
        if (.not. validate_linear_response_cross_metadata(metadata, status)) then
            call clear_metadata(metadata)
            return
        end if
    end subroutine initialize_linear_response_cross_metadata

    logical function validate_linear_response_cross_metadata(metadata, status) &
            result(valid)
        type(linear_response_cross_metadata_t), intent(in) :: metadata
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (metadata%schema_version /= "fortfem-linear-cross-1") return
        if (len_trim(metadata%normalization_label) < 1 .or. &
            len_trim(metadata%provenance) < 1) return
        if (.not. ieee_is_finite(metadata%normalization_scale) .or. &
            metadata%normalization_scale <= 0.0_dp) return
        if (metadata%retained_response_count < 1 .or. &
            metadata%shielding_trace_count < 0) return
        if (.not. allocated(metadata%retained_indices) .or. &
            .not. allocated(metadata%response_labels)) return
        if (size(metadata%retained_indices) /= metadata%retained_response_count .or. &
            size(metadata%response_labels) /= metadata%retained_response_count) return
        if (any(metadata%retained_indices < 1) .or. &
            .not. unique_indices(metadata%retained_indices) .or. &
            .not. unique_labels(metadata%response_labels)) return
        valid = .true.
        status = 0
    end function validate_linear_response_cross_metadata

    subroutine assemble_linear_response_eigen_cross_residual( &
            stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, metadata, residual, status)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenvalue
        complex(dp), intent(in) :: eigenstate(:), response_matrix(:, :)
        complex(dp), intent(in) :: response_coupling(:, :), response_state(:)
        complex(dp), intent(in) :: response_source(:), shielding_matrix(:, :)
        complex(dp), intent(in) :: shielding_target(:)
        type(linear_response_cross_metadata_t), intent(in) :: metadata
        complex(dp), intent(out) :: residual(:)
        integer, intent(out) :: status
        integer :: eigen_count, response_count, shielding_count

        residual = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, metadata, residual)) return
        eigen_count = size(eigenstate)
        response_count = size(response_state)
        shielding_count = size(shielding_target)
        residual(1:eigen_count) = matmul(stiffness, eigenstate) - &
            eigenvalue*matmul(mass, eigenstate)
        residual(eigen_count + 1:eigen_count + response_count) = &
            matmul(response_matrix, response_state) - &
            matmul(response_coupling, eigenstate) - response_source
        if (shielding_count > 0) residual(eigen_count + response_count + 1:) = &
            matmul(shielding_matrix, response_state) - shielding_target
        status = 0
    end subroutine assemble_linear_response_eigen_cross_residual

    subroutine assemble_linear_response_eigen_cross_residual_jvp( &
            stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, metadata, stiffness_dot, mass_dot, eigenvalue_dot, &
            eigenstate_dot, response_matrix_dot, response_coupling_dot, &
            response_state_dot, response_source_dot, shielding_matrix_dot, &
            shielding_target_dot, residual_dot, status)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenvalue
        complex(dp), intent(in) :: eigenstate(:), response_matrix(:, :)
        complex(dp), intent(in) :: response_coupling(:, :), response_state(:)
        complex(dp), intent(in) :: response_source(:), shielding_matrix(:, :)
        complex(dp), intent(in) :: shielding_target(:)
        type(linear_response_cross_metadata_t), intent(in) :: metadata
        complex(dp), intent(in) :: stiffness_dot(:, :), mass_dot(:, :)
        complex(dp), intent(in) :: eigenvalue_dot, eigenstate_dot(:)
        complex(dp), intent(in) :: response_matrix_dot(:, :)
        complex(dp), intent(in) :: response_coupling_dot(:, :)
        complex(dp), intent(in) :: response_state_dot(:), response_source_dot(:)
        complex(dp), intent(in) :: shielding_matrix_dot(:, :)
        complex(dp), intent(in) :: shielding_target_dot(:)
        complex(dp), intent(out) :: residual_dot(:)
        integer, intent(out) :: status
        integer :: eigen_count, response_count, shielding_count

        residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, metadata, residual_dot)) return
        if (.not. valid_tangent_inputs( &
            stiffness, mass, eigenstate, response_matrix, response_coupling, &
            response_state, response_source, shielding_matrix, shielding_target, &
            stiffness_dot, mass_dot, eigenstate_dot, response_matrix_dot, &
            response_coupling_dot, response_state_dot, response_source_dot, &
            shielding_matrix_dot, shielding_target_dot) .or. &
            .not. finite_complex(eigenvalue_dot)) return
        eigen_count = size(eigenstate)
        response_count = size(response_state)
        shielding_count = size(shielding_target)
        residual_dot(1:eigen_count) = matmul(stiffness_dot, eigenstate) + &
            matmul(stiffness, eigenstate_dot) - &
            eigenvalue_dot*matmul(mass, eigenstate) - &
            eigenvalue*matmul(mass_dot, eigenstate) - &
            eigenvalue*matmul(mass, eigenstate_dot)
        residual_dot(eigen_count + 1:eigen_count + response_count) = &
            matmul(response_matrix_dot, response_state) + &
            matmul(response_matrix, response_state_dot) - &
            matmul(response_coupling_dot, eigenstate) - &
            matmul(response_coupling, eigenstate_dot) - response_source_dot
        if (shielding_count > 0) residual_dot(eigen_count + response_count + 1:) = &
            matmul(shielding_matrix_dot, response_state) + &
            matmul(shielding_matrix, response_state_dot) - shielding_target_dot
        status = 0
    end subroutine assemble_linear_response_eigen_cross_residual_jvp

    subroutine assemble_linear_response_eigen_cross_residual_vjp( &
            stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, metadata, residual_bar, stiffness_bar, mass_bar, &
            eigenvalue_bar, eigenstate_bar, response_matrix_bar, &
            response_coupling_bar, response_state_bar, response_source_bar, &
            shielding_matrix_bar, shielding_target_bar, status)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenvalue
        complex(dp), intent(in) :: eigenstate(:), response_matrix(:, :)
        complex(dp), intent(in) :: response_coupling(:, :), response_state(:)
        complex(dp), intent(in) :: response_source(:), shielding_matrix(:, :)
        complex(dp), intent(in) :: shielding_target(:)
        type(linear_response_cross_metadata_t), intent(in) :: metadata
        complex(dp), intent(in) :: residual_bar(:)
        complex(dp), intent(out) :: stiffness_bar(:, :), mass_bar(:, :)
        complex(dp), intent(out) :: eigenvalue_bar, eigenstate_bar(:)
        complex(dp), intent(out) :: response_matrix_bar(:, :)
        complex(dp), intent(out) :: response_coupling_bar(:, :)
        complex(dp), intent(out) :: response_state_bar(:), response_source_bar(:)
        complex(dp), intent(out) :: shielding_matrix_bar(:, :)
        complex(dp), intent(out) :: shielding_target_bar(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: mass_state(:)
        integer :: eigen_count, response_count, shielding_count

        stiffness_bar = cmplx(0.0_dp, 0.0_dp, dp)
        mass_bar = cmplx(0.0_dp, 0.0_dp, dp)
        eigenvalue_bar = cmplx(0.0_dp, 0.0_dp, dp)
        eigenstate_bar = cmplx(0.0_dp, 0.0_dp, dp)
        response_matrix_bar = cmplx(0.0_dp, 0.0_dp, dp)
        response_coupling_bar = cmplx(0.0_dp, 0.0_dp, dp)
        response_state_bar = cmplx(0.0_dp, 0.0_dp, dp)
        response_source_bar = cmplx(0.0_dp, 0.0_dp, dp)
        shielding_matrix_bar = cmplx(0.0_dp, 0.0_dp, dp)
        shielding_target_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, metadata, residual_bar) .or. &
            .not. valid_vjp_outputs( &
            stiffness, mass, eigenstate, response_matrix, response_coupling, &
            response_state, response_source, shielding_matrix, shielding_target, &
            stiffness_bar, mass_bar, eigenstate_bar, response_matrix_bar, &
            response_coupling_bar, response_state_bar, response_source_bar, &
            shielding_matrix_bar, shielding_target_bar)) return

        eigen_count = size(eigenstate)
        response_count = size(response_state)
        shielding_count = size(shielding_target)
        allocate(mass_state(eigen_count))
        mass_state = matmul(mass, eigenstate)
        call rank_one_product(stiffness_bar, residual_bar(1:eigen_count), eigenstate)
        call rank_one_product(mass_bar, residual_bar(1:eigen_count), eigenstate)
        mass_bar = -conjg(eigenvalue)*mass_bar
        eigenvalue_bar = -sum(residual_bar(1:eigen_count)*conjg(mass_state))
        eigenstate_bar = matmul(conjg(transpose(stiffness)), residual_bar(1:eigen_count)) - &
            conjg(eigenvalue)*matmul(conjg(transpose(mass)), residual_bar(1:eigen_count))

        call rank_one_product(response_matrix_bar, &
            residual_bar(eigen_count + 1:eigen_count + response_count), response_state)
        call rank_one_product(response_coupling_bar, &
            residual_bar(eigen_count + 1:eigen_count + response_count), eigenstate)
        response_coupling_bar = -response_coupling_bar
        response_state_bar = matmul(conjg(transpose(response_matrix)), &
            residual_bar(eigen_count + 1:eigen_count + response_count))
        eigenstate_bar = eigenstate_bar - matmul(conjg(transpose(response_coupling)), &
            residual_bar(eigen_count + 1:eigen_count + response_count))
        response_source_bar = -residual_bar(eigen_count + 1:eigen_count + response_count)

        if (shielding_count > 0) then
            call rank_one_product(shielding_matrix_bar, &
                residual_bar(eigen_count + response_count + 1:), response_state)
            response_state_bar = response_state_bar + &
                matmul(conjg(transpose(shielding_matrix)), &
                residual_bar(eigen_count + response_count + 1:))
            shielding_target_bar = -residual_bar(eigen_count + response_count + 1:)
        end if
        status = 0
    end subroutine assemble_linear_response_eigen_cross_residual_vjp

    logical function valid_inputs( &
            stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, metadata, residual) result(valid)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenvalue
        complex(dp), intent(in) :: eigenstate(:), response_matrix(:, :)
        complex(dp), intent(in) :: response_coupling(:, :), response_state(:)
        complex(dp), intent(in) :: response_source(:), shielding_matrix(:, :)
        complex(dp), intent(in) :: shielding_target(:)
        type(linear_response_cross_metadata_t), intent(in) :: metadata
        complex(dp), intent(in) :: residual(:)
        integer :: eigen_count, response_count, shielding_count, metadata_status

        valid = .false.
        if (.not. validate_linear_response_cross_metadata(metadata, metadata_status)) return
        eigen_count = size(eigenstate)
        response_count = size(response_state)
        shielding_count = size(shielding_target)
        if (eigen_count < 1 .or. response_count < 1 .or. shielding_count < 0) return
        if (size(stiffness, 1) /= eigen_count .or. size(stiffness, 2) /= eigen_count .or. &
            any(shape(mass) /= shape(stiffness))) return
        if (size(response_matrix, 1) /= response_count .or. &
            size(response_matrix, 2) /= response_count .or. &
            size(response_coupling, 1) /= response_count .or. &
            size(response_coupling, 2) /= eigen_count .or. &
            size(response_source) /= response_count) return
        if (size(shielding_matrix, 1) /= shielding_count .or. &
            size(shielding_matrix, 2) /= response_count) return
        if (size(residual) /= eigen_count + response_count + shielding_count) return
        if (size(metadata%retained_indices) /= response_count .or. &
            maxval(metadata%retained_indices) > response_count .or. &
            metadata%shielding_trace_count /= shielding_count) return
        valid = finite_complex(stiffness) .and. finite_complex(mass) .and. &
            finite_complex(eigenvalue) .and. finite_complex(eigenstate) .and. &
            finite_complex(response_matrix) .and. finite_complex(response_coupling) .and. &
            finite_complex(response_state) .and. finite_complex(response_source) .and. &
            finite_complex(shielding_matrix) .and. finite_complex(shielding_target)
    end function valid_inputs

    logical function valid_tangent_inputs( &
            stiffness, mass, eigenstate, response_matrix, response_coupling, &
            response_state, response_source, shielding_matrix, shielding_target, &
            stiffness_dot, mass_dot, eigenstate_dot, response_matrix_dot, &
            response_coupling_dot, response_state_dot, response_source_dot, &
            shielding_matrix_dot, shielding_target_dot) result(valid)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenstate(:)
        complex(dp), intent(in) :: response_matrix(:, :), response_coupling(:, :)
        complex(dp), intent(in) :: response_state(:), response_source(:)
        complex(dp), intent(in) :: shielding_matrix(:, :), shielding_target(:)
        complex(dp), intent(in) :: stiffness_dot(:, :), mass_dot(:, :)
        complex(dp), intent(in) :: eigenstate_dot(:), response_matrix_dot(:, :)
        complex(dp), intent(in) :: response_coupling_dot(:, :)
        complex(dp), intent(in) :: response_state_dot(:), response_source_dot(:)
        complex(dp), intent(in) :: shielding_matrix_dot(:, :), shielding_target_dot(:)

        valid = all(shape(stiffness_dot) == shape(stiffness)) .and. &
            all(shape(mass_dot) == shape(mass)) .and. size(eigenstate_dot) == size(eigenstate) .and. &
            all(shape(response_matrix_dot) == shape(response_matrix)) .and. &
            all(shape(response_coupling_dot) == shape(response_coupling)) .and. &
            size(response_state_dot) == size(response_state) .and. &
            size(response_source_dot) == size(response_source) .and. &
            all(shape(shielding_matrix_dot) == shape(shielding_matrix)) .and. &
            size(shielding_target_dot) == size(shielding_target) .and. &
            finite_complex(stiffness_dot) .and. finite_complex(mass_dot) .and. &
            finite_complex(eigenstate_dot) .and. finite_complex(response_matrix_dot) .and. &
            finite_complex(response_coupling_dot) .and. finite_complex(response_state_dot) .and. &
            finite_complex(response_source_dot) .and. finite_complex(shielding_matrix_dot) .and. &
            finite_complex(shielding_target_dot)
    end function valid_tangent_inputs

    logical function valid_vjp_outputs( &
            stiffness, mass, eigenstate, response_matrix, response_coupling, &
            response_state, response_source, shielding_matrix, shielding_target, &
            stiffness_bar, mass_bar, eigenstate_bar, response_matrix_bar, &
            response_coupling_bar, response_state_bar, response_source_bar, &
            shielding_matrix_bar, shielding_target_bar) result(valid)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenstate(:)
        complex(dp), intent(in) :: response_matrix(:, :), response_coupling(:, :)
        complex(dp), intent(in) :: response_state(:), response_source(:)
        complex(dp), intent(in) :: shielding_matrix(:, :), shielding_target(:)
        complex(dp), intent(in) :: stiffness_bar(:, :), mass_bar(:, :)
        complex(dp), intent(in) :: eigenstate_bar(:), response_matrix_bar(:, :)
        complex(dp), intent(in) :: response_coupling_bar(:, :)
        complex(dp), intent(in) :: response_state_bar(:), response_source_bar(:)
        complex(dp), intent(in) :: shielding_matrix_bar(:, :), shielding_target_bar(:)

        valid = all(shape(stiffness_bar) == shape(stiffness)) .and. &
            all(shape(mass_bar) == shape(mass)) .and. size(eigenstate_bar) == size(eigenstate) .and. &
            all(shape(response_matrix_bar) == shape(response_matrix)) .and. &
            all(shape(response_coupling_bar) == shape(response_coupling)) .and. &
            size(response_state_bar) == size(response_state) .and. &
            size(response_source_bar) == size(response_source) .and. &
            all(shape(shielding_matrix_bar) == shape(shielding_matrix)) .and. &
            size(shielding_target_bar) == size(shielding_target)
    end function valid_vjp_outputs

    subroutine rank_one_product(matrix, left, right)
        complex(dp), intent(out) :: matrix(:, :)
        complex(dp), intent(in) :: left(:), right(:)
        integer :: row, column

        do column = 1, size(right)
            do row = 1, size(left)
                matrix(row, column) = left(row)*conjg(right(column))
            end do
        end do
    end subroutine rank_one_product

    subroutine clear_metadata(metadata)
        type(linear_response_cross_metadata_t), intent(inout) :: metadata

        if (allocated(metadata%retained_indices)) deallocate(metadata%retained_indices)
        if (allocated(metadata%response_labels)) deallocate(metadata%response_labels)
        metadata%retained_response_count = 0
        metadata%shielding_trace_count = 0
    end subroutine clear_metadata

    subroutine assign_linear_response_cross_metadata(lhs, rhs)
        class(linear_response_cross_metadata_t), intent(out) :: lhs
        type(linear_response_cross_metadata_t), intent(in) :: rhs

        call clear_metadata(lhs)
        lhs%schema_version = rhs%schema_version
        lhs%normalization_label = rhs%normalization_label
        lhs%provenance = rhs%provenance
        lhs%normalization_scale = rhs%normalization_scale
        lhs%ideal_shielding = rhs%ideal_shielding
        lhs%retained_response_count = rhs%retained_response_count
        lhs%shielding_trace_count = rhs%shielding_trace_count
        if (allocated(rhs%retained_indices)) allocate(lhs%retained_indices, source=rhs%retained_indices)
        if (allocated(rhs%response_labels)) allocate(lhs%response_labels, source=rhs%response_labels)
    end subroutine assign_linear_response_cross_metadata

    logical function unique_indices(values) result(unique)
        integer, intent(in) :: values(:)
        integer :: first, second

        unique = .true.
        do first = 1, size(values)
            do second = first + 1, size(values)
                if (values(first) == values(second)) unique = .false.
            end do
        end do
    end function unique_indices

    logical function unique_labels(values) result(unique)
        character(len=*), intent(in) :: values(:)
        integer :: first, second

        unique = .true.
        do first = 1, size(values)
            if (len_trim(values(first)) < 1) unique = .false.
            do second = first + 1, size(values)
                if (trim(values(first)) == trim(values(second))) unique = .false.
            end do
        end do
    end function unique_labels

    logical function finite_complex_scalar(value) result(valid)
        complex(dp), intent(in) :: value

        valid = ieee_is_finite(real(value, dp)) .and. ieee_is_finite(aimag(value))
    end function finite_complex_scalar

    logical function finite_complex_vector(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_vector

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

end module fortfem_linear_response_cross
