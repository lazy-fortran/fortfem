module fortfem_boundary_operator_parity
    !! Common physical-sample parity diagnostics for exterior operators.
    !!
    !! A caller supplies one manufactured reference field and the values
    !! reconstructed by FEM, BEM, DtN, PML, or another boundary client on the
    !! same physical sample set.  The contract checks that all metadata refer
    !! to one fixed topology and reports weighted absolute/relative errors.
    !! It deliberately contains no mesh, kernel, solver, or application file
    !! format and therefore can be used by neutral open-boundary clients.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_boundary_operator_contract, only: &
        boundary_operator_contract_t, validate_boundary_operator_contract
    implicit none
    private

    type, public :: boundary_operator_parity_t
        character(len=32) :: schema_version = "fortfem-boundary-parity-1"
        integer :: component_count = 0
        integer :: sample_count = 0
        integer :: backend_count = 0
        real(dp) :: reference_norm = 0.0_dp
        real(dp) :: absolute_tolerance = 0.0_dp
        real(dp) :: relative_tolerance = 0.0_dp
        character(len=64) :: topology_id = ""
        character(len=128) :: provenance = ""
        integer, allocatable :: backend_kind(:)
        real(dp), allocatable :: absolute_error(:)
        real(dp), allocatable :: relative_error(:)
        logical, allocatable :: within_tolerance(:)
    end type boundary_operator_parity_t

    public :: compare_boundary_operator_parity
    public :: compare_boundary_operator_parity_jvp
    public :: compare_boundary_operator_parity_vjp
    public :: validate_boundary_operator_parity

contains

    subroutine compare_boundary_operator_parity( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, report, status)
        complex(dp), intent(in) :: reference(:, :), candidates(:, :, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: contracts(:)
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        type(boundary_operator_parity_t), intent(out) :: report
        integer, intent(out) :: status

        integer :: backend, component, sample, contract_status
        real(dp) :: reference_squared, error_squared
        complex(dp) :: difference

        report = boundary_operator_parity_t()
        status = 1
        if (size(reference, 1) < 1 .or. size(reference, 2) < 1) return
        if (size(candidates, 1) /= size(reference, 1) .or. &
            size(candidates, 2) /= size(reference, 2)) return
        if (size(candidates, 3) /= size(contracts) .or. size(contracts) < 1) return
        if (size(weights) /= size(reference, 2)) return
        if (.not. finite_complex_2d(reference) .or. .not. finite_complex_3d(candidates)) return
        if (.not. finite_real(weights) .or. any(weights <= 0.0_dp)) return
        if (.not. ieee_is_finite(absolute_tolerance) .or. &
            .not. ieee_is_finite(relative_tolerance) .or. &
            absolute_tolerance < 0.0_dp .or. relative_tolerance < 0.0_dp) return

        do backend = 1, size(contracts)
            if (.not. validate_boundary_operator_contract(contracts(backend), contract_status)) return
            if (backend > 1) then
                if (trim(contracts(backend)%topology_id) /= &
                    trim(contracts(1)%topology_id)) return
                if (trim(contracts(backend)%equation) /= trim(contracts(1)%equation)) return
                if (trim(contracts(backend)%space) /= trim(contracts(1)%space)) return
                if (trim(contracts(backend)%units) /= trim(contracts(1)%units)) return
                if (trim(contracts(backend)%normalization) /= &
                    trim(contracts(1)%normalization)) return
            end if
            if (backend > 1) then
                if (any(contracts(backend)%backend_kind == &
                    contracts(1:backend-1)%backend_kind)) return
            end if
        end do

        report%component_count = size(reference, 1)
        report%sample_count = size(reference, 2)
        report%backend_count = size(contracts)
        report%absolute_tolerance = absolute_tolerance
        report%relative_tolerance = relative_tolerance
        report%topology_id = contracts(1)%topology_id
        report%provenance = contracts(1)%provenance
        allocate(report%backend_kind(report%backend_count), &
            report%absolute_error(report%backend_count), &
            report%relative_error(report%backend_count), &
            report%within_tolerance(report%backend_count))
        do backend = 1, report%backend_count
            report%backend_kind(backend) = contracts(backend)%backend_kind
        end do

        reference_squared = 0.0_dp
        do sample = 1, report%sample_count
            do component = 1, report%component_count
                reference_squared = reference_squared + weights(sample)*abs( &
                    reference(component, sample))**2
            end do
        end do
        report%reference_norm = sqrt(reference_squared)
        do backend = 1, report%backend_count
            error_squared = 0.0_dp
            do sample = 1, report%sample_count
                do component = 1, report%component_count
                    difference = candidates(component, sample, backend) - &
                        reference(component, sample)
                    error_squared = error_squared + weights(sample)*abs(difference)**2
                end do
            end do
            report%absolute_error(backend) = sqrt(error_squared)
            report%relative_error(backend) = report%absolute_error(backend)/max( &
                report%reference_norm, tiny(1.0_dp))
            report%within_tolerance(backend) = &
                report%absolute_error(backend) <= absolute_tolerance .or. &
                report%relative_error(backend) <= relative_tolerance
        end do
        status = 0
    end subroutine compare_boundary_operator_parity

    subroutine compare_boundary_operator_parity_jvp( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, reference_dot, candidates_dot, weights_dot, &
            reference_norm_dot, absolute_error_dot, relative_error_dot, status)
        !! Fixed-topology complex JVP of weighted parity metrics.
        !!
        !! The metadata and tolerance decisions remain fixed.  The common
        !! reference and every FEM/BEM/DtN/PML candidate may vary, as may the
        !! supplied quadrature weights.  Zero reference or candidate error
        !! norms are rejected because their Euclidean norm has no unique JVP.
        complex(dp), intent(in) :: reference(:, :), candidates(:, :, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: contracts(:)
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        complex(dp), intent(in) :: reference_dot(:, :), candidates_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: reference_norm_dot
        real(dp), intent(out) :: absolute_error_dot(:), relative_error_dot(:)
        integer, intent(out) :: status

        type(boundary_operator_parity_t) :: report
        integer :: backend, component, sample
        real(dp) :: reference_squared_dot, error_squared_dot
        complex(dp) :: difference, difference_dot

        reference_norm_dot = 0.0_dp
        absolute_error_dot = 0.0_dp
        relative_error_dot = 0.0_dp
        call compare_boundary_operator_parity( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, report, status)
        if (status /= 0) return
        if (.not. valid_jvp_directions( &
            reference, candidates, weights, reference_dot, candidates_dot, &
            weights_dot, absolute_error_dot, relative_error_dot)) then
            status = 1
            return
        end if
        if (report%reference_norm <= tiny(1.0_dp) .or. &
            any(report%absolute_error <= tiny(1.0_dp))) then
            status = 1
            return
        end if

        reference_squared_dot = 0.0_dp
        do sample = 1, size(reference, 2)
            do component = 1, size(reference, 1)
                reference_squared_dot = reference_squared_dot + &
                    weights_dot(sample)*abs(reference(component, sample))**2 + &
                    2.0_dp*weights(sample)*real(conjg(reference(component, sample))* &
                    reference_dot(component, sample), dp)
            end do
        end do
        reference_norm_dot = reference_squared_dot/(2.0_dp*report%reference_norm)
        do backend = 1, size(contracts)
            error_squared_dot = 0.0_dp
            do sample = 1, size(reference, 2)
                do component = 1, size(reference, 1)
                    difference = candidates(component, sample, backend) - &
                        reference(component, sample)
                    difference_dot = candidates_dot(component, sample, backend) - &
                        reference_dot(component, sample)
                    error_squared_dot = error_squared_dot + &
                        weights_dot(sample)*abs(difference)**2 + &
                        2.0_dp*weights(sample)*real(conjg(difference)*difference_dot, dp)
                end do
            end do
            absolute_error_dot(backend) = error_squared_dot/( &
                2.0_dp*report%absolute_error(backend))
            relative_error_dot(backend) = (absolute_error_dot(backend)* &
                report%reference_norm - report%absolute_error(backend)* &
                reference_norm_dot)/report%reference_norm**2
        end do
        if (.not. ieee_is_finite(reference_norm_dot) .or. &
            .not. finite_real(absolute_error_dot) .or. &
            .not. finite_real(relative_error_dot)) then
            reference_norm_dot = 0.0_dp
            absolute_error_dot = 0.0_dp
            relative_error_dot = 0.0_dp
            status = 1
            return
        end if
        status = 0
    end subroutine compare_boundary_operator_parity_jvp

    subroutine compare_boundary_operator_parity_vjp( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, reference_norm_bar, absolute_error_bar, &
            relative_error_bar, reference_bar, candidates_bar, weights_bar, status)
        !! Fixed-topology real-part complex adjoint of parity metrics.
        complex(dp), intent(in) :: reference(:, :), candidates(:, :, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: contracts(:)
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        real(dp), intent(in) :: reference_norm_bar
        real(dp), intent(in) :: absolute_error_bar(:), relative_error_bar(:)
        complex(dp), intent(out) :: reference_bar(:, :), candidates_bar(:, :, :)
        real(dp), intent(out) :: weights_bar(:)
        integer, intent(out) :: status

        type(boundary_operator_parity_t) :: report
        integer :: backend, component, sample
        real(dp) :: reference_bar_effective, error_bar_effective
        complex(dp) :: difference

        reference_bar = cmplx(0.0_dp, 0.0_dp, dp)
        candidates_bar = cmplx(0.0_dp, 0.0_dp, dp)
        weights_bar = 0.0_dp
        call compare_boundary_operator_parity( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, report, status)
        if (status /= 0) return
        if (.not. valid_vjp_inputs( &
            reference, candidates, weights, reference_norm_bar, &
            absolute_error_bar, relative_error_bar, reference_bar, &
            candidates_bar, weights_bar)) then
            status = 1
            return
        end if
        if (report%reference_norm <= tiny(1.0_dp) .or. &
            any(report%absolute_error <= tiny(1.0_dp))) then
            status = 1
            return
        end if
        reference_bar_effective = reference_norm_bar
        do backend = 1, size(contracts)
            reference_bar_effective = reference_bar_effective - &
                relative_error_bar(backend)*report%absolute_error(backend)/ &
                report%reference_norm**2
        end do
        do sample = 1, size(reference, 2)
            do component = 1, size(reference, 1)
                reference_bar(component, sample) = &
                    reference_bar_effective*weights(sample)*reference(component, sample)/ &
                    report%reference_norm
                weights_bar(sample) = weights_bar(sample) + &
                    reference_bar_effective*abs(reference(component, sample))**2/ &
                    (2.0_dp*report%reference_norm)
            end do
        end do
        do backend = 1, size(contracts)
            error_bar_effective = absolute_error_bar(backend) + &
                relative_error_bar(backend)/report%reference_norm
            do sample = 1, size(reference, 2)
                do component = 1, size(reference, 1)
                    difference = candidates(component, sample, backend) - &
                        reference(component, sample)
                    candidates_bar(component, sample, backend) = &
                        error_bar_effective*weights(sample)*difference/ &
                        report%absolute_error(backend)
                    reference_bar(component, sample) = reference_bar(component, sample) - &
                        error_bar_effective*weights(sample)*difference/ &
                        report%absolute_error(backend)
                    weights_bar(sample) = weights_bar(sample) + &
                        error_bar_effective*abs(difference)**2/ &
                        (2.0_dp*report%absolute_error(backend))
                end do
            end do
        end do
        if (.not. finite_complex_2d(reference_bar) .or. &
            .not. finite_complex_3d(candidates_bar) .or. &
            .not. finite_real(weights_bar)) then
            reference_bar = cmplx(0.0_dp, 0.0_dp, dp)
            candidates_bar = cmplx(0.0_dp, 0.0_dp, dp)
            weights_bar = 0.0_dp
            status = 1
            return
        end if
        status = 0
    end subroutine compare_boundary_operator_parity_vjp

    logical function validate_boundary_operator_parity(report, status) result(valid)
        type(boundary_operator_parity_t), intent(in) :: report
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (report%schema_version /= "fortfem-boundary-parity-1") return
        if (report%component_count < 1 .or. report%sample_count < 1 .or. &
            report%backend_count < 1) return
        if (.not. ieee_is_finite(report%reference_norm) .or. &
            .not. ieee_is_finite(report%absolute_tolerance) .or. &
            .not. ieee_is_finite(report%relative_tolerance)) return
        if (report%reference_norm < 0.0_dp .or. report%absolute_tolerance < 0.0_dp .or. &
            report%relative_tolerance < 0.0_dp) return
        if (len_trim(report%topology_id) == 0 .or. len_trim(report%provenance) == 0) return
        if (.not. allocated(report%backend_kind) .or. &
            .not. allocated(report%absolute_error) .or. &
            .not. allocated(report%relative_error) .or. &
            .not. allocated(report%within_tolerance)) return
        if (size(report%backend_kind) /= report%backend_count .or. &
            size(report%absolute_error) /= report%backend_count .or. &
            size(report%relative_error) /= report%backend_count .or. &
            size(report%within_tolerance) /= report%backend_count) return
        if (any(report%backend_kind < 1) .or. any(report%absolute_error < 0.0_dp) .or. &
            any(report%relative_error < 0.0_dp)) return
        if (.not. finite_real(report%absolute_error) .or. &
            .not. finite_real(report%relative_error)) return
        status = 0
        valid = .true.
    end function validate_boundary_operator_parity

    logical function valid_jvp_directions( &
            reference, candidates, weights, reference_dot, candidates_dot, &
            weights_dot, absolute_error_dot, relative_error_dot) result(valid)
        complex(dp), intent(in) :: reference(:, :), candidates(:, :, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: reference_dot(:, :), candidates_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(in) :: absolute_error_dot(:), relative_error_dot(:)

        valid = size(reference_dot, 1) == size(reference, 1) .and. &
            size(reference_dot, 2) == size(reference, 2) .and. &
            size(candidates_dot, 1) == size(candidates, 1) .and. &
            size(candidates_dot, 2) == size(candidates, 2) .and. &
            size(candidates_dot, 3) == size(candidates, 3) .and. &
            size(weights_dot) == size(weights) .and. &
            size(absolute_error_dot) == size(candidates, 3) .and. &
            size(relative_error_dot) == size(candidates, 3) .and. &
            finite_complex_2d(reference_dot) .and. &
            finite_complex_3d(candidates_dot) .and. finite_real(weights_dot)
    end function valid_jvp_directions

    logical function valid_vjp_inputs( &
            reference, candidates, weights, reference_norm_bar, absolute_error_bar, &
            relative_error_bar, reference_bar, candidates_bar, weights_bar) result(valid)
        complex(dp), intent(in) :: reference(:, :), candidates(:, :, :)
        real(dp), intent(in) :: weights(:), reference_norm_bar
        real(dp), intent(in) :: absolute_error_bar(:), relative_error_bar(:)
        complex(dp), intent(in) :: reference_bar(:, :), candidates_bar(:, :, :)
        real(dp), intent(in) :: weights_bar(:)

        valid = size(absolute_error_bar) == size(candidates, 3) .and. &
            size(relative_error_bar) == size(candidates, 3) .and. &
            size(reference_bar, 1) == size(reference, 1) .and. &
            size(reference_bar, 2) == size(reference, 2) .and. &
            size(candidates_bar, 1) == size(candidates, 1) .and. &
            size(candidates_bar, 2) == size(candidates, 2) .and. &
            size(candidates_bar, 3) == size(candidates, 3) .and. &
            size(weights_bar) == size(weights) .and. &
            ieee_is_finite(reference_norm_bar) .and. &
            finite_real(absolute_error_bar) .and. finite_real(relative_error_bar)
    end function valid_vjp_inputs

    logical function finite_complex_2d(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        integer :: component, sample

        valid = .true.
        do sample = 1, size(values, 2)
            do component = 1, size(values, 1)
                if (.not. ieee_is_finite(real(values(component, sample), dp)) .or. &
                    .not. ieee_is_finite(aimag(values(component, sample)))) then
                    valid = .false.
                    return
                end if
            end do
        end do
    end function finite_complex_2d

    logical function finite_complex_3d(values) result(valid)
        complex(dp), intent(in) :: values(:, :, :)

        integer :: component, sample, backend

        valid = .true.
        do backend = 1, size(values, 3)
            do sample = 1, size(values, 2)
                do component = 1, size(values, 1)
                    if (.not. ieee_is_finite(real(values(component, sample, backend), dp)) .or. &
                        .not. ieee_is_finite(aimag(values(component, sample, backend)))) then
                        valid = .false.
                        return
                    end if
                end do
            end do
        end do
    end function finite_complex_3d

    logical function finite_real(values) result(valid)
        real(dp), intent(in) :: values(:)

        integer :: index

        valid = .true.
        do index = 1, size(values)
            if (.not. ieee_is_finite(values(index))) then
                valid = .false.
                return
            end if
        end do
    end function finite_real

end module fortfem_boundary_operator_parity
