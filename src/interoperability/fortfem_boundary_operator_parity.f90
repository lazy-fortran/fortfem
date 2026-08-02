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

    public :: evaluate_boundary_operator_parity
    public :: validate_boundary_operator_parity

contains

    subroutine evaluate_boundary_operator_parity( &
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
    end subroutine evaluate_boundary_operator_parity

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
