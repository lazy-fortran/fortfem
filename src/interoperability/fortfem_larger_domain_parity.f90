module fortfem_larger_domain_parity
    !! Common-interior comparison for two open-boundary domain truncations.
    !!
    !! A caller supplies two finite complex fields evaluated at the same
    !! physical samples.  The fields may come from FEM, BEM, DtN, PML, or
    !! another boundary client and may even use the same backend at two
    !! different artificial-boundary distances.  No mesh, kernel, solver, or
    !! application file format is owned here.  The contract only validates
    !! shared physical metadata and reports a weighted field difference plus
    !! the distance increment needed to interpret a larger-domain control.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_boundary_operator_contract, only: &
        boundary_operator_contract_t, validate_boundary_operator_contract
    implicit none
    private

    type, public :: larger_domain_parity_t
        character(len=32) :: schema_version = "fortfem-larger-domain-parity-1"
        integer :: component_count = 0
        integer :: sample_count = 0
        integer :: inner_backend_kind = 0
        integer :: outer_backend_kind = 0
        real(dp) :: inner_boundary_distance = 0.0_dp
        real(dp) :: outer_boundary_distance = 0.0_dp
        real(dp) :: distance_increase = 0.0_dp
        real(dp) :: distance_ratio = 0.0_dp
        real(dp) :: comparison_norm = 0.0_dp
        real(dp) :: absolute_difference = 0.0_dp
        real(dp) :: relative_difference = 0.0_dp
        real(dp) :: relative_difference_per_distance = 0.0_dp
        real(dp) :: absolute_tolerance = 0.0_dp
        real(dp) :: relative_tolerance = 0.0_dp
        logical :: farther_boundary = .false.
        logical :: within_tolerance = .false.
        character(len=64) :: topology_id = ""
        character(len=64) :: equation = ""
        character(len=64) :: space = ""
        character(len=32) :: units = ""
        character(len=64) :: normalization = ""
        character(len=128) :: provenance = ""
    end type larger_domain_parity_t

    public :: evaluate_larger_domain_parity
    public :: validate_larger_domain_parity

contains

    subroutine evaluate_larger_domain_parity( &
            inner_field, outer_field, weights, inner_contract, outer_contract, &
            inner_boundary_distance, outer_boundary_distance, absolute_tolerance, &
            relative_tolerance, report, status)
        !! Compare two fields at a shared physical sample set.
        !!
        !! The relative difference is symmetric: its denominator is the
        !! maximum of the two weighted field norms.  This avoids designating
        !! either domain as an exact solution.  `relative_difference_per_distance`
        !! is a finite-difference trend diagnostic, not an observed order; an
        !! order requires at least three domain distances or an exact oracle.
        complex(dp), intent(in) :: inner_field(:, :), outer_field(:, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: inner_contract
        type(boundary_operator_contract_t), intent(in) :: outer_contract
        real(dp), intent(in) :: inner_boundary_distance, outer_boundary_distance
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        type(larger_domain_parity_t), intent(out) :: report
        integer, intent(out) :: status

        integer :: sample, component, contract_status
        real(dp) :: inner_squared, outer_squared, difference_squared
        complex(dp) :: difference

        report = larger_domain_parity_t()
        status = 1
        if (size(inner_field, 1) < 1 .or. size(inner_field, 2) < 1) return
        if (size(outer_field, 1) /= size(inner_field, 1) .or. &
            size(outer_field, 2) /= size(inner_field, 2)) return
        if (size(weights) /= size(inner_field, 2)) return
        if (.not. finite_complex_2d(inner_field) .or. &
            .not. finite_complex_2d(outer_field) .or. .not. finite_real(weights)) return
        if (any(weights <= 0.0_dp)) return
        if (.not. ieee_is_finite(inner_boundary_distance) .or. &
            .not. ieee_is_finite(outer_boundary_distance) .or. &
            inner_boundary_distance <= 0.0_dp .or. &
            outer_boundary_distance <= inner_boundary_distance) return
        if (.not. ieee_is_finite(absolute_tolerance) .or. &
            .not. ieee_is_finite(relative_tolerance) .or. &
            absolute_tolerance < 0.0_dp .or. relative_tolerance < 0.0_dp) return
        if (.not. validate_boundary_operator_contract(inner_contract, contract_status)) return
        if (.not. validate_boundary_operator_contract(outer_contract, contract_status)) return
        if (.not. inner_contract%complex_valued .or. .not. outer_contract%complex_valued) return
        if (trim(inner_contract%topology_id) /= trim(outer_contract%topology_id) .or. &
            trim(inner_contract%equation) /= trim(outer_contract%equation) .or. &
            trim(inner_contract%space) /= trim(outer_contract%space) .or. &
            trim(inner_contract%units) /= trim(outer_contract%units) .or. &
            trim(inner_contract%normalization) /= trim(outer_contract%normalization)) return

        report%component_count = size(inner_field, 1)
        report%sample_count = size(inner_field, 2)
        report%inner_backend_kind = inner_contract%backend_kind
        report%outer_backend_kind = outer_contract%backend_kind
        report%inner_boundary_distance = inner_boundary_distance
        report%outer_boundary_distance = outer_boundary_distance
        report%distance_increase = outer_boundary_distance - inner_boundary_distance
        report%distance_ratio = outer_boundary_distance / inner_boundary_distance
        report%absolute_tolerance = absolute_tolerance
        report%relative_tolerance = relative_tolerance
        report%farther_boundary = .true.
        report%topology_id = inner_contract%topology_id
        report%equation = inner_contract%equation
        report%space = inner_contract%space
        report%units = inner_contract%units
        report%normalization = inner_contract%normalization
        report%provenance = inner_contract%provenance

        inner_squared = 0.0_dp
        outer_squared = 0.0_dp
        difference_squared = 0.0_dp
        do sample = 1, report%sample_count
            do component = 1, report%component_count
                inner_squared = inner_squared + weights(sample)*abs( &
                    inner_field(component, sample))**2
                outer_squared = outer_squared + weights(sample)*abs( &
                    outer_field(component, sample))**2
                difference = outer_field(component, sample) - &
                    inner_field(component, sample)
                difference_squared = difference_squared + weights(sample)*abs(difference)**2
            end do
        end do
        report%comparison_norm = max(sqrt(inner_squared), sqrt(outer_squared))
        report%absolute_difference = sqrt(difference_squared)
        report%relative_difference = report%absolute_difference / max( &
            report%comparison_norm, tiny(1.0_dp))
        report%relative_difference_per_distance = report%relative_difference / &
            report%distance_increase
        report%within_tolerance = report%absolute_difference <= absolute_tolerance .or. &
            report%relative_difference <= relative_tolerance
        if (.not. ieee_is_finite(report%comparison_norm) .or. &
            .not. ieee_is_finite(report%absolute_difference) .or. &
            .not. ieee_is_finite(report%relative_difference) .or. &
            .not. ieee_is_finite(report%relative_difference_per_distance)) then
            report = larger_domain_parity_t()
            return
        end if
        status = 0
    end subroutine evaluate_larger_domain_parity

    logical function validate_larger_domain_parity(report, status) result(valid)
        type(larger_domain_parity_t), intent(in) :: report
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (report%schema_version /= "fortfem-larger-domain-parity-1") return
        if (report%component_count < 1 .or. report%sample_count < 1) return
        if (report%inner_backend_kind < 1 .or. report%inner_backend_kind > 8 .or. &
            report%outer_backend_kind < 1 .or. report%outer_backend_kind > 8) return
        if (.not. ieee_is_finite(report%inner_boundary_distance) .or. &
            .not. ieee_is_finite(report%outer_boundary_distance) .or. &
            .not. ieee_is_finite(report%distance_increase) .or. &
            .not. ieee_is_finite(report%distance_ratio) .or. &
            .not. ieee_is_finite(report%comparison_norm) .or. &
            .not. ieee_is_finite(report%absolute_difference) .or. &
            .not. ieee_is_finite(report%relative_difference) .or. &
            .not. ieee_is_finite(report%relative_difference_per_distance) .or. &
            .not. ieee_is_finite(report%absolute_tolerance) .or. &
            .not. ieee_is_finite(report%relative_tolerance)) return
        if (report%inner_boundary_distance <= 0.0_dp .or. &
            report%outer_boundary_distance <= report%inner_boundary_distance .or. &
            report%distance_increase <= 0.0_dp .or. report%distance_ratio <= 1.0_dp .or. &
            report%comparison_norm < 0.0_dp .or. report%absolute_difference < 0.0_dp .or. &
            report%relative_difference < 0.0_dp .or. &
            report%relative_difference_per_distance < 0.0_dp .or. &
            report%absolute_tolerance < 0.0_dp .or. report%relative_tolerance < 0.0_dp) return
        if (.not. report%farther_boundary .or. len_trim(report%topology_id) == 0 .or. &
            len_trim(report%equation) == 0 .or. len_trim(report%space) == 0 .or. &
            len_trim(report%units) == 0 .or. len_trim(report%normalization) == 0 .or. &
            len_trim(report%provenance) == 0) return
        status = 0
        valid = .true.
    end function validate_larger_domain_parity

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

end module fortfem_larger_domain_parity
