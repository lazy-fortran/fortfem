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
    public :: evaluate_larger_domain_parity_jvp
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

    subroutine evaluate_larger_domain_parity_jvp( &
            inner_field, outer_field, weights, inner_contract, outer_contract, &
            inner_boundary_distance, outer_boundary_distance, absolute_tolerance, &
            relative_tolerance, inner_field_dot, outer_field_dot, weights_dot, &
            inner_boundary_distance_dot, outer_boundary_distance_dot, &
            comparison_norm_dot, absolute_difference_dot, relative_difference_dot, &
            relative_difference_per_distance_dot, distance_increase_dot, &
            distance_ratio_dot, status)
        !! Fixed-topology JVP of the common-interior parity metrics.
        !!
        !! The supplied fields may be nonzero scattering states from distinct
        !! FEM, BEM, DtN, or PML clients.  Integer/backend metadata and the
        !! tolerance decision remain fixed.  A JVP is rejected at the
        !! nondifferentiable equal-norm or zero-difference branches.
        complex(dp), intent(in) :: inner_field(:, :), outer_field(:, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: inner_contract
        type(boundary_operator_contract_t), intent(in) :: outer_contract
        real(dp), intent(in) :: inner_boundary_distance, outer_boundary_distance
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        complex(dp), intent(in) :: inner_field_dot(:, :), outer_field_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(in) :: inner_boundary_distance_dot
        real(dp), intent(in) :: outer_boundary_distance_dot
        real(dp), intent(out) :: comparison_norm_dot, absolute_difference_dot
        real(dp), intent(out) :: relative_difference_dot
        real(dp), intent(out) :: relative_difference_per_distance_dot
        real(dp), intent(out) :: distance_increase_dot, distance_ratio_dot
        integer, intent(out) :: status

        type(larger_domain_parity_t) :: report
        integer :: sample, component
        real(dp) :: inner_squared, outer_squared, difference_squared
        real(dp) :: inner_squared_dot, outer_squared_dot, difference_squared_dot
        real(dp) :: comparison_norm, absolute_difference
        complex(dp) :: inner_value, outer_value, difference, difference_dot

        comparison_norm_dot = 0.0_dp
        absolute_difference_dot = 0.0_dp
        relative_difference_dot = 0.0_dp
        relative_difference_per_distance_dot = 0.0_dp
        distance_increase_dot = 0.0_dp
        distance_ratio_dot = 0.0_dp
        call evaluate_larger_domain_parity( &
            inner_field, outer_field, weights, inner_contract, outer_contract, &
            inner_boundary_distance, outer_boundary_distance, absolute_tolerance, &
            relative_tolerance, report, status)
        if (status /= 0) return
        if (.not. valid_jvp_directions( &
                inner_field, outer_field, weights, inner_field_dot, &
                outer_field_dot, weights_dot)) then
            status = 1
            return
        end if

        inner_squared = 0.0_dp
        outer_squared = 0.0_dp
        difference_squared = 0.0_dp
        inner_squared_dot = 0.0_dp
        outer_squared_dot = 0.0_dp
        difference_squared_dot = 0.0_dp
        do sample = 1, size(inner_field, 2)
            do component = 1, size(inner_field, 1)
                inner_value = inner_field(component, sample)
                outer_value = outer_field(component, sample)
                difference = outer_value - inner_value
                difference_dot = outer_field_dot(component, sample) - &
                    inner_field_dot(component, sample)
                inner_squared = inner_squared + weights(sample)*abs(inner_value)**2
                outer_squared = outer_squared + weights(sample)*abs(outer_value)**2
                difference_squared = difference_squared + weights(sample)*abs(difference)**2
                inner_squared_dot = inner_squared_dot + weights_dot(sample)*abs(inner_value)**2 + &
                    2.0_dp*weights(sample)*real(conjg(inner_value)* &
                    inner_field_dot(component, sample), dp)
                outer_squared_dot = outer_squared_dot + weights_dot(sample)*abs(outer_value)**2 + &
                    2.0_dp*weights(sample)*real(conjg(outer_value)* &
                    outer_field_dot(component, sample), dp)
                difference_squared_dot = difference_squared_dot + &
                    weights_dot(sample)*abs(difference)**2 + &
                    2.0_dp*weights(sample)*real(conjg(difference)*difference_dot, dp)
            end do
        end do
        if (abs(inner_squared - outer_squared) <= tiny(1.0_dp) .or. &
            difference_squared <= tiny(1.0_dp)) then
            status = 1
            return
        end if
        comparison_norm = report%comparison_norm
        absolute_difference = report%absolute_difference
        if (outer_squared > inner_squared) then
            comparison_norm_dot = outer_squared_dot/(2.0_dp*comparison_norm)
        else
            comparison_norm_dot = inner_squared_dot/(2.0_dp*comparison_norm)
        end if
        absolute_difference_dot = difference_squared_dot/(2.0_dp*absolute_difference)
        relative_difference_dot = (absolute_difference_dot*comparison_norm - &
            absolute_difference*comparison_norm_dot)/comparison_norm**2
        distance_increase_dot = outer_boundary_distance_dot - inner_boundary_distance_dot
        distance_ratio_dot = outer_boundary_distance_dot/inner_boundary_distance - &
            outer_boundary_distance*inner_boundary_distance_dot/ &
            inner_boundary_distance**2
        relative_difference_per_distance_dot = relative_difference_dot/ &
            report%distance_increase - report%relative_difference*distance_increase_dot/ &
            report%distance_increase**2
        if (.not. ieee_is_finite(comparison_norm_dot) .or. &
            .not. ieee_is_finite(absolute_difference_dot) .or. &
            .not. ieee_is_finite(relative_difference_dot) .or. &
            .not. ieee_is_finite(relative_difference_per_distance_dot) .or. &
            .not. ieee_is_finite(distance_increase_dot) .or. &
            .not. ieee_is_finite(distance_ratio_dot)) then
            comparison_norm_dot = 0.0_dp
            absolute_difference_dot = 0.0_dp
            relative_difference_dot = 0.0_dp
            relative_difference_per_distance_dot = 0.0_dp
            distance_increase_dot = 0.0_dp
            distance_ratio_dot = 0.0_dp
            status = 1
            return
        end if
        status = 0
    end subroutine evaluate_larger_domain_parity_jvp

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

    logical function valid_jvp_directions( &
            inner_field, outer_field, weights, inner_field_dot, &
            outer_field_dot, weights_dot) result(valid)
        complex(dp), intent(in) :: inner_field(:, :), outer_field(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: inner_field_dot(:, :), outer_field_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)

        valid = size(inner_field_dot, 1) == size(inner_field, 1) .and. &
            size(inner_field_dot, 2) == size(inner_field, 2) .and. &
            size(outer_field_dot, 1) == size(outer_field, 1) .and. &
            size(outer_field_dot, 2) == size(outer_field, 2) .and. &
            size(weights_dot) == size(weights) .and. &
            finite_complex_2d(inner_field_dot) .and. &
            finite_complex_2d(outer_field_dot) .and. finite_real(weights_dot)
    end function valid_jvp_directions

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
