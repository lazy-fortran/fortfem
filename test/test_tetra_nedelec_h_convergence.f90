program test_tetra_nedelec_h_convergence
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_structured_tetra_box_mesh
    use fortfem_feec, only: evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, map_tetra_nedelec_covariant, &
        tetra_duffy_quadrature, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_interpolation, only: &
        interpolate_physical_tetra_nedelec
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    implicit none

    real(dp), parameter :: unit_bounds(3, 2) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [3, 2])
    real(dp) :: errors(2, 2), rates(2)
    integer :: order
    logical :: all_passed

    all_passed = .true.
    do order = 1, 4
        call measure_error(2, order, errors(:, 1))
        call measure_error(4, order, errors(:, 2))
        rates = log(errors(:, 1)/errors(:, 2))/log(2.0_dp)
        write (*, '(a,i0,a,2(es12.4,1x),a,2(f8.3,1x))') &
            "Hcurl order ", order, " fine errors ", errors(:, 2), &
            " rates ", rates
        call record_condition( &
            all(errors(:, 2) < errors(:, 1)), &
            "tetrahedral Hcurl interpolation improves under refinement")
        call record_condition( &
            rates(1) > real(order, dp) - 0.5_dp, &
            "tetrahedral Hcurl field reaches its expected h order")
        call record_condition( &
            rates(2) > real(order, dp) - 0.6_dp, &
            "tetrahedral Hcurl curl reaches its expected h order")
    end do

    call check_summary("tetrahedral Hcurl h convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine measure_error(refinement, order, error)
        integer, intent(in) :: refinement, order
        real(dp), intent(out) :: error(2)

        type(tetra_nedelec_first_kind_t) :: basis
        integer, allocatable :: tetrahedra(:, :)
        real(dp), allocatable :: dofs(:), physical_curls(:, :)
        real(dp), allocatable :: physical_values(:, :), reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :), vertices(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: curl_value(3), jacobian(3, 3), physical_point(3)
        real(dp) :: reference_point(3), value(3)
        integer :: cell, dof_count, node, status

        call generate_structured_tetra_box_mesh( &
            unit_bounds, [refinement, refinement, refinement], vertices, &
            tetrahedra, status)
        if (status /= 0) error stop "Hcurl box mesh generation failed"
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) error stop "Hcurl basis initialization failed"
        dof_count = tetra_nedelec_dof_count(basis)
        allocate ( &
            reference_values(3, dof_count), reference_curls(3, dof_count), &
            physical_values(3, dof_count), physical_curls(3, dof_count))
        call tetra_duffy_quadrature( &
            2*order + 6, x, y, z, weights, status)
        if (status /= 0) error stop "Hcurl error quadrature failed"

        error = 0.0_dp
        do cell = 1, size(tetrahedra, 2)
            jacobian(:, 1) = vertices(:, tetrahedra(2, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 2) = vertices(:, tetrahedra(3, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 3) = vertices(:, tetrahedra(4, cell)) - &
                vertices(:, tetrahedra(1, cell))
            call interpolate_physical_tetra_nedelec( &
                basis, vertices(:, tetrahedra(:, cell)), exact_field, dofs, &
                status)
            if (status /= 0) error stop "physical Hcurl interpolation failed"
            do node = 1, size(weights)
                reference_point = [x(node), y(node), z(node)]
                call evaluate_tetra_nedelec_first_kind( &
                    basis, reference_point, reference_values, &
                    reference_curls, status)
                if (status /= 0) error stop "Hcurl basis evaluation failed"
                call map_tetra_nedelec_covariant( &
                    jacobian, reference_values, reference_curls, &
                    physical_values, physical_curls, status)
                if (status /= 0) error stop "Hcurl Piola evaluation failed"
                physical_point = vertices(:, tetrahedra(1, cell)) + &
                    matmul(jacobian, reference_point)
                value = matmul(physical_values, dofs)
                curl_value = matmul(physical_curls, dofs)
                error(1) = error(1) + det3(jacobian)*weights(node)* &
                    sum((value - exact_field_value(physical_point))**2)
                error(2) = error(2) + det3(jacobian)*weights(node)* &
                    sum((curl_value - exact_curl(physical_point))**2)
            end do
        end do
        error = sqrt(error)
    end subroutine measure_error

    pure subroutine exact_field(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)

        value = exact_field_value([x, y, z])
    end subroutine exact_field

    pure function exact_field_value(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [ &
            sin(pi*point(2)), sin(pi*point(3)), sin(pi*point(1))]
    end function exact_field_value

    pure function exact_curl(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = -pi*[ &
            cos(pi*point(3)), cos(pi*point(1)), cos(pi*point(2))]
    end function exact_curl

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_h_convergence
