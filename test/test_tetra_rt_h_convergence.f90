program test_tetra_rt_h_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_rt, &
        generate_structured_tetra_box_mesh, initialize_tetra_rt, &
        interpolate_physical_tetra_rt, map_tetra_rt_contravariant, &
        tetra_duffy_quadrature, tetra_rt_dof_count, tetra_rt_t
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    implicit none

    real(dp), parameter :: unit_bounds(3, 2) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [3, 2])
    real(dp) :: errors(2, 2), rates(2)
    integer :: degree
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 4
        call measure_error(2, degree, errors(:, 1))
        call measure_error(4, degree, errors(:, 2))
        rates = log(errors(:, 1)/errors(:, 2))/log(2.0_dp)
        write (*, '(a,i0,a,2(es12.4,1x),a,2(f8.3,1x))') &
            "Hdiv degree ", degree, " fine errors ", errors(:, 2), &
            " rates ", rates
        call record_condition( &
            all(errors(:, 2) < errors(:, 1)), &
            "tetrahedral Hdiv interpolation improves under refinement")
        call record_condition( &
            rates(1) > real(degree + 1, dp) - 0.5_dp, &
            "tetrahedral Hdiv field reaches its expected h order")
        call record_condition( &
            rates(2) > real(degree + 1, dp) - 0.6_dp, &
            "tetrahedral Hdiv divergence reaches its expected h order")
    end do

    call check_summary("tetrahedral Hdiv h convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine measure_error(refinement, degree, error)
        integer, intent(in) :: refinement, degree
        real(dp), intent(out) :: error(2)

        type(tetra_rt_t) :: basis
        integer, allocatable :: tetrahedra(:, :)
        real(dp), allocatable :: divergences(:), dofs(:)
        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :), values(:, :)
        real(dp), allocatable :: vertices(:, :), weights(:), x(:), y(:), z(:)
        real(dp) :: divergence_value, jacobian(3, 3), physical_point(3)
        real(dp) :: reference_point(3), value(3)
        integer :: cell, dof_count, node, status

        call generate_structured_tetra_box_mesh( &
            unit_bounds, [refinement, refinement, refinement], vertices, &
            tetrahedra, status)
        if (status /= 0) error stop "Hdiv box mesh generation failed"
        call initialize_tetra_rt(degree, basis, status)
        if (status /= 0) error stop "Hdiv basis initialization failed"
        dof_count = tetra_rt_dof_count(basis)
        allocate ( &
            values(3, dof_count), divergences(dof_count), &
            physical_values(3, dof_count), physical_divergences(dof_count))
        call tetra_duffy_quadrature( &
            2*degree + 8, x, y, z, weights, status)
        if (status /= 0) error stop "Hdiv error quadrature failed"

        error = 0.0_dp
        do cell = 1, size(tetrahedra, 2)
            jacobian(:, 1) = vertices(:, tetrahedra(2, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 2) = vertices(:, tetrahedra(3, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 3) = vertices(:, tetrahedra(4, cell)) - &
                vertices(:, tetrahedra(1, cell))
            call interpolate_physical_tetra_rt( &
                basis, vertices(:, tetrahedra(:, cell)), exact_field, dofs, &
                status)
            if (status /= 0) error stop "physical Hdiv interpolation failed"
            do node = 1, size(weights)
                reference_point = [x(node), y(node), z(node)]
                call evaluate_tetra_rt( &
                    basis, reference_point, values, divergences, status)
                if (status /= 0) error stop "Hdiv basis evaluation failed"
                call map_tetra_rt_contravariant( &
                    jacobian, values, divergences, physical_values, &
                    physical_divergences, status)
                if (status /= 0) error stop "Hdiv Piola evaluation failed"
                physical_point = vertices(:, tetrahedra(1, cell)) + &
                    matmul(jacobian, reference_point)
                value = matmul(physical_values, dofs)
                divergence_value = dot_product(physical_divergences, dofs)
                error(1) = error(1) + det3(jacobian)*weights(node)* &
                    sum((value - exact_field_value(physical_point))**2)
                error(2) = error(2) + det3(jacobian)*weights(node)* &
                    (divergence_value - exact_divergence(physical_point))**2
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

        value = sin(pi*point)
    end function exact_field_value

    pure function exact_divergence(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = pi*sum(cos(pi*point))
    end function exact_divergence

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_rt_h_convergence
