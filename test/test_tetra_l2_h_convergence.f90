program test_tetra_l2_h_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_discontinuous, &
        generate_structured_tetra_box_mesh, initialize_tetra_discontinuous, &
        project_physical_tetra_discontinuous, tetra_discontinuous_dof_count, &
        tetra_discontinuous_t, tetra_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    implicit none

    real(dp), parameter :: unit_bounds(3, 2) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [3, 2])
    real(dp) :: errors(2), rate
    integer :: degree
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 4
        call measure_error(2, degree, errors(1))
        call measure_error(4, degree, errors(2))
        rate = log(errors(1)/errors(2))/log(2.0_dp)
        write (*, '(a,i0,a,es12.4,a,f8.3)') &
            "L2 degree ", degree, " fine error ", errors(2), " rate ", rate
        call record_condition( &
            errors(2) < errors(1), &
            "tetrahedral L2 projection improves under refinement")
        call record_condition( &
            rate > real(degree + 1, dp) - 0.5_dp, &
            "tetrahedral L2 projection reaches its expected h order")
    end do

    call check_summary("tetrahedral L2 h convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine measure_error(refinement, degree, error)
        integer, intent(in) :: refinement, degree
        real(dp), intent(out) :: error

        type(tetra_discontinuous_t) :: basis
        integer, allocatable :: tetrahedra(:, :)
        real(dp), allocatable :: coefficients(:), values(:), vertices(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: jacobian(3, 3), physical_point(3)
        real(dp) :: reference_point(3), projected_value
        integer :: cell, dof_count, node, status

        call generate_structured_tetra_box_mesh( &
            unit_bounds, [refinement, refinement, refinement], vertices, &
            tetrahedra, status)
        if (status /= 0) error stop "L2 box mesh generation failed"
        call initialize_tetra_discontinuous(degree, basis, status)
        if (status /= 0) error stop "L2 basis initialization failed"
        dof_count = tetra_discontinuous_dof_count(basis)
        allocate(values(dof_count))
        call tetra_duffy_quadrature( &
            2*degree + 8, x, y, z, weights, status)
        if (status /= 0) error stop "L2 error quadrature failed"

        error = 0.0_dp
        do cell = 1, size(tetrahedra, 2)
            jacobian(:, 1) = vertices(:, tetrahedra(2, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 2) = vertices(:, tetrahedra(3, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 3) = vertices(:, tetrahedra(4, cell)) - &
                vertices(:, tetrahedra(1, cell))
            call project_physical_tetra_discontinuous( &
                basis, vertices(:, tetrahedra(:, cell)), exact_field, &
                coefficients, status)
            if (status /= 0) error stop "physical L2 projection failed"
            do node = 1, size(weights)
                reference_point = [x(node), y(node), z(node)]
                call evaluate_tetra_discontinuous( &
                    basis, reference_point(1), reference_point(2), &
                    reference_point(3), values, status)
                if (status /= 0) error stop "L2 basis evaluation failed"
                physical_point = vertices(:, tetrahedra(1, cell)) + &
                    matmul(jacobian, reference_point)
                projected_value = dot_product(values, coefficients)
                error = error + det3(jacobian)*weights(node)* &
                    (projected_value - exact_field( &
                    physical_point(1), physical_point(2), physical_point(3)))**2
            end do
        end do
        error = sqrt(error)
    end subroutine measure_error

    pure real(dp) function exact_field(x, y, z) result(value)
        real(dp), intent(in) :: x, y, z
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = sin(pi*x)*sin(pi*y)*sin(pi*z)
    end function exact_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_l2_h_convergence
