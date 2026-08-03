program test_tetra_lagrange_poisson_p_convergence
    use check, only: check_condition, check_summary
    use fortfem_feec, only: evaluate_tetra_lagrange_solution, &
        solve_tetra_lagrange_diffusion_reaction, &
        solve_tetra_lagrange_poisson, tetra_duffy_quadrature
    use fortfem_generated_tetra_h1_oracle, only: generated_tetra_h1_oracle
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: tetrahedra(4, 4) = reshape( &
        [1, 2, 3, 5, 1, 4, 2, 5, 1, 3, 4, 5, 2, 4, 3, 5], [4, 4])
    real(dp), parameter :: vertices(3, 5) = reshape( &
        [0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        0.25_dp, 0.25_dp, 0.25_dp], [3, 5])

    real(dp) :: h1_errors(4), l2_errors(4)
    logical :: all_passed
    integer :: degree

    all_passed = .true.
    do degree = 1, 4
        call solve_and_measure(degree, l2_errors(degree), h1_errors(degree))
    end do
    call record_condition( &
        all(l2_errors(2:) < l2_errors(:3)), &
        "tetrahedral H1 L2 error decreases under p refinement")
    call record_condition( &
        all(h1_errors(2:) < h1_errors(:3)), &
        "tetrahedral H1 gradient error decreases under p refinement")
    call record_condition( &
        l2_errors(4) < 1.0e-12_dp, &
        "degree-four solve reproduces the quartic field")
    call record_condition( &
        h1_errors(4) < 1.0e-11_dp, &
        "degree-four solve reproduces the quartic gradient")
    call check_affine_nonzero_boundary()
    call check_reaction_without_essential_boundary()

    call check_summary("tetrahedral H1 Poisson p convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine solve_and_measure(degree, l2_error, h1_error)
        integer, intent(in) :: degree
        real(dp), intent(out) :: l2_error, h1_error

        type(fortsparse_status_t) :: status
        real(dp), allocatable :: solution(:)

        call solve_tetra_lagrange_poisson( &
            vertices, tetrahedra, degree, bubble_source, zero_boundary, &
            solution, status)
        call record_condition( &
            status%code == 0, "public tetrahedral H1 Poisson solve succeeds")
        if (status%code /= 0) then
            l2_error = huge(1.0_dp)
            h1_error = huge(1.0_dp)
            return
        end if
        call measure_error(degree, solution, l2_error, h1_error)
    end subroutine solve_and_measure

    subroutine measure_error(degree, solution, l2_error, h1_error)
        integer, intent(in) :: degree
        real(dp), intent(in) :: solution(:)
        real(dp), intent(out) :: l2_error, h1_error

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: exact(5), gradient(3), jacobian(3, 3), point(3)
        real(dp) :: value, weight
        integer :: cell, node, quadrature_status, sample, status

        call tetra_duffy_quadrature(10, x, y, z, weights, quadrature_status)
        call record_condition( &
            quadrature_status == 0, "error quadrature is available")
        l2_error = 0.0_dp
        h1_error = 0.0_dp
        do cell = 1, size(tetrahedra, 2)
            do node = 1, 3
                jacobian(:, node) = &
                    vertices(:, tetrahedra(node + 1, cell)) - &
                    vertices(:, tetrahedra(1, cell))
            end do
            do sample = 1, size(weights)
                point = vertices(:, tetrahedra(1, cell)) + &
                    matmul(jacobian, [x(sample), y(sample), z(sample)])
                call generated_tetra_h1_oracle( &
                    point(1), point(2), point(3), exact)
                call evaluate_tetra_lagrange_solution( &
                    vertices, tetrahedra, degree, solution, cell, &
                    [x(sample), y(sample), z(sample)], value, gradient, status)
                call record_condition( &
                    status == 0, "tetrahedral H1 reconstruction succeeds")
                weight = det3(jacobian)*weights(sample)
                l2_error = l2_error + weight*(value - exact(1))**2
                h1_error = h1_error + weight*sum((gradient - exact(2:4))**2)
            end do
        end do
        l2_error = sqrt(l2_error)
        h1_error = sqrt(h1_error)
    end subroutine measure_error

    subroutine check_affine_nonzero_boundary()
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: solution(:)
        real(dp) :: gradient(3), value
        integer :: evaluation_status

        call solve_tetra_lagrange_poisson( &
            vertices, tetrahedra, 1, zero_source, affine_boundary, solution, &
            status)
        call record_condition( &
            status%code == 0, "nonzero H1 Dirichlet solve succeeds")
        if (status%code /= 0) return
        call evaluate_tetra_lagrange_solution( &
            vertices, tetrahedra, 1, solution, 1, &
            [0.2_dp, 0.3_dp, 0.1_dp], value, gradient, evaluation_status)
        call record_condition( &
            evaluation_status == 0, "affine H1 reconstruction succeeds")
        call record_condition( &
            abs(value - affine_value(0.225_dp, 0.325_dp, 0.025_dp)) < &
            1.0e-12_dp, &
            "nonzero Dirichlet solve reproduces an affine harmonic field")
        call record_condition( &
            maxval(abs(gradient - [1.0_dp, 2.0_dp, -3.0_dp])) < 1.0e-12_dp, &
            "nonzero Dirichlet solve reproduces the affine gradient")
    end subroutine check_affine_nonzero_boundary

    subroutine check_reaction_without_essential_boundary()
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: solution(:)
        real(dp) :: gradient(3), value
        integer :: evaluation_status

        call solve_tetra_lagrange_diffusion_reaction( &
            vertices, tetrahedra, 1, constant_source, 0.0_dp, 1.0_dp, &
            solution, status)
        call record_condition( &
            status%code == 0, "unconstrained reaction solve succeeds")
        if (status%code /= 0) return
        call evaluate_tetra_lagrange_solution( &
            vertices, tetrahedra, 1, solution, 4, &
            [0.2_dp, 0.1_dp, 0.3_dp], value, gradient, evaluation_status)
        call record_condition( &
            evaluation_status == 0, "reaction field reconstruction succeeds")
        call record_condition( &
            abs(value - 2.0_dp) < 1.0e-12_dp, &
            "reaction solve reproduces a constant field")
        call record_condition( &
            maxval(abs(gradient)) < 1.0e-12_dp, &
            "reaction solve reproduces a zero gradient")
    end subroutine check_reaction_without_essential_boundary

    pure subroutine bubble_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        real(dp) :: oracle(5)

        call generated_tetra_h1_oracle(x, y, z, oracle)
        value = oracle(5)
    end subroutine bubble_source

    pure subroutine zero_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        value = 0.0_dp*(x + y + z)
    end subroutine zero_source

    pure subroutine constant_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        value = 2.0_dp + 0.0_dp*(x + y + z)
    end subroutine constant_source

    pure subroutine zero_boundary(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        value = 0.0_dp*(x + y + z)
    end subroutine zero_boundary

    pure subroutine affine_boundary(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        value = affine_value(x, y, z)
    end subroutine affine_boundary

    pure real(dp) function affine_value(x, y, z) result(value)
        real(dp), intent(in) :: x, y, z

        value = 1.0_dp + x + 2.0_dp*y - 3.0_dp*z
    end function affine_value

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_poisson_p_convergence
