program test_maxwell_torus_curved_rbc_pairing
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_rwg_rbc_pairing, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgesv
    end interface

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    integer, allocatable :: pivots(:), triangles(:, :)
    real(dp), allocatable :: coarse(:, :), matrix(:, :), parameters(:, :)
    real(dp), allocatable :: right_hand_side(:, :), scaled(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :), work(:, :)
    real(dp) :: error, residual
    integer :: basis, info, status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 4, 6, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
        vertices, triangles, parameters, major_radius, minor_radius, 6, &
        coarse, status)
    call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
        vertices, triangles, parameters, major_radius, minor_radius, 12, &
        matrix, status)
    error = maxval(abs(matrix - coarse))/maxval(abs(matrix))
    call record_condition(status == 0 .and. error < 3.0e-5_dp, &
        "exact-torus RWG-RBC pairing converges under quadrature refinement")

    allocate( &
        work(size(matrix, 1), size(matrix, 2)), &
        right_hand_side(size(matrix, 1), 1), pivots(size(matrix, 1)))
    work = matrix
    do basis = 1, size(matrix, 1)
        right_hand_side(basis, 1) = sin(real(2*basis + 1, dp))
    end do
    call dgesv( &
        size(matrix, 1), 1, work, size(matrix, 1), pivots, right_hand_side, &
        size(matrix, 1), info)
    residual = norm2(matmul(matrix, right_hand_side(:, 1)) - &
        [(sin(real(2*basis + 1, dp)), basis=1, size(matrix, 1))])
    call record_condition(info == 0 .and. residual < 3.0e-12_dp, &
        "exact-torus RWG-RBC pairing is full rank")

    scaled_vertices = 3.0_dp*vertices
    call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
        scaled_vertices, triangles, parameters, 3.0_dp*major_radius, &
        3.0_dp*minor_radius, 12, scaled, status)
    error = maxval(abs(scaled - 3.0_dp*matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. error < 2.0e-12_dp, &
        "exact-torus RWG-RBC pairing obeys analytical BC length scaling")

    call check_summary("Exact-curved torus RWG-RBC duality")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_rbc_pairing
