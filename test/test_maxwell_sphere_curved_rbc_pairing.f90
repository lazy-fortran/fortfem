program test_maxwell_sphere_curved_rbc_pairing
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing, &
        generate_sphere_surface_mesh
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

    integer, allocatable :: pivots(:), triangles(:, :)
    real(dp), allocatable :: coarse(:, :), matrix(:, :), right_hand_side(:, :)
    real(dp), allocatable :: original_right_hand_side(:), residual_vector(:)
    real(dp), allocatable :: scaled(:, :), scaled_vertices(:, :), vertices(:, :)
    real(dp), allocatable :: work(:, :)
    real(dp) :: error, residual
    integer :: basis, info, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
        vertices, triangles, 1.0_dp, 6, coarse, status)
    call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
        vertices, triangles, 1.0_dp, 12, matrix, status)
    error = maxval(abs(matrix - coarse))/maxval(abs(matrix))
    call record_condition(status == 0 .and. error < 2.0e-5_dp, &
        "curved RWG-RBC pairing converges under quadrature refinement")

    allocate( &
        work(size(matrix, 1), size(matrix, 2)), &
        right_hand_side(size(matrix, 1), 1), pivots(size(matrix, 1)), &
        original_right_hand_side(size(matrix, 1)), &
        residual_vector(size(matrix, 1)))
    work = matrix
    do basis = 1, size(matrix, 1)
        right_hand_side(basis, 1) = sin(real(2*basis + 1, dp))
    end do
    original_right_hand_side = right_hand_side(:, 1)
    call dgesv( &
        size(matrix, 1), 1, work, size(matrix, 1), pivots, right_hand_side, &
        size(matrix, 1), info)
    residual_vector = matmul(matrix, right_hand_side(:, 1)) - &
        original_right_hand_side
    residual = sqrt(sum(residual_vector**2))
    call record_condition(info == 0 .and. residual < 2.0e-12_dp, &
        "curved RWG-RBC pairing is full rank")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
        scaled_vertices, triangles, 2.0_dp, 12, scaled, status)
    error = maxval(abs(scaled - 2.0_dp*matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. error < 8.0e-13_dp, &
        "curved RWG-RBC pairing obeys analytical BC length scaling")

    call check_summary("Curved-sphere RWG-RBC duality")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_rbc_pairing
