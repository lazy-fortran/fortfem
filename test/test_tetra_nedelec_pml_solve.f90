program test_tetra_nedelec_pml_solve
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_tetra_nedelec_pml_csc, &
        solve_tetra_nedelec_pml
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_z_t, fortsparse_status_t
    implicit none

    type(csc_z_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: vertices(3, 5)
    integer :: tetrahedra(4, 2)
    complex(dp) :: stretch(3, 2)
    complex(dp), allocatable :: exact(:), load(:), solution(:)
    complex(dp) :: boundary_form(3, 3)
    complex(dp) :: prescribed_values(2)
    integer :: boundary_operator_dofs(3), prescribed_dofs(2), dof, status
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    vertices(:, 5) = [0.0_dp, 0.0_dp, -1.0_dp]
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [1, 3, 2, 5]
    stretch(:, 1) = [ &
        cmplx(1.0_dp, 0.4_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)]
    stretch(:, 2) = stretch(:, 1)
    call assemble_tetra_nedelec_pml_csc( &
        vertices, tetrahedra, 2, stretch, 1.1_dp, matrix, sparse_status)
    if (sparse_status%code /= 0) error stop "PML manufactured matrix failed"
    allocate(exact(matrix%nrow))
    do dof = 1, size(exact)
        exact(dof) = cmplx( &
            sin(0.3_dp*real(dof, dp)), cos(0.2_dp*real(dof, dp)), dp)
    end do
    load = csc_matvec(matrix, exact)
    boundary_operator_dofs = [2, 3, 4]
    boundary_form = reshape([ &
        cmplx(0.5_dp, -0.2_dp, dp), cmplx(0.1_dp, 0.05_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.1_dp, 0.05_dp, dp), &
        cmplx(0.7_dp, -0.1_dp, dp), cmplx(0.08_dp, -0.04_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.08_dp, -0.04_dp, dp), &
        cmplx(0.6_dp, -0.3_dp, dp)], [3, 3])
    load(boundary_operator_dofs) = load(boundary_operator_dofs) + &
        matmul(boundary_form, exact(boundary_operator_dofs))
    prescribed_dofs = [1, matrix%nrow]
    prescribed_values = exact(prescribed_dofs)
    call solve_tetra_nedelec_pml( &
        vertices, tetrahedra, 2, stretch, 1.1_dp, load, prescribed_dofs, &
        prescribed_values, solution, status, &
        boundary_operator_dofs=boundary_operator_dofs, &
        boundary_operator=boundary_form)
    call record_condition(status == 0, &
        "complex Maxwell PML boundary-value solve succeeds")
    call record_condition(maxval(abs(solution - exact)) < 2.0e-10_dp, &
        "Maxwell solve recovers a field with an active DtN boundary form")
    call check_summary("Tetrahedral Nedelec PML solve")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_pml_solve
