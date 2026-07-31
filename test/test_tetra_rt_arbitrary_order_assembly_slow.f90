program test_tetra_rt_arbitrary_order_assembly_slow
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_tetra_rt_div_mass_csc, assemble_tetra_rt_div_mass_element, &
        assemble_tetra_rt_divergence_csc, build_tetra_rt_dof_map
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: expected_counts(0:6) = &
        [4, 15, 36, 70, 120, 189, 280]
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    integer, allocatable :: face_orientations(:, :)
    integer, allocatable :: face_permutations(:, :, :), faces(:, :)
    integer, allocatable :: global_dofs(:, :)
    real(dp), allocatable :: coefficients(:), element_matrix(:, :)
    real(dp), allocatable :: product(:)
    real(dp) :: field(3), mesh_vertices(3, 5), vertices(3, 4)
    integer :: degree, face, status, tetrahedra(4, 2)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.2_dp, -0.3_dp, 0.1_dp]
    vertices(:, 2) = [2.1_dp, -0.5_dp, 0.4_dp]
    vertices(:, 3) = [0.6_dp, 1.4_dp, -0.2_dp]
    vertices(:, 4) = [-0.1_dp, 0.2_dp, 1.6_dp]
    do degree = 0, 6
        call assemble_tetra_rt_div_mass_element( &
            vertices, degree, 2*degree + 3, element_matrix, status)
        call record_condition(status == 0 .and. &
            size(element_matrix, 1) == expected_counts(degree) .and. &
            maxval(abs(element_matrix - transpose(element_matrix))) < &
            2.0e-11_dp, &
            "Tetrahedral RT element assembly is symmetric at every degree")
    end do

    mesh_vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    mesh_vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    mesh_vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    mesh_vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    mesh_vertices(:, 5) = [0.0_dp, 0.0_dp, -1.0_dp]
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [1, 3, 2, 5]
    call build_tetra_rt_dof_map( &
        0, tetrahedra, faces, global_dofs, face_orientations, &
        face_permutations, status)
    allocate(coefficients(size(faces, 2)))
    field = [2.0_dp, -1.0_dp, 0.5_dp]
    do face = 1, size(faces, 2)
        coefficients(face) = 0.5_dp*dot_product( &
            field, cross_product( &
            mesh_vertices(:, faces(2, face)) - &
            mesh_vertices(:, faces(1, face)), &
            mesh_vertices(:, faces(3, face)) - &
            mesh_vertices(:, faces(1, face))))
    end do
    call assemble_tetra_rt_div_mass_csc( &
        mesh_vertices, tetrahedra, 0, 3, matrix, sparse_status, &
        divergence_coefficient=1.0_dp, mass_coefficient=1.0_dp)
    product = csc_matvec(matrix, coefficients)
    call record_condition(sparse_status%code == 0 .and. &
        abs(dot_product(coefficients, product) - &
        dot_product(field, field)/3.0_dp) < 2.0e-12_dp, &
        "Global RT assembly preserves constant-field energy across a face")

    call build_tetra_rt_dof_map( &
        1, tetrahedra, faces, global_dofs, face_orientations, &
        face_permutations, status)
    deallocate(coefficients)
    allocate(coefficients(maxval(global_dofs)))
    coefficients = 0.0_dp
    do face = 1, size(faces, 2)
        call linear_field_face_moments( &
            mesh_vertices(:, faces(:, face)), &
            coefficients(3*face - 2:3*face))
    end do
    coefficients(global_dofs(13, 1)) = 1.0_dp/24.0_dp
    coefficients(global_dofs(14, 2)) = 1.0_dp/24.0_dp
    call assemble_tetra_rt_divergence_csc( &
        mesh_vertices, tetrahedra, 1, 4, matrix, sparse_status)
    product = csc_matvec(matrix, coefficients)
    call record_condition(sparse_status%code == 0 .and. &
        maxval(abs(product - [ &
        1.0_dp/6.0_dp, 1.0_dp/24.0_dp, 1.0_dp/24.0_dp, 1.0_dp/24.0_dp, &
        1.0_dp/6.0_dp, 1.0_dp/24.0_dp, 1.0_dp/24.0_dp, 1.0_dp/24.0_dp])) < &
        3.0e-12_dp, &
        "RT-to-DG divergence commutes for the exact field v=(x,0,0)")

    call check_summary("Tetrahedral RT arbitrary-order assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine linear_field_face_moments(face_vertices, moments)
        real(dp), intent(in) :: face_vertices(3, 3)
        real(dp), intent(out) :: moments(3)

        integer, parameter :: quadrature_order = 5
        real(dp) :: area_normal(3), nodes(quadrature_order), point(3)
        real(dp) :: s, t, weights(quadrature_order)
        integer :: first, second

        area_normal = cross_product( &
            face_vertices(:, 2) - face_vertices(:, 1), &
            face_vertices(:, 3) - face_vertices(:, 1))
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)
        moments = 0.0_dp
        do first = 1, quadrature_order
            s = nodes(first)
            do second = 1, quadrature_order
                t = (1.0_dp - s)*nodes(second)
                point = face_vertices(:, 1) + &
                    s*(face_vertices(:, 2) - face_vertices(:, 1)) + &
                    t*(face_vertices(:, 3) - face_vertices(:, 1))
                moments = moments + weights(first)*weights(second)* &
                    (1.0_dp - s)*point(1)*area_normal(1)*[1.0_dp, t, s]
            end do
        end do
    end subroutine linear_field_face_moments

    pure function cross_product(first, second) result(product_)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product_(3)

        product_ = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_rt_arbitrary_order_assembly_slow
