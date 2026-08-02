program test_solid_torus_tetra_mesh
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_solid_torus_tetra_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    integer, parameter :: radial_layers = 2, poloidal_count = 8
    integer, parameter :: toroidal_count = 10
    integer, parameter :: face_nodes(3, 4) = reshape( &
        [2, 3, 4, 1, 4, 3, 1, 2, 4, 1, 3, 2], [3, 4])
    integer, allocatable :: boundary_triangles(:, :), tetrahedra(:, :)
    integer, allocatable :: faces(:, :)
    real(dp), allocatable :: boundary_parameters(:, :), vertices(:, :)
    real(dp) :: radial_distance, signed_volume, torus_residual
    integer :: boundary_face_count, face, face_count, local_face, tetrahedron
    integer :: vertex
    logical :: all_passed

    all_passed = .true.
    call generate_solid_torus_tetra_mesh( &
        major_radius, minor_radius, radial_layers, poloidal_count, &
        toroidal_count, vertices, tetrahedra, boundary_triangles, &
        boundary_parameters)

    call record_condition(size(boundary_triangles, 2) == &
        2*poloidal_count*toroidal_count, &
        "solid torus exposes the expected closed triangular boundary")
    call record_condition(size(boundary_parameters, 1) == 2 .and. &
        size(boundary_parameters, 2) == size(vertices, 2), &
        "solid torus supplies toroidal parameters for every volume vertex")
    do tetrahedron = 1, size(tetrahedra, 2)
        signed_volume = determinant3( &
            vertices(:, tetrahedra(2, tetrahedron)) - &
            vertices(:, tetrahedra(1, tetrahedron)), &
            vertices(:, tetrahedra(3, tetrahedron)) - &
            vertices(:, tetrahedra(1, tetrahedron)), &
            vertices(:, tetrahedra(4, tetrahedron)) - &
            vertices(:, tetrahedra(1, tetrahedron)))/6.0_dp
        if (signed_volume <= 1.0e-12_dp) all_passed = .false.
    end do
    call record_condition(all_passed, &
        "every structured solid-torus tetrahedron has positive volume")
    allocate(faces(3, 4*size(tetrahedra, 2)))
    face = 0
    do tetrahedron = 1, size(tetrahedra, 2)
        do local_face = 1, 4
            face = face + 1
            faces(:, face) = tetrahedra(face_nodes(:, local_face), tetrahedron)
            call sort_three(faces(:, face))
        end do
    end do
    boundary_face_count = 0
    do face = 1, size(faces, 2)
        face_count = count(all(faces == &
            spread(faces(:, face), 2, size(faces, 2)), dim=1))
        if (face_count == 1) boundary_face_count = boundary_face_count + 1
        if (face_count > 2) all_passed = .false.
    end do
    call record_condition(all_passed .and. boundary_face_count == &
        size(boundary_triangles, 2), &
        "tetrahedra meet face-to-face with exactly the declared outer boundary")
    torus_residual = 0.0_dp
    do vertex = 1, size(boundary_triangles, 2)
        radial_distance = sqrt(sum( &
            vertices(:2, boundary_triangles(1, vertex))**2))
        torus_residual = max(torus_residual, abs( &
            (radial_distance - major_radius)**2 + &
            vertices(3, boundary_triangles(1, vertex))**2 - minor_radius**2))
    end do
    call record_condition(torus_residual < 2.0e-14_dp, &
        "every boundary vertex lies on the exact implicit torus")

    call check_summary("Solid torus tetrahedral mesh")
    if (.not. all_passed) error stop 1

contains

    pure function determinant3(first, second, third) result(value)
        real(dp), intent(in) :: first(3), second(3), third(3)
        real(dp) :: value

        value = dot_product(first, cross_product(second, third))
    end function determinant3

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    pure subroutine sort_three(values)
        integer, intent(inout) :: values(3)

        integer :: first, second, temporary

        do first = 1, 2
            do second = first + 1, 3
                if (values(second) < values(first)) then
                    temporary = values(first)
                    values(first) = values(second)
                    values(second) = temporary
                end if
            end do
        end do
    end subroutine sort_three

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_solid_torus_tetra_mesh
