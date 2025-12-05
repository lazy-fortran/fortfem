program test_6pt_step
    use fortfem_kinds
    use delaunay_types
    use bowyer_watson
    use point_location, only: build_adjacency
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 6)
    integer :: i, super_tri_idx, point_idx, valid_count
    logical :: contains_super

    write(*,*) "=== 6-Point Step-by-Step Debug ==="

    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [2.0_dp, 0.0_dp]
    points(:, 3) = [2.0_dp, 1.0_dp]
    points(:, 4) = [1.0_dp, 1.0_dp]
    points(:, 5) = [1.0_dp, 2.0_dp]
    points(:, 6) = [0.0_dp, 2.0_dp]

    call create_mesh(mesh, 9, 30, 15)
    call create_super_triangle(points, mesh, super_tri_idx)

    write(*,*) "Super-triangle vertices:", mesh%super_vertices

    do i = 1, 6
        point_idx = add_point(mesh, points(1, i), points(2, i), i)
        write(*,*) "Added point", i, "as mesh point", point_idx
    end do

    call build_adjacency(mesh)

    write(*,*) ""
    write(*,*) "Before point insertion:"
    valid_count = 0
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            valid_count = valid_count + 1
            write(*,*) "  Triangle", i, "v=", mesh%triangles(i)%vertices,           &
                "neighbors=", mesh%triangles(i)%neighbors
        end if
    end do
    write(*,*) "Valid before insert:", valid_count

    ! Insert points one at a time
    do i = 4, mesh%npoints
        write(*,*) ""
        write(*,*) "Inserting point", i
        call insert_point(mesh, i)

        valid_count = 0
        do point_idx = 1, mesh%ntriangles
            if (mesh%triangles(point_idx)%valid) valid_count = valid_count + 1
        end do
        write(*,*) "  Valid triangles after:", valid_count
    end do

    write(*,*) ""
    write(*,*) "Before super-triangle removal:"
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            contains_super = .false.
            do point_idx = 1, 3
                if (any(mesh%triangles(i)%vertices == mesh%super_vertices(point_idx))) then
                    contains_super = .true.
                end if
            end do
            if (.not. contains_super) then
                write(*,*) "  Interior triangle", i, "v=", mesh%triangles(i)%vertices
            end if
        end if
    end do
end program test_6pt_step
