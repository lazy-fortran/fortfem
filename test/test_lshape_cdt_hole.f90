program test_lshape_cdt_hole
    use fortfem_kinds, only: dp
    use fortfem_api, only: boundary_t, mesh_t, l_shape_boundary,              &
        mesh_from_boundary
    use check, only: check_condition, check_summary
    implicit none

    type(boundary_t) :: boundary
    type(mesh_t) :: mesh
    real(dp) :: size, resolution
    integer :: t, v1, v2, v3
    real(dp) :: cx, cy
    integer :: n_bad
    real(dp), parameter :: eps = 1.0e-8_dp

    size = 1.0_dp
    resolution = 0.08_dp

    boundary = l_shape_boundary(size, 24)
    mesh = mesh_from_boundary(boundary, resolution=resolution)

    call check_condition(mesh%data%n_triangles > 0, &
        "L-shape CDT: mesh has triangles")

    ! Reentrant cavity occupies rectangle (x > size, y < size).
    n_bad = 0
    do t = 1, mesh%data%n_triangles
        v1 = mesh%data%triangles(1, t)
        v2 = mesh%data%triangles(2, t)
        v3 = mesh%data%triangles(3, t)

        cx = (mesh%data%vertices(1, v1) + mesh%data%vertices(1, v2) +          &
              mesh%data%vertices(1, v3)) / 3.0_dp
        cy = (mesh%data%vertices(2, v1) + mesh%data%vertices(2, v2) +          &
              mesh%data%vertices(2, v3)) / 3.0_dp

        if (cx > size + eps .and. cy < size - eps) then
            n_bad = n_bad + 1
        end if
    end do

    call check_condition(n_bad == 0, &
        "L-shape CDT: reentrant square remains empty")
    call check_summary("L-shape CDT hole removal")
end program test_lshape_cdt_hole
