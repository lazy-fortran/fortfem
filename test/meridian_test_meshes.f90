module meridian_test_meshes
    !! Curved meridional test meshes: an annulus about (r0, 0) and a warped
    !! rectangle.  Both use polynomial isoparametric geometry of order geo.
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, &
        meridian_mesh_from_map
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), save :: ann_r0 = 3.0_dp

    public :: make_annulus_mesh, make_warped_rectangle_mesh
    public :: make_bubble_rectangle_mesh

contains

    subroutine annulus_map(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)

        x = [ann_r0 + p(1)*cos(p(2)), p(1)*sin(p(2))]
    end subroutine annulus_map

    subroutine warp_map(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)

        x(1) = p(1) + 0.12_dp*sin(2.0_dp*p(2)) + 0.05_dp*p(2)*p(1)
        x(2) = p(2) + 0.10_dp*sin(3.0_dp*p(1)) - 0.04_dp*p(1)*p(2)
    end subroutine warp_map

    subroutine bubble_map(p, x)
        !! Interior warp of [2, 3] x [-0.5, 0.5] that is the identity on the
        !! boundary, so the curved-element mesh covers the rectangle exactly.
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)
        real(dp) :: b

        b = (p(1) - 2.0_dp)*(3.0_dp - p(1))*(p(2) + 0.5_dp)*(0.5_dp - p(2))
        x(1) = p(1) + 0.9_dp*b*sin(3.0_dp*p(2) + p(1))
        x(2) = p(2) + 0.9_dp*b*cos(2.0_dp*p(1) - p(2))
    end subroutine bubble_map

    subroutine make_annulus_mesh(nr, nth, r1, r2, r0, geo, mm)
        !! Annulus r1 < |x - (r0, 0)| < r2 with nr radial and nth angular
        !! cells; r2 < r0.
        integer, intent(in) :: nr, nth, geo
        real(dp), intent(in) :: r1, r2, r0
        type(meridian_mesh_t), intent(out) :: mm
        type(mesh_2d_t) :: topo
        real(dp), allocatable :: pv(:, :, :)
        real(dp) :: dr, dth, x(2)
        integer :: i, j, jn, k, status

        ann_r0 = r0
        dr = (r2 - r1)/real(nr, dp)
        dth = 2.0_dp*pi/real(nth, dp)
        topo%n_vertices = (nr + 1)*nth
        topo%n_triangles = 2*nr*nth
        topo%has_triangles = .true.
        allocate (topo%vertices(2, topo%n_vertices))
        allocate (topo%triangles(3, topo%n_triangles))
        allocate (pv(2, 3, topo%n_triangles))
        do i = 0, nr
            do j = 0, nth - 1
                call annulus_map([r1 + i*dr, j*dth], x)
                topo%vertices(:, i*nth + j + 1) = x
            end do
        end do
        k = 0
        do i = 0, nr - 1
            do j = 0, nth - 1
                jn = mod(j + 1, nth)
                k = k + 1
                topo%triangles(:, k) = [i*nth + j + 1, (i + 1)*nth + j + 1, &
                    (i + 1)*nth + jn + 1]
                pv(:, :, k) = reshape([r1 + i*dr, j*dth, r1 + (i + 1)*dr, &
                    j*dth, r1 + (i + 1)*dr, (j + 1)*dth], [2, 3])
                k = k + 1
                topo%triangles(:, k) = [i*nth + j + 1, (i + 1)*nth + jn + 1, &
                    i*nth + jn + 1]
                pv(:, :, k) = reshape([r1 + i*dr, j*dth, r1 + (i + 1)*dr, &
                    (j + 1)*dth, r1 + i*dr, (j + 1)*dth], [2, 3])
            end do
        end do
        call meridian_mesh_from_map(topo, pv, geo, mm, status, annulus_map)
        if (status /= 0) error stop "annulus mesh failed"
    end subroutine make_annulus_mesh

    subroutine make_bubble_rectangle_mesh(nx, geo, mm)
        !! Exactly the rectangle [2, 3] x [-0.5, 0.5], curved interior
        !! elements of geometry order geo (see bubble_map).
        integer, intent(in) :: nx, geo
        type(meridian_mesh_t), intent(out) :: mm
        type(mesh_2d_t) :: topo, pmesh
        real(dp), allocatable :: pv(:, :, :)
        integer :: t, v, status

        call pmesh%create_rectangular(nx + 1, nx + 1, 2.0_dp, 3.0_dp, &
            -0.5_dp, 0.5_dp)
        topo = pmesh
        allocate (pv(2, 3, topo%n_triangles))
        do t = 1, topo%n_triangles
            pv(:, :, t) = pmesh%vertices(:, pmesh%triangles(:, t))
        end do
        do v = 1, topo%n_vertices
            call bubble_map(pmesh%vertices(:, v), topo%vertices(:, v))
        end do
        call meridian_mesh_from_map(topo, pv, geo, mm, status, bubble_map)
        if (status /= 0) error stop "bubble rectangle mesh failed"
    end subroutine make_bubble_rectangle_mesh

    subroutine make_warped_rectangle_mesh(nx, geo, mm)
        !! Warped image of [2, 3] x [-0.5, 0.5] with nx x nx cells.
        integer, intent(in) :: nx, geo
        type(meridian_mesh_t), intent(out) :: mm
        type(mesh_2d_t) :: topo, pmesh
        real(dp), allocatable :: pv(:, :, :)
        integer :: t, v, status

        call pmesh%create_rectangular(nx + 1, nx + 1, 2.0_dp, 3.0_dp, &
            -0.5_dp, 0.5_dp)
        topo = pmesh
        allocate (pv(2, 3, topo%n_triangles))
        do t = 1, topo%n_triangles
            pv(:, :, t) = pmesh%vertices(:, pmesh%triangles(:, t))
        end do
        do v = 1, topo%n_vertices
            call warp_map(pmesh%vertices(:, v), topo%vertices(:, v))
        end do
        call meridian_mesh_from_map(topo, pv, geo, mm, status, warp_map)
        if (status /= 0) error stop "warped mesh failed"
    end subroutine make_warped_rectangle_mesh

end module meridian_test_meshes
