program test_scalar_helmholtz_outgoing_waves
    !! Time convention exp(-i omega t): a plane wave exp(i k x) driven by u = 1 at
    !! x = 0 travels to +x and must leave through the DtN boundary without
    !! reflection, and must be damped by a PML with stretch s = 1 + i sigma/k so
    !! that the wave reflected at its outer Dirichlet wall has the predicted
    !! amplitude exp(-2 int sigma dx) = exp(-2 sigma_max L/3) (sigma quadratic over
    !! the layer width L). The opposite sign of either boundary would reflect
    !! fully or amplify (|R| = exp(+2 sigma_max L/3)). Lateral faces are natural,
    !! so the exact solution is one-dimensional; R is fitted from u = exp(i k x) +
    !! R exp(-i k x) (up to normalization) at two nodes of the physical region.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_scalar_helmholtz_planar_dtn_p1, &
        solve_scalar_helmholtz_pml_p1_2d, solve_scalar_helmholtz_pml_p1_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: sigma_max = 4.0_dp
    real(dp) :: dtn_reflection(2)
    logical :: all_passed

    all_passed = .true.
    call dtn_case(33, dtn_reflection(1))
    call dtn_case(65, dtn_reflection(2))
    ! The DtN map samples a periodic boundary (period n h) while the strip is
    ! (n - 1) h high, so the end samples carry an O(h) excess weight.
    call record(dtn_reflection(2) < 1.0e-2_dp .and. &
        dtn_reflection(2) < 0.6_dp*dtn_reflection(1), &
        "DtN boundary lets exp(i k x) leave; reflection O(h) -> 0")
    call pml_2d_case()
    call pml_3d_case()
    call check_summary("Outgoing scalar Helmholtz waves under exp(-i omega t)")
    if (.not. all_passed) error stop 1

contains

    subroutine dtn_case(n, reflection)
        integer, intent(in) :: n
        real(dp), intent(out) :: reflection
        real(dp), parameter :: k = 2.0_dp
        real(dp), allocatable :: vertices(:, :)
        integer, allocatable :: triangles(:, :), dtn(:), left(:)
        complex(dp), allocatable :: load(:), solution(:)
        complex(dp) :: transmitted, reflected
        integer :: row, status, middle

        call strip(n, n, 1.0_dp, vertices, triangles)
        dtn = [(n*row, row=1, n)]
        left = [(1 + n*(row - 1), row=1, n)]
        allocate(load(n*n), solution(n*n))
        load = (0.0_dp, 0.0_dp)
        call solve_scalar_helmholtz_planar_dtn_p1(vertices, triangles, dtn, k, &
            real(n, dp)/real(n - 1, dp), load, left, &
            spread((1.0_dp, 0.0_dp), 1, n), solution, status)
        middle = 1 + n*(n/2)
        call fit_waves(k, vertices(1, middle + n/2), solution(middle), &
            solution(middle + n/2), transmitted, reflected)
        reflection = abs(reflected/transmitted)
        write (*, "(a,i0,a,es10.3)") "DtN reflection, n = ", n, ": ", reflection
        call record(status == 0, "DtN solve succeeds")
    end subroutine dtn_case

    subroutine pml_2d_case()
        integer, parameter :: nx = 81, ny = 5
        real(dp), parameter :: k = 3.0_dp
        real(dp), allocatable :: vertices(:, :)
        integer, allocatable :: triangles(:, :), walls(:)
        complex(dp), allocatable :: stretch(:, :), load(:), solution(:)
        complex(dp) :: transmitted, reflected
        integer :: e, row, status

        call strip(nx, ny, 2.0_dp, vertices, triangles)
        allocate(stretch(2, size(triangles, 2)))
        stretch = (1.0_dp, 0.0_dp)
        do e = 1, size(triangles, 2)
            stretch(1, e) = layer_stretch(sum(vertices(1, triangles(:, e)))/3.0_dp, &
                k)
        end do
        walls = [[(1 + nx*(row - 1), row=1, ny)], [(nx*row, row=1, ny)]]
        allocate(load(nx*ny))
        load = (0.0_dp, 0.0_dp)
        call solve_scalar_helmholtz_pml_p1_2d(vertices, triangles, stretch, k, &
            load, walls, [spread((1.0_dp, 0.0_dp), 1, ny), &
            spread((0.0_dp, 0.0_dp), 1, ny)], solution, status)
        call fit_waves(k, vertices(1, 11), solution(1), solution(11), &
            transmitted, reflected)
        write (*, "(a,2es10.3)") "2D PML reflection, predicted: ", &
            abs(reflected/transmitted), expected_reflection()
        call record(status == 0, "2D PML solve succeeds")
        call record(abs(log(abs(reflected/transmitted)/expected_reflection())) &
            < 8.0e-2_dp, "2D PML damps exp(i k x) as predicted")
    end subroutine pml_2d_case

    subroutine pml_3d_case()
        integer, parameter :: nx = 81, ny = 2
        real(dp), parameter :: k = 3.0_dp
        real(dp), allocatable :: vertices(:, :)
        integer, allocatable :: tetrahedra(:, :), walls(:)
        complex(dp), allocatable :: stretch(:, :), load(:), solution(:), values(:)
        complex(dp) :: transmitted, reflected
        integer :: e, status, v

        call slab_box(nx, ny, vertices, tetrahedra)
        allocate(stretch(3, size(tetrahedra, 2)))
        stretch = (1.0_dp, 0.0_dp)
        do e = 1, size(tetrahedra, 2)
            stretch(1, e) = layer_stretch(sum(vertices(1, tetrahedra(:, e)))/ &
                4.0_dp, k)
        end do
        ! The walls x = 0 and x = 2.
        walls = pack([(v, v=1, size(vertices, 2))], &
            abs(vertices(1, :) - 1.0_dp) > 1.0_dp - 1.0e-12_dp)
        values = merge((1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
            vertices(1, walls) < 1.0_dp)
        allocate(load(size(vertices, 2)))
        load = (0.0_dp, 0.0_dp)
        call solve_scalar_helmholtz_pml_p1_3d(vertices, tetrahedra, stretch, k, &
            load, walls, values, solution, status)
        ! Vertex 1 + 10 (ny + 1)**2 lies at x = 0.25 on the edge y = z = 0.
        v = 1 + 10*(ny + 1)**2
        call fit_waves(k, vertices(1, v), solution(1), solution(v), transmitted, &
            reflected)
        write (*, "(a,2es10.3)") "3D PML reflection, predicted: ", &
            abs(reflected/transmitted), expected_reflection()
        call record(status == 0, "3D PML solve succeeds")
        call record(abs(log(abs(reflected/transmitted)/expected_reflection())) &
            < 8.0e-2_dp, "3D PML damps exp(i k x) as predicted")
    end subroutine pml_3d_case

    pure complex(dp) function layer_stretch(x, k) result(s)
        !! s = 1 + i sigma/k in the layer 1 < x < 2, sigma = sigma_max (x - 1)**2.
        real(dp), intent(in) :: x, k

        s = (1.0_dp, 0.0_dp)
        if (x > 1.0_dp) s = cmplx(1.0_dp, sigma_max*(x - 1.0_dp)**2/k, dp)
    end function layer_stretch

    pure real(dp) function expected_reflection()
        expected_reflection = exp(-2.0_dp*sigma_max/3.0_dp)
    end function expected_reflection

    pure subroutine fit_waves(k, x, u0, u, transmitted, reflected)
        !! u(0) = T + R and u(x) = T exp(i k x) + R exp(-i k x).
        real(dp), intent(in) :: k, x
        complex(dp), intent(in) :: u0, u
        complex(dp), intent(out) :: transmitted, reflected
        complex(dp) :: forward, backward

        forward = exp(cmplx(0.0_dp, k*x, dp))
        backward = exp(cmplx(0.0_dp, -k*x, dp))
        reflected = (u - u0*forward)/(backward - forward)
        transmitted = u0 - reflected
    end subroutine fit_waves

    subroutine strip(nx, ny, width, vertices, triangles)
        !! [0, width] x [0, 1], vertex 1 + i + nx*j at (i, j) of the grid.
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: width
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: triangles(:, :)
        integer :: c, i, j, t

        allocate(vertices(2, nx*ny), triangles(3, 2*(nx - 1)*(ny - 1)))
        do j = 0, ny - 1
            do i = 0, nx - 1
                vertices(:, 1 + i + nx*j) = [width*real(i, dp)/real(nx - 1, dp), &
                    real(j, dp)/real(ny - 1, dp)]
            end do
        end do
        t = 0
        do j = 0, ny - 2
            do i = 1, nx - 1
                c = i + nx*j
                triangles(:, t + 1) = [c, c + 1, c + nx + 1]
                triangles(:, t + 2) = [c, c + nx + 1, c + nx]
                t = t + 2
            end do
        end do
    end subroutine strip

    subroutine slab_box(nx, ny, vertices, tetrahedra)
        !! [0, 2] x [0, 0.25]**2 with nx x (ny + 1) x (ny + 1) vertices, x slowest:
        !! vertex 1 + l + (ny + 1)*(j + (ny + 1)*i) at grid point (i, j, l).
        integer, intent(in) :: nx, ny
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: tetrahedra(:, :)
        integer, parameter :: kuhn(4, 6) = reshape([1, 2, 4, 8, 1, 4, 3, 8, &
            1, 3, 7, 8, 1, 7, 5, 8, 1, 5, 6, 8, 1, 6, 2, 8], [4, 6])
        integer :: corner(8), i, j, l, m, t, swap
        real(dp) :: e(3, 3)

        m = ny + 1
        allocate(vertices(3, nx*m*m), tetrahedra(4, 6*(nx - 1)*ny*ny))
        do i = 0, nx - 1
            do j = 0, ny
                do l = 0, ny
                    vertices(:, 1 + l + m*(j + m*i)) = [2.0_dp*real(i, dp)/ &
                        real(nx - 1, dp), 0.25_dp*real(j, dp)/real(ny, dp), &
                        0.25_dp*real(l, dp)/real(ny, dp)]
                end do
            end do
        end do
        t = 0
        do i = 0, nx - 2
            do j = 0, ny - 1
                do l = 0, ny - 1
                    corner = [node(i, j, l, m), node(i + 1, j, l, m), &
                        node(i, j + 1, l, m), node(i + 1, j + 1, l, m), &
                        node(i, j, l + 1, m), node(i + 1, j, l + 1, m), &
                        node(i, j + 1, l + 1, m), node(i + 1, j + 1, l + 1, m)]
                    do swap = 1, 6
                        t = t + 1
                        tetrahedra(:, t) = corner(kuhn(:, swap))
                        e(:, 1) = vertices(:, tetrahedra(2, t)) - &
                            vertices(:, tetrahedra(1, t))
                        e(:, 2) = vertices(:, tetrahedra(3, t)) - &
                            vertices(:, tetrahedra(1, t))
                        e(:, 3) = vertices(:, tetrahedra(4, t)) - &
                            vertices(:, tetrahedra(1, t))
                        if (dot_product(e(:, 1), [e(2, 2)*e(3, 3) - e(3, 2)* &
                            e(2, 3), e(3, 2)*e(1, 3) - e(1, 2)*e(3, 3), &
                            e(1, 2)*e(2, 3) - e(2, 2)*e(1, 3)]) < 0.0_dp) &
                            tetrahedra(3:4, t) = tetrahedra([4, 3], t)
                    end do
                end do
            end do
        end do
    end subroutine slab_box

    pure integer function node(a, b, c, m)
        integer, intent(in) :: a, b, c, m

        node = 1 + c + m*(b + m*a)
    end function node

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_scalar_helmholtz_outgoing_waves
