program test_bspline_3d_sparse_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_bspline_h1_operator_3d_csc, &
        assemble_bspline_hcurl_operator_3d_csc, &
        assemble_bspline_hdiv_operator_3d_csc, &
        assemble_bspline_hdiv_l2_divergence_3d_csc, &
        assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc, &
        assemble_bspline_l2_mass_3d_csc, &
        build_bspline_feec_3d_operators_csc
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: knots_x(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.35_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_y(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    type(csc_t) :: adjoint_divergence, curl, divergence, gradient, hcurl_mass
    type(csc_t) :: hcurl_stiffness, hdiv_mass, hdiv_stiffness
    type(csc_t) :: l2_mass, mass, stiffness, weak_divergence
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: coefficients(:), controls(:, :, :, :)
    real(dp), allocatable :: hcurl_coefficients(:), hdiv_coefficients(:)
    real(dp), allocatable :: l2_coefficients(:)
    real(dp), allocatable :: product(:), weights(:, :, :)
    real(dp), allocatable :: x_points(:), y_points(:), z_points(:)
    real(dp) :: energy
    integer :: bx_count, by_count, ix, iy, iz, nx, ny, nz

    x_points = greville_abscissae(knots_x, 2)
    y_points = greville_abscissae(knots_y, 2)
    z_points = greville_abscissae(knots_y, 2)
    allocate( &
        controls(3, size(x_points), size(y_points), size(z_points)), &
        weights(size(x_points), size(y_points), size(z_points)))
    weights = 1.0_dp
    do iz = 1, size(z_points)
        do iy = 1, size(y_points)
            do ix = 1, size(x_points)
                controls(:, ix, iy, iz) = [ &
                    1.0_dp + 2.0_dp*x_points(ix), &
                    -1.0_dp + 3.0_dp*y_points(iy), &
                    0.5_dp + 4.0_dp*z_points(iz)]
            end do
        end do
    end do

    call assemble_bspline_h1_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, mass, &
        sparse_status, stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    allocate(coefficients(size(x_points)*size(y_points)*size(z_points)))
    coefficients = 1.0_dp
    product = csc_matvec(mass, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 24.0_dp) < 2.0e-11_dp, &
        "Sparse 3D spline mass integrates affine-box volume exactly")

    call assemble_bspline_h1_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, stiffness, &
        sparse_status, stiffness_coefficient=1.0_dp, mass_coefficient=0.0_dp)
    do iz = 1, size(z_points)
        do iy = 1, size(y_points)
            do ix = 1, size(x_points)
                coefficients(index_3d( &
                    ix, iy, iz, size(x_points), size(y_points))) = &
                    1.0_dp + 2.0_dp*x_points(ix)
            end do
        end do
    end do
    product = csc_matvec(stiffness, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(abs(energy - 24.0_dp) < 3.0e-10_dp, &
        "Sparse 3D spline diffusion gives exact coordinate-field energy")

    call assemble_bspline_l2_mass_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, l2_mass, &
        sparse_status)
    nx = size(x_points)
    ny = size(y_points)
    nz = size(z_points)
    allocate(l2_coefficients((nx - 1)*(ny - 1)*(nz - 1)))
    l2_coefficients = 24.0_dp
    product = csc_matvec(l2_mass, l2_coefficients)
    energy = dot_product(l2_coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 24.0_dp) < 3.0e-10_dp, &
        "Sparse 3D L2 mass preserves the physical constant field")

    call build_bspline_feec_3d_operators_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, gradient, curl, divergence, &
        sparse_status)
    bx_count = nx*(ny - 1)*(nz - 1)
    by_count = (nx - 1)*ny*(nz - 1)
    allocate(hdiv_coefficients( &
        bx_count + by_count + (nx - 1)*(ny - 1)*nz))
    do iz = 1, nz - 1
        do iy = 1, ny - 1
            do ix = 1, nx
                hdiv_coefficients(index_3d(ix, iy, iz, nx, ny - 1)) = &
                    12.0_dp*(1.0_dp + 2.0_dp*x_points(ix))
            end do
        end do
    end do
    do iz = 1, nz - 1
        do iy = 1, ny
            do ix = 1, nx - 1
                hdiv_coefficients(bx_count + &
                    index_3d(ix, iy, iz, nx - 1, ny)) = &
                    8.0_dp*(-1.0_dp + 3.0_dp*y_points(iy))
            end do
        end do
    end do
    do iz = 1, nz
        do iy = 1, ny - 1
            do ix = 1, nx - 1
                hdiv_coefficients(bx_count + by_count + &
                    index_3d(ix, iy, iz, nx - 1, ny - 1)) = &
                    6.0_dp*(0.5_dp + 4.0_dp*z_points(iz))
            end do
        end do
    end do
    product = csc_matvec(divergence, hdiv_coefficients)
    call check_condition(maxval(abs(product - 72.0_dp)) < 2.0e-12_dp, &
        "Sparse 3D divergence reproduces div(x,y,z) under Piola mapping")

    call assemble_bspline_hdiv_l2_divergence_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, &
        weak_divergence, sparse_status)
    product = csc_matvec(weak_divergence, hdiv_coefficients)
    energy = dot_product(csc_matvec(divergence, hdiv_coefficients), product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 216.0_dp) < 3.0e-9_dp, &
        "Weak 3D divergence has exact affine-field physical energy")

    call assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, &
        adjoint_divergence, sparse_status)
    product = csc_matvec(adjoint_divergence, l2_coefficients)
    energy = dot_product(hdiv_coefficients, product)
    call check_condition(abs(energy - 72.0_dp) < 3.0e-9_dp, &
        "Adjoint 3D divergence satisfies the physical L2 duality pairing")

    hdiv_coefficients = 0.0_dp
    hdiv_coefficients(:bx_count) = 12.0_dp
    call assemble_bspline_hdiv_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, hdiv_mass, &
        sparse_status, divergence_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    product = csc_matvec(hdiv_mass, hdiv_coefficients)
    energy = dot_product(hdiv_coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 24.0_dp) < 3.0e-9_dp, &
        "Sparse 3D Hdiv mass preserves a constant physical vector")

    call set_affine_hdiv_coefficients( &
        hdiv_coefficients, x_points, y_points, z_points)
    call assemble_bspline_hdiv_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, &
        hdiv_stiffness, sparse_status, divergence_coefficient=1.0_dp, &
        mass_coefficient=0.0_dp)
    product = csc_matvec(hdiv_stiffness, hdiv_coefficients)
    energy = dot_product(hdiv_coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 216.0_dp) < 3.0e-9_dp, &
        "Sparse 3D div-div form has exact affine-field energy")

    allocate(hcurl_coefficients( &
        (nx - 1)*ny*nz + nx*(ny - 1)*nz + nx*ny*(nz - 1)))
    hcurl_coefficients = 0.0_dp
    hcurl_coefficients(:(nx - 1)*ny*nz) = 2.0_dp
    call assemble_bspline_hcurl_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, hcurl_mass, &
        sparse_status, curl_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    product = csc_matvec(hcurl_mass, hcurl_coefficients)
    energy = dot_product(hcurl_coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 24.0_dp) < 3.0e-9_dp, &
        "Sparse 3D Hcurl mass preserves a constant physical vector")

    call set_rigid_rotation_hcurl_coefficients( &
        hcurl_coefficients, x_points, y_points, z_points)
    call assemble_bspline_hcurl_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, &
        hcurl_stiffness, sparse_status, curl_coefficient=1.0_dp, &
        mass_coefficient=0.0_dp)
    product = csc_matvec(hcurl_stiffness, hcurl_coefficients)
    energy = dot_product(hcurl_coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 24.0_dp) < 3.0e-9_dp, &
        "Sparse 3D curl-curl form has exact rigid-rotation energy")

    call check_summary("Sparse 3D isogeometric de Rham assembly")

contains

    pure function greville_abscissae(knots, degree) result(points)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), allocatable :: points(:)
        integer :: basis

        allocate(points(size(knots) - degree - 1))
        do basis = 1, size(points)
            points(basis) = sum(knots(basis + 1:basis + degree))/real(degree, dp)
        end do
    end function greville_abscissae

    pure integer function index_3d( &
            ix, iy, iz, count_x, count_y) result(index)
        integer, intent(in) :: ix, iy, iz, count_x, count_y

        index = ix + (iy - 1)*count_x + (iz - 1)*count_x*count_y
    end function index_3d

    subroutine set_affine_hdiv_coefficients(values, x, y, z)
        real(dp), intent(out) :: values(:)
        real(dp), intent(in) :: x(:), y(:), z(:)
        integer :: bx_size, by_size, i, j, k

        bx_size = size(x)*(size(y) - 1)*(size(z) - 1)
        by_size = (size(x) - 1)*size(y)*(size(z) - 1)
        do k = 1, size(z) - 1
            do j = 1, size(y) - 1
                do i = 1, size(x)
                    values(index_3d(i, j, k, size(x), size(y) - 1)) = &
                        12.0_dp*(1.0_dp + 2.0_dp*x(i))
                end do
            end do
        end do
        do k = 1, size(z) - 1
            do j = 1, size(y)
                do i = 1, size(x) - 1
                    values(bx_size + &
                        index_3d(i, j, k, size(x) - 1, size(y))) = &
                        8.0_dp*(-1.0_dp + 3.0_dp*y(j))
                end do
            end do
        end do
        do k = 1, size(z)
            do j = 1, size(y) - 1
                do i = 1, size(x) - 1
                    values(bx_size + by_size + &
                        index_3d(i, j, k, size(x) - 1, size(y) - 1)) = &
                        6.0_dp*(0.5_dp + 4.0_dp*z(k))
                end do
            end do
        end do
    end subroutine set_affine_hdiv_coefficients

    subroutine set_rigid_rotation_hcurl_coefficients(values, x, y, z)
        real(dp), intent(out) :: values(:)
        real(dp), intent(in) :: x(:), y(:), z(:)
        integer :: ex_size, i, j, k

        values = 0.0_dp
        ex_size = (size(x) - 1)*size(y)*size(z)
        do k = 1, size(z)
            do j = 1, size(y)
                do i = 1, size(x) - 1
                    values(index_3d( &
                        i, j, k, size(x) - 1, size(y))) = &
                        1.0_dp - 3.0_dp*y(j)
                end do
            end do
        end do
        do k = 1, size(z)
            do j = 1, size(y) - 1
                do i = 1, size(x)
                    values(ex_size + &
                        index_3d(i, j, k, size(x), size(y) - 1)) = &
                        1.5_dp*(1.0_dp + 2.0_dp*x(i))
                end do
            end do
        end do
    end subroutine set_rigid_rotation_hcurl_coefficients

end program test_bspline_3d_sparse_assembly
