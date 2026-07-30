program test_bspline_h1_sparse_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_bspline_h1_operator_csc, &
        assemble_bspline_hcurl_operator_csc, &
        assemble_bspline_hdiv_operator_csc, &
        assemble_bspline_h1_hcurl_gradient_csc, &
        assemble_bspline_hcurl_h1_adjoint_gradient_csc, &
        assemble_bspline_hcurl_l2_curl_csc, &
        assemble_bspline_l2_hcurl_adjoint_curl_csc, &
        assemble_bspline_l2_mass_csc, &
        assemble_bspline_grad_shafranov_csc, &
        assemble_bspline_toroidal_fourier_laplacian_csc, &
        assemble_bspline_poloidal_bracket_csc, &
        build_bspline_feec_2d_operators_csc
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: knots_x(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.35_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_y(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    type(csc_t) :: anisotropic, curl_incidence, gradient_incidence
    type(csc_t) :: grad_shafranov, toroidal_fourier
    type(csc_t) :: poloidal_bracket
    type(csc_t) :: hcurl_mass, curl_curl, mass, stiffness, weighted_stiffness
    type(csc_t) :: hdiv_mass, div_div
    type(csc_t) :: weak_gradient
    type(csc_t) :: adjoint_gradient
    type(csc_t) :: adjoint_curl, l2_mass, weak_curl
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: coefficients(:), control_points(:, :, :)
    real(dp), allocatable :: edge_values(:), product(:), weights(:, :)
    real(dp), allocatable :: hcurl_coefficients(:)
    real(dp), allocatable :: hdiv_coefficients(:)
    real(dp), allocatable :: l2_coefficients(:)
    real(dp), allocatable :: bracket_field(:), bracket_test(:), bracket_trial(:)
    real(dp), allocatable :: x_points(:), y_points(:)
    real(dp) :: energy, integral
    integer :: ix, iy

    x_points = greville_abscissae(knots_x, 2)
    y_points = greville_abscissae(knots_y, 2)
    allocate(control_points(2, size(x_points), size(y_points)))
    allocate(weights(size(x_points), size(y_points)))
    weights = 1.0_dp
    do iy = 1, size(y_points)
        do ix = 1, size(x_points)
            control_points(:, ix, iy) = [1.0_dp + x_points(ix), y_points(iy)]
        end do
    end do

    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 4, mass, &
        sparse_status, stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    allocate(coefficients(size(x_points)*size(y_points)))
    coefficients = 1.0_dp
    product = csc_matvec(mass, coefficients)
    integral = dot_product(coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(integral - 1.0_dp) < 3.0e-13_dp, &
        "Sparse spline mass integrates a constant exactly")

    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 4, stiffness, &
        sparse_status, stiffness_coefficient=1.0_dp, mass_coefficient=0.0_dp)
    do iy = 1, size(y_points)
        coefficients(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = x_points
    end do
    product = csc_matvec(stiffness, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 1.0_dp) < 4.0e-12_dp, &
        "Sparse isogeometric stiffness gives exact affine-field energy")
    call check_condition(abs(sum(product)) < 3.0e-13_dp, &
        "Spline stiffness annihilates constants in the weak form")

    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 6, &
        weighted_stiffness, sparse_status, stiffness_coefficient=1.0_dp, &
        mass_coefficient=0.0_dp, stiffness_weight_function=one_over_radius)
    do iy = 1, size(y_points)
        coefficients(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = 1.0_dp + x_points
    end do
    product = csc_matvec(weighted_stiffness, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(abs(energy - log(2.0_dp)) < 2.0e-11_dp, &
        "Cylindrical 1/R spline diffusion matches analytical log(2) energy")
    call assemble_bspline_grad_shafranov_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 6, grad_shafranov, &
        sparse_status)
    product = csc_matvec(grad_shafranov, coefficients)
    call check_condition(abs(dot_product(coefficients, product) - &
        log(2.0_dp)) < 2.0e-11_dp, &
        "Named isogeometric Grad-Shafranov operator matches log(2) energy")
    call assemble_bspline_toroidal_fourier_laplacian_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 6, 2, &
        toroidal_fourier, sparse_status)
    coefficients = 1.0_dp
    product = csc_matvec(toroidal_fourier, coefficients)
    call check_condition(abs(dot_product(coefficients, product) - &
        4.0_dp*log(2.0_dp)) < 3.0e-11_dp, &
        "Toroidal Fourier mode energy matches analytical 4 log(2)")
    allocate( &
        bracket_field(size(coefficients)), bracket_test(size(coefficients)), &
        bracket_trial(size(coefficients)))
    do iy = 1, size(y_points)
        bracket_field(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = 1.0_dp + x_points
        bracket_test(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = 1.0_dp + x_points
        bracket_trial(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = y_points(iy)
    end do
    call assemble_bspline_poloidal_bracket_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, bracket_field, 5, &
        poloidal_bracket, sparse_status)
    product = csc_matvec(poloidal_bracket, bracket_trial)
    energy = dot_product(bracket_test, product)
    call check_condition(abs(energy - 0.75_dp) < 3.0e-11_dp, &
        "Poloidal spline bracket reproduces exact affine weak action")
    call check_condition(abs(energy + dot_product(bracket_trial, &
        csc_matvec(poloidal_bracket, bracket_test))) < 3.0e-13_dp, &
        "Poloidal spline bracket is exactly skew in the Galerkin pairing")

    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, anisotropic, &
        sparse_status, stiffness_coefficient=1.0_dp, mass_coefficient=0.0_dp, &
        stiffness_tensor_function=constant_anisotropy)
    do iy = 1, size(y_points)
        coefficients(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = 1.0_dp + x_points + 2.0_dp*y_points(iy)
    end do
    product = csc_matvec(anisotropic, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(abs(energy - 14.0_dp) < 3.0e-11_dp, &
        "Anisotropic spline diffusion reproduces exact tensor energy")

    call build_bspline_feec_2d_operators_csc( &
        knots_x, knots_y, 2, 2, gradient_incidence, curl_incidence, &
        sparse_status)
    coefficients = [(sin(real(3*ix + 1, dp)), ix=1, size(coefficients))]
    edge_values = csc_matvec(gradient_incidence, coefficients)
    product = csc_matvec(curl_incidence, edge_values)
    call check_condition(sparse_status%code == 0 .and. &
        maxval(abs(product)) < 4.0e-14_dp, &
        "FortSparse isogeometric incidence satisfies curl grad equals zero")

    call assemble_bspline_h1_hcurl_gradient_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, weak_gradient, &
        sparse_status)
    do iy = 1, size(y_points)
        coefficients(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = 1.0_dp + x_points
    end do
    edge_values = csc_matvec(gradient_incidence, coefficients)
    product = csc_matvec(weak_gradient, coefficients)
    call check_condition(abs(dot_product(edge_values, product) - 1.0_dp) < &
        3.0e-12_dp, &
        "Mixed spline gradient block reproduces exact affine Hcurl energy")
    call assemble_bspline_hcurl_h1_adjoint_gradient_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, &
        adjoint_gradient, sparse_status)
    call check_condition(abs( &
        dot_product(edge_values, product) - dot_product(coefficients, &
        csc_matvec(adjoint_gradient, edge_values))) < 3.0e-13_dp, &
        "Sparse mixed gradient and adjoint preserve weak duality")

    call assemble_bspline_hcurl_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, hcurl_mass, &
        sparse_status, curl_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    allocate(hcurl_coefficients((size(x_points) - 1)*size(y_points) + &
        size(x_points)*(size(y_points) - 1)))
    hcurl_coefficients = 0.0_dp
    hcurl_coefficients(1:(size(x_points) - 1)*size(y_points)) = 1.0_dp
    product = csc_matvec(hcurl_mass, hcurl_coefficients)
    energy = dot_product(hcurl_coefficients, product)
    call check_condition(abs(energy - 1.0_dp) < 3.0e-12_dp, &
        "Covariant spline Hcurl mass gives exact constant-field energy")

    call assemble_bspline_hcurl_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, curl_curl, &
        sparse_status, curl_coefficient=1.0_dp, mass_coefficient=0.0_dp)
    call set_rotation_coefficients( &
        x_points, y_points, hcurl_coefficients)
    product = csc_matvec(curl_curl, hcurl_coefficients)
    energy = dot_product(hcurl_coefficients, product)
    call check_condition(abs(energy - 4.0_dp) < 3.0e-11_dp, &
        "Spline curl-curl gives exact rigid-rotation curl energy")

    call assemble_bspline_l2_mass_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, l2_mass, &
        sparse_status)
    allocate(l2_coefficients((size(x_points) - 1)*(size(y_points) - 1)))
    l2_coefficients = 1.0_dp
    product = csc_matvec(l2_mass, l2_coefficients)
    call check_condition(abs(dot_product(l2_coefficients, product) - &
        1.0_dp) < 3.0e-12_dp, &
        "Determinant-scaled spline L2 mass integrates unit density exactly")
    call assemble_bspline_hcurl_l2_curl_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, weak_curl, &
        sparse_status)
    edge_values = csc_matvec(curl_incidence, hcurl_coefficients)
    product = csc_matvec(weak_curl, hcurl_coefficients)
    call check_condition(abs(dot_product(edge_values, product) - 4.0_dp) < &
        3.0e-11_dp, "Mixed spline curl block reproduces exact curl energy")
    call assemble_bspline_l2_hcurl_adjoint_curl_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, adjoint_curl, &
        sparse_status)
    call check_condition(abs( &
        dot_product(edge_values, product) - dot_product(hcurl_coefficients, &
        csc_matvec(adjoint_curl, edge_values))) < 3.0e-13_dp, &
        "Sparse mixed curl and adjoint preserve weak duality")

    call assemble_bspline_hdiv_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, hdiv_mass, &
        sparse_status, divergence_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    allocate(hdiv_coefficients(size(x_points)*(size(y_points) - 1) + &
        (size(x_points) - 1)*size(y_points)))
    hdiv_coefficients = 0.0_dp
    hdiv_coefficients(1:size(x_points)*(size(y_points) - 1)) = 1.0_dp
    product = csc_matvec(hdiv_mass, hdiv_coefficients)
    energy = dot_product(hdiv_coefficients, product)
    call check_condition(abs(energy - 1.0_dp) < 3.0e-12_dp, &
        "Contravariant spline Hdiv mass gives exact constant-flux energy")

    call assemble_bspline_hdiv_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 5, div_div, &
        sparse_status, divergence_coefficient=1.0_dp, mass_coefficient=0.0_dp)
    call set_affine_flux_coefficients(x_points, y_points, hdiv_coefficients)
    product = csc_matvec(div_div, hdiv_coefficients)
    energy = dot_product(hdiv_coefficients, product)
    call check_condition(abs(energy - 4.0_dp) < 3.0e-11_dp, &
        "Spline div-div gives exact affine-flux divergence energy")

    call check_summary("Sparse isogeometric H1 assembly")

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

    pure subroutine one_over_radius(point, value)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: value

        value = 1.0_dp/point(1)
    end subroutine one_over_radius

    pure subroutine constant_anisotropy(point, value)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: value(2, 2)

        value = reshape([2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2])
        associate (unused => point)
        end associate
    end subroutine constant_anisotropy

    pure subroutine set_rotation_coefficients(x, y, coefficients_local)
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: coefficients_local(:)
        integer :: basis_x, basis_y, offset

        offset = (size(x) - 1)*size(y)
        do basis_y = 1, size(y)
            do basis_x = 1, size(x) - 1
                coefficients_local(basis_x + &
                    (basis_y - 1)*(size(x) - 1)) = -y(basis_y)
            end do
        end do
        do basis_y = 1, size(y) - 1
            do basis_x = 1, size(x)
                coefficients_local(offset + basis_x + &
                    (basis_y - 1)*size(x)) = 1.0_dp + x(basis_x)
            end do
        end do
    end subroutine set_rotation_coefficients

    pure subroutine set_affine_flux_coefficients(x, y, coefficients_local)
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: coefficients_local(:)
        integer :: basis_x, basis_y, offset

        offset = size(x)*(size(y) - 1)
        do basis_y = 1, size(y) - 1
            do basis_x = 1, size(x)
                coefficients_local(basis_x + (basis_y - 1)*size(x)) = &
                    1.0_dp + x(basis_x)
            end do
        end do
        do basis_y = 1, size(y)
            do basis_x = 1, size(x) - 1
                coefficients_local(offset + basis_x + &
                    (basis_y - 1)*(size(x) - 1)) = y(basis_y)
            end do
        end do
    end subroutine set_affine_flux_coefficients

end program test_bspline_h1_sparse_assembly
