program test_maxwell_surface_rt_basis
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_maxwell_surface_rt_basis, &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_rt_basis_t) :: basis
    real(dp) :: divergences(3), flux_moments(3, 3), values(3, 3)
    real(dp), parameter :: edge_midpoints(2, 3) = reshape([ &
        0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp], [2, 3])
    real(dp), parameter :: outward_conormals(3, 3) = reshape([ &
        0.0_dp, -1.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, &
        -1.0_dp, 0.0_dp, 0.0_dp], [3, 3])
    integer :: edge, status

    call initialize_triangle_raviart_thomas(0, basis, status)
    if (status /= 0) error stop "surface RT basis initialization failed"
    do edge = 1, 3
        call evaluate_maxwell_surface_rt_basis( &
            basis, edge_midpoints(1, edge), edge_midpoints(2, edge), &
            [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], &
            1.0_dp, values, divergences, status)
        if (status /= 0) error stop "surface RT basis evaluation failed"
        flux_moments(edge, :) = matmul( &
            transpose(values), outward_conormals(:, edge))
    end do
    write (*, "(a,9(1x,es10.2))") "surface RT0 edge moments:", flux_moments
    call check_condition( &
        maxval(abs(flux_moments - identity3())) < 2.0e-12_dp, &
        "Surface RT0 basis has Kronecker normal-edge fluxes")
    call check_summary("Higher-order Maxwell surface RT basis")

contains

    pure function identity3() result(matrix)
        real(dp) :: matrix(3, 3)

        matrix = 0.0_dp
        matrix(1, 1) = 1.0_dp
        matrix(2, 2) = 1.0_dp
        matrix(3, 3) = 1.0_dp
    end function identity3

end program test_maxwell_surface_rt_basis
