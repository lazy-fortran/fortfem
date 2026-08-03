program test_curvilinear_pml_geometry_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        build_curvilinear_normal_pml_element_stretch, &
        build_curvilinear_normal_pml_element_stretch_jvp, &
        build_curvilinear_normal_pml_element_stretch_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp), parameter :: vertices(3, 5) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    integer, parameter :: cells(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: layer_origins(3, 2) = reshape([ &
        0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp], [3, 2])
    real(dp), parameter :: layer_normals(3, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.6_dp, 0.8_dp, 0.0_dp], [3, 2])
    real(dp), parameter :: layer_width(2) = [0.6_dp, 0.8_dp]
    real(dp), parameter :: wave_number = 2.0_dp
    real(dp), parameter :: sigma_max = 5.0_dp
    integer, parameter :: polynomial_degree = 2
    real(dp) :: vertices_dot(3, 5), origins_dot(3, 2)
    real(dp) :: normals_dot(3, 2), width_dot(2)
    real(dp) :: wave_number_dot, sigma_max_dot
    real(dp) :: vertices_bar(3, 5), origins_bar(3, 2)
    real(dp) :: normals_bar(3, 2), width_bar(2)
    real(dp) :: wave_number_bar, sigma_max_bar
    complex(dp), allocatable :: stretch(:, :, :), stretch_dot(:, :, :)
    complex(dp), allocatable :: stretch_minus(:, :, :), stretch_plus(:, :, :)
    complex(dp), allocatable :: stretch_bar(:, :, :)
    real(dp) :: forward_pairing, reverse_pairing
    integer :: status, status_minus, status_plus, i
    logical :: all_passed

    all_passed = .true.
    vertices_dot = reshape([ &
        (0.01_dp*cos(0.2_dp*real(i, dp)), i=1, 15)], [3, 5])
    origins_dot = reshape([ &
        (0.01_dp*sin(0.3_dp*real(i, dp)), i=1, 6)], [3, 2])
    normals_dot = reshape([ &
        0.0_dp, 0.03_dp, 0.0_dp, 0.08_dp, -0.06_dp, 0.0_dp], [3, 2])
    width_dot = [0.04_dp, -0.03_dp]
    wave_number_dot = 0.07_dp
    sigma_max_dot = -0.11_dp

    call build_curvilinear_normal_pml_element_stretch( &
        vertices, cells, layer_origins, layer_normals, layer_width, &
        wave_number, sigma_max, polynomial_degree, stretch, status)
    call record_condition(status == 0, &
        "curved normal-frame PML geometry is accepted")
    call record_condition(abs(stretch(2, 2, 1) - &
        cmplx(1.0_dp, 0.0_dp, dp)) < 1.0e-14_dp, &
        "normal-frame PML preserves the tangential diagonal")
    call record_condition(abs(stretch(1, 2, 2) - &
        cmplx(0.0_dp, 0.0_dp, dp)) > 1.0e-8_dp .and. &
        abs(stretch(1, 2, 2) - stretch(2, 1, 2)) < 1.0e-14_dp, &
        "curved normal-frame PML produces symmetric off-diagonal stretch")

    call build_curvilinear_normal_pml_element_stretch_jvp( &
        vertices, cells, layer_origins, layer_normals, layer_width, &
        wave_number, sigma_max, polynomial_degree, vertices_dot, origins_dot, &
        normals_dot, width_dot, wave_number_dot, sigma_max_dot, stretch_dot, &
        status)
    call build_curvilinear_normal_pml_element_stretch( &
        vertices - difference_step*vertices_dot, cells, &
        layer_origins - difference_step*origins_dot, &
        layer_normals - difference_step*normals_dot, &
        layer_width - difference_step*width_dot, &
        wave_number - difference_step*wave_number_dot, &
        sigma_max - difference_step*sigma_max_dot, polynomial_degree, &
        stretch_minus, status_minus)
    call build_curvilinear_normal_pml_element_stretch( &
        vertices + difference_step*vertices_dot, cells, &
        layer_origins + difference_step*origins_dot, &
        layer_normals + difference_step*normals_dot, &
        layer_width + difference_step*width_dot, &
        wave_number + difference_step*wave_number_dot, &
        sigma_max + difference_step*sigma_max_dot, polynomial_degree, &
        stretch_plus, status_plus)
    call record_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "curved normal-frame PML JVP inputs are accepted")
    call record_condition(maxval(abs(stretch_dot - &
        (stretch_plus - stretch_minus)/(2.0_dp*difference_step))) < 3.0e-9_dp, &
        "curved normal-frame PML JVP matches central difference")

    allocate(stretch_bar(3, 3, 2))
    stretch_bar = reshape([ &
        (cmplx(0.01_dp*real(i, dp), -0.02_dp*real(i + 1, dp), dp), &
        i=1, 18)], [3, 3, 2])
    call build_curvilinear_normal_pml_element_stretch_vjp( &
        vertices, cells, layer_origins, layer_normals, layer_width, &
        wave_number, sigma_max, polynomial_degree, stretch_bar, vertices_bar, &
        origins_bar, normals_bar, width_bar, wave_number_bar, sigma_max_bar, &
        status)
    forward_pairing = real(sum(conjg(stretch_bar)*stretch_dot), dp)
    reverse_pairing = sum(vertices_bar*vertices_dot) + &
        sum(origins_bar*origins_dot) + sum(normals_bar*normals_dot) + &
        dot_product(width_bar, width_dot) + wave_number_bar*wave_number_dot + &
        sigma_max_bar*sigma_max_dot
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-11_dp, &
        "curved normal-frame PML JVP and VJP satisfy the dot identity")

    call build_curvilinear_normal_pml_element_stretch( &
        vertices, cells, layer_origins, reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.6_dp, 0.7_dp, 0.0_dp], [3, 2]), &
        layer_width, wave_number, sigma_max, polynomial_degree, &
        stretch_minus, status)
    call record_condition(status /= 0, "non-unit curved layer normals are rejected")
    call check_summary("Curvilinear PML geometry derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_curvilinear_pml_geometry_ad
