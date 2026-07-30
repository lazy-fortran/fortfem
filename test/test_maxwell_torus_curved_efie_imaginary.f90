program test_maxwell_torus_curved_efie_imaginary
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_efie_imaginary_rwg_3d, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: decay_rate = 0.6_dp, impedance = 1.7_dp
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    complex(dp), allocatable :: matrix(:, :), scaled(:, :)
    complex(dp), allocatable :: vector(:)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: scaled_vertices(:, :)
    complex(dp) :: energy
    real(dp) :: scaling_error, symmetry_error
    integer :: basis, status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 4, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        decay_rate, impedance, 3, 1.0e-5_dp, 1, matrix, status)
    allocate(vector(size(matrix, 1)))
    do basis = 1, size(vector)
        vector(basis) = cmplx(cos(real(2*basis, dp)), 0.0_dp, dp)
    end do
    energy = sum(vector*matmul(matrix, vector))
    symmetry_error = maxval(abs(matrix - transpose(matrix)))
    call record_condition(status == 0 .and. symmetry_error < 3.0e-13_dp .and. &
        maxval(abs(aimag(matrix))) < 3.0e-14_dp, &
        "toroidal imaginary-wavenumber EFIE is real and reciprocal")
    call record_condition(real(energy, dp) < 0.0_dp .and. &
        abs(aimag(energy)) < 3.0e-13_dp*abs(real(energy, dp)), &
        "toroidal Yukawa EFIE has coercive signed energy")

    scaled_vertices = 2.0_dp*vertices
    call assemble_maxwell_torus_curved_efie_imaginary_rwg_3d( &
        scaled_vertices, triangles, parameters, 2.0_dp*major_radius, &
        2.0_dp*minor_radius, decay_rate/2.0_dp, impedance, 3, 1.0e-5_dp, 1, &
        scaled, status)
    scaling_error = maxval(abs(scaled - 4.0_dp*matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. scaling_error < 2.0e-13_dp, &
        "toroidal Yukawa EFIE obeys analytical area scaling")

    call check_summary("Exact-curved torus imaginary-wavenumber EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_efie_imaginary
