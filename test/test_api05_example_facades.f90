program test_api05_example_facades
    !! Behavioral oracle for the four migrated gallery capability families.
    !!
    !! The expected values below are independent analytical identities: a
    !! constant Nedelec field, the exact logarithmic BEM self integral, a
    !! circular DtN eigenmode, and a reversible harmonic midpoint step.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        apply_circular_helmholtz_dtn, &
        assemble_laplace_single_layer_constant, &
        circular_helmholtz_dtn_eigenvalue
    use fortfem_feec, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, &
        interpolate_reference_tetra_nedelec, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    use fortfem_time, only: advance_mixed_wave_midpoint
    use fortsparse, only: fortsparse_status_t
    implicit none

    call check_nedelec_solution()
    call check_bem_solution()
    call check_dtn_solution()
    call check_wave_solution()
    call check_summary("API-05 representative example facades")

contains

    subroutine check_nedelec_solution()
        type(tetra_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: dofs(:)
        real(dp), allocatable :: values(:, :), curls(:, :)
        real(dp) :: point(3), field_value(3)
        integer :: status, basis_status

        call initialize_tetra_nedelec_first_kind(1, basis, basis_status)
        call check_condition(basis_status == 0 .and. &
            tetra_nedelec_dof_count(basis) == 6, &
            "FEEC facade exposes the order-one Nedelec space")
        call interpolate_reference_tetra_nedelec( &
            basis, constant_vector, dofs, status)
        allocate(values(3, size(dofs)), curls(3, size(dofs)))
        point = [0.2_dp, 0.3_dp, 0.1_dp]
        call evaluate_tetra_nedelec_first_kind( &
            basis, point, values, curls, status)
        field_value = matmul(values, dofs)
        call check_condition(status == 0 .and. &
            maxval(abs(field_value - [1.0_dp, 0.0_dp, 0.0_dp])) < 2.0e-12_dp, &
            "Nedelec facade reproduces an independent constant vector field")
    end subroutine check_nedelec_solution

    subroutine check_bem_solution()
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: starts(2, 1), ends(2, 1), matrix(1, 1), expected
        integer :: status

        starts(:, 1) = [0.0_dp, 0.0_dp]
        ends(:, 1) = [2.0_dp, 0.0_dp]
        call assemble_laplace_single_layer_constant( &
            starts, ends, 24, matrix, status)
        expected = 4.0_dp*(1.5_dp - log(2.0_dp))/(2.0_dp*pi)
        call check_condition(status == 0 .and. &
            abs(matrix(1, 1) - expected) < 2.0e-14_dp, &
            "boundary facade BEM agrees with the logarithmic self integral")
    end subroutine check_bem_solution

    subroutine check_dtn_solution()
        integer, parameter :: point_count = 8
        real(dp), parameter :: wavenumber = 1.7_dp, radius = 1.25_dp
        complex(dp) :: trace(point_count), derivative(point_count)
        complex(dp) :: eigenvalue
        real(dp) :: discarded
        integer :: status, eigenvalue_status

        trace = cmplx(1.0_dp, 0.0_dp, dp)
        call circular_helmholtz_dtn_eigenvalue( &
            0, wavenumber, radius, eigenvalue, eigenvalue_status)
        call apply_circular_helmholtz_dtn( &
            trace, wavenumber, radius, 0, derivative, discarded, status)
        call check_condition(status == 0 .and. eigenvalue_status == 0 .and. &
            maxval(abs(derivative - eigenvalue*trace)) < 2.0e-12_dp .and. &
            discarded < 2.0e-14_dp, &
            "boundary facade DtN maps the constant circular mode exactly")
    end subroutine check_dtn_solution

    subroutine check_wave_solution()
        real(dp), parameter :: time_step = 0.1_dp
        real(dp) :: mass(1, 1), coupling(1, 1), q(1), v(1)
        real(dp) :: q_initial(1), v_initial(1), energy_initial, energy
        type(fortsparse_status_t) :: status

        mass = 1.0_dp
        coupling = 1.4_dp
        q = [1.0_dp]
        v = [-0.3_dp]
        q_initial = q
        v_initial = v
        energy_initial = 0.5_dp*(q(1)**2 + v(1)**2)
        call advance_mixed_wave_midpoint( &
            mass, mass, coupling, time_step, q, v, status)
        energy = 0.5_dp*(q(1)**2 + v(1)**2)
        call check_condition(status%code == 0 .and. &
            abs(energy - energy_initial) < 2.0e-13_dp, &
            "time facade midpoint preserves the wave energy")
        call advance_mixed_wave_midpoint( &
            mass, mass, coupling, -time_step, q, v, status)
        call check_condition(status%code == 0 .and. &
            maxval(abs(q - q_initial)) < 2.0e-13_dp .and. &
            maxval(abs(v - v_initial)) < 2.0e-13_dp, &
            "time facade midpoint is reversible for the wave mode")
    end subroutine check_wave_solution

    pure subroutine constant_vector(point, value)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: value(3)

        value = [1.0_dp, 0.0_dp, 0.0_dp]
    end subroutine constant_vector

end program test_api05_example_facades
