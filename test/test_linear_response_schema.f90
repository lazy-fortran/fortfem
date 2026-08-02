program test_linear_response_schema
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        initialize_linear_response_interchange, &
        linear_response_schema_magic, &
        linear_response_interchange_t, &
        read_linear_response_interchange, &
        validate_linear_response_interchange, &
        write_linear_response_interchange
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: n = 2, mode_count = 2, response_count = 1
    character(len=*), parameter :: filename = "linear_response_schema_test.dat"
    character(len=*), parameter :: bad_filename = "linear_response_schema_bad.dat"
    type(linear_response_interchange_t) :: original, restored
    complex(dp) :: equilibrium(n, n), inertia(n, n), resistive(n, n)
    complex(dp) :: vacuum(n, n), wall(n, n), response(response_count, response_count)
    complex(dp) :: frequency
    integer :: mode_numbers(2, mode_count), status, unit
    character(len=32) :: mode_labels(mode_count), response_labels(response_count)
    logical :: valid

    equilibrium = reshape([ &
        cmplx(2.0_dp, 0.1_dp, dp), cmplx(0.3_dp, -0.2_dp, dp), &
        cmplx(-0.4_dp, 0.5_dp, dp), cmplx(1.2_dp, 0.0_dp, dp)], shape(equilibrium))
    inertia = reshape([ &
        cmplx(0.7_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.1_dp, dp), &
        cmplx(0.0_dp, -0.1_dp, dp), cmplx(0.9_dp, 0.0_dp, dp)], shape(inertia))
    resistive = 0.2_dp*inertia
    vacuum = 0.1_dp*equilibrium
    wall = reshape([ &
        cmplx(0.4_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp), cmplx(0.3_dp, 0.0_dp, dp)], shape(wall))
    response(1, 1) = cmplx(0.8_dp, -0.15_dp, dp)
    frequency = cmplx(1.2_dp, -0.08_dp, dp)
    mode_numbers = reshape([0, 1, 1, -1], shape(mode_numbers))
    mode_labels = [character(len=32) :: "m0_n1", "m1_nminus1"]
    response_labels = [character(len=32) :: "wall_trace"]

    call initialize_linear_response_interchange( &
        original, frequency, mode_numbers, mode_labels, equilibrium, inertia, &
        resistive, vacuum, wall, response, response_labels, 3.5_dp, status)
    original%producer = "neutral-wall-fixture"
    original%provenance = "independent-schema-oracle"
    call check_condition(status == 0, "schema fixture initializes")
    call write_linear_response_interchange(filename, original, status)
    call check_condition(status == 0, "versioned response schema writes")
    call read_linear_response_interchange(filename, restored, status)
    call check_condition(status == 0, "versioned response schema reads")
    valid = validate_linear_response_interchange(restored, status)
    call check_condition(valid .and. status == 0, "round-tripped response validates")
    call check_condition(restored%schema_version == original%schema_version, &
        "schema version is retained")
    call check_condition(restored%producer == original%producer .and. &
        restored%provenance == original%provenance, &
        "provenance metadata is retained")
    call check_condition(restored%mode_count == original%mode_count .and. &
        restored%response_count == original%response_count .and. &
        maxval(abs(restored%response_matrix - original%response_matrix)) < 1.0e-14_dp, &
        "response dimensions and matrix are retained")
    call check_condition(all(restored%mode_numbers == original%mode_numbers) .and. &
        all(restored%mode_labels == original%mode_labels) .and. &
        all(restored%response_labels == original%response_labels), &
        "labels and Fourier mode numbers are retained")
    call check_condition(maxval(abs(restored%equilibrium_block - &
        original%equilibrium_block)) < 1.0e-14_dp .and. &
        maxval(abs(restored%wall_block - original%wall_block)) < 1.0e-14_dp .and. &
        abs(restored%frequency - original%frequency) < 1.0e-14_dp .and. &
        abs(restored%normalization_scale - original%normalization_scale) < 1.0e-14_dp, &
        "complex blocks and normalization are retained")

    open(newunit=unit, file=filename, status="old", action="read")
    close(unit, status="delete")
    open(newunit=unit, file=bad_filename, status="replace", action="write")
    write(unit, '(A)') linear_response_schema_magic
    write(unit, '(A)') "fortfem-linear-response-1"
    write(unit, '(A)') "producer"
    write(unit, '(A)') "provenance"
    write(unit, *) 0.0_dp, 0.0_dp, 1.0_dp
    write(unit, *) 1, 1, 100000000
    close(unit)
    call read_linear_response_interchange(bad_filename, restored, status)
    call check_condition(status /= 0, "oversized schema is rejected before allocation")
    open(newunit=unit, file=bad_filename, status="old", action="read")
    close(unit, status="delete")
    call check_summary("linear response schema")
end program test_linear_response_schema
