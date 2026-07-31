program test_maxwell_cfie_resonance_3d_slow
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_efie_rwg_3d, &
        assemble_maxwell_regularized_cfie_rwg_3d, &
        assemble_maxwell_sphere_curved_efie_rwg_3d, &
        assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: efie(:, :), matrix(:, :), mfie(:, :)
    complex(dp), allocatable :: product(:, :), regularizer(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: cfie_condition, condition, curved_cfie_condition
    real(dp) :: curved_efie_condition, curved_peak_wave_number
    real(dp) :: efie_condition, peak_wave_number
    real(dp) :: wave_number
    integer :: sample, status
    logical :: all_passed

    interface
        subroutine zgesvd( &
                jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
                lwork, rwork, info)
            import :: dp
            character, intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            complex(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: s(*), rwork(*)
            complex(dp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
            integer, intent(out) :: info
        end subroutine zgesvd
    end interface

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 1, vertices, triangles)
    efie_condition = 0.0_dp
    peak_wave_number = 0.0_dp
    do sample = -4, 4
        wave_number = 2.743707269992269_dp + 0.1_dp*real(sample, dp)
        call assemble_maxwell_efie_rwg_3d( &
            vertices, triangles, wave_number, 1.0_dp, 2, 1.0e-5_dp, 1, &
            efie, status)
        condition = spectral_condition_number(efie, status)
        if (condition > efie_condition) then
            efie_condition = condition
            peak_wave_number = wave_number
        end if
    end do
    call assemble_maxwell_regularized_cfie_rwg_3d( &
        vertices, triangles, peak_wave_number, 1.0_dp, 2, 1.0e-5_dp, 1, &
        matrix, efie, mfie, regularizer, product, status)
    cfie_condition = spectral_condition_number(matrix, status)
    call record_condition(status == 0 .and. &
        abs(peak_wave_number - 2.743707269992269_dp) <= 0.4_dp .and. &
        cfie_condition < 0.2_dp*efie_condition .and. cfie_condition < 30.0_dp, &
        "regularized CFIE suppresses the first sphere-resonance condition peak")

    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    curved_efie_condition = 0.0_dp
    curved_peak_wave_number = 0.0_dp
    do sample = -4, 4
        wave_number = 2.743707269992269_dp + 0.1_dp*real(sample, dp)
        call assemble_maxwell_sphere_curved_efie_rwg_3d( &
            vertices, triangles, 1.0_dp, wave_number, 1.0_dp, 2, 2.0e-4_dp, 1, &
            efie, status)
        condition = spectral_condition_number(efie, status)
        if (condition > curved_efie_condition) then
            curved_efie_condition = condition
            curved_peak_wave_number = wave_number
        end if
    end do
    call assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d( &
        vertices, triangles, 1.0_dp, curved_peak_wave_number, 1.0_dp, 2, &
        2.0e-4_dp, 1, 0.18_dp, matrix, efie, mfie, regularizer, product, status)
    curved_cfie_condition = spectral_condition_number(matrix, status)
    if (status /= 0 .or. curved_cfie_condition >= &
        0.8_dp*curved_efie_condition) write (*, *) &
        "curved EFIE/CFIE conditions", curved_peak_wave_number, &
        curved_efie_condition, curved_cfie_condition
    call record_condition(status == 0 .and. &
        abs(curved_peak_wave_number - 2.743707269992269_dp) <= 0.4_dp .and. &
        curved_cfie_condition < 0.8_dp*curved_efie_condition, &
        "curved regularized CFIE suppresses the sphere-resonance condition peak")

    call check_summary("Maxwell CFIE sphere resonance")
    if (.not. all_passed) error stop 1

contains

    function spectral_condition_number(input, local_status) result(condition)
        complex(dp), intent(in) :: input(:, :)
        integer, intent(out) :: local_status
        real(dp) :: condition

        complex(dp), allocatable :: matrix_copy(:, :), work(:)
        complex(dp) :: unused_u(1, 1), unused_vt(1, 1)
        real(dp), allocatable :: singular_values(:), real_work(:)
        integer :: info, n

        local_status = 1
        condition = huge(1.0_dp)
        if (size(input, 1) /= size(input, 2)) return
        n = size(input, 1)
        allocate( &
            matrix_copy(n, n), singular_values(n), real_work(5*n), work(4*n))
        matrix_copy = input
        call zgesvd( &
            "N", "N", n, n, matrix_copy, n, singular_values, unused_u, 1, &
            unused_vt, 1, work, size(work), real_work, info)
        if (info /= 0) return
        if (singular_values(n) <= tiny(1.0_dp)) return
        condition = singular_values(1)/singular_values(n)
        local_status = 0
    end function spectral_condition_number

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_cfie_resonance_3d_slow
