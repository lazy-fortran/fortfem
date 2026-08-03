program test_api07_tensor_consumer
    !! Downstream-style smoke for tensor-valued pressure and constitutive APIs.
    !!
    !! The client imports only the canonical ``fortfem_core`` and
    !! ``fortfem_feec`` facades.  The value, directional derivative, and
    !! transpose checks below are written from the defining contractions, so
    !! this test does not merely compare two implementation paths.
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, unit_square_mesh
    use fortfem_feec, only: &
        evaluate_cgl_pressure_tensor, evaluate_cgl_pressure_tensor_jvp, &
        evaluate_cgl_pressure_traction, evaluate_cgl_pressure_work, &
        evaluate_cgl_pressure_work_jvp, evaluate_cgl_pressure_work_vjp, &
        evaluate_field_aligned_constitutive_tensor, &
        evaluate_field_aligned_constitutive_tensor_jvp, &
        evaluate_field_aligned_constitutive_tensor_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    call check_core_facade()
    call check_pressure_stress()
    call check_field_aligned_anisotropy()
    call check_summary("API-07 tensor canonical consumer")

contains

    subroutine check_core_facade()
        type(mesh_t) :: mesh

        mesh = unit_square_mesh(2)
        call check_condition(mesh%data%n_vertices == 4, &
            "tensor consumer can construct a canonical core mesh")
    end subroutine check_core_facade

    subroutine check_pressure_stress()
        real(dp), parameter :: p_parallel = 3.0_dp
        real(dp), parameter :: p_perpendicular = 1.2_dp
        real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
        real(dp), parameter :: direction_dot(3) = [0.2_dp, -0.15_dp, 0.0_dp]
        real(dp), parameter :: p_parallel_dot = 0.25_dp
        real(dp), parameter :: p_perpendicular_dot = -0.1_dp
        real(dp), parameter :: normal(3) = [0.5_dp, -0.2_dp, 0.7_dp]
        real(dp), parameter :: gradient(3, 3) = reshape([ &
            1.0_dp, -0.5_dp, 0.7_dp, 0.2_dp, 0.4_dp, -0.3_dp, &
            -0.6_dp, 0.8_dp, 0.9_dp], [3, 3])
        real(dp), parameter :: gradient_dot(3, 3) = reshape([ &
            -0.2_dp, 0.1_dp, 0.4_dp, -0.3_dp, 0.6_dp, 0.2_dp, &
            0.5_dp, -0.4_dp, 0.3_dp], [3, 3])
        real(dp), parameter :: work_bar = 0.75_dp
        real(dp) :: pressure(3, 3), pressure_dot(3, 3), expected(3, 3)
        real(dp) :: expected_dot(3, 3), traction(3), expected_traction(3)
        real(dp) :: work, work_dot, expected_work, expected_work_dot
        real(dp) :: p_parallel_bar, p_perpendicular_bar, direction_bar(3)
        real(dp) :: gradient_bar(3, 3), lhs, rhs
        integer :: row, column
        type(fortsparse_status_t) :: status

        call pressure_oracle(p_parallel, p_perpendicular, direction, expected)
        call pressure_oracle_jvp( &
            p_parallel, p_perpendicular, direction, p_parallel_dot, &
            p_perpendicular_dot, direction_dot, expected_dot)

        call evaluate_cgl_pressure_tensor( &
            p_parallel, p_perpendicular, direction, pressure, status)
        call check_condition(status%code == 0 .and. &
            maxval(abs(pressure - expected)) < 2.0e-14_dp, &
            "tensor consumer matches the independent gyrotropic pressure oracle")

        call evaluate_cgl_pressure_tensor_jvp( &
            p_parallel, p_perpendicular, direction, p_parallel_dot, &
            p_perpendicular_dot, direction_dot, pressure_dot, status)
        call check_condition(status%code == 0 .and. &
            maxval(abs(pressure_dot - expected_dot)) < 2.0e-13_dp, &
            "tensor consumer exposes the independent pressure JVP")

        call evaluate_cgl_pressure_traction( &
            p_parallel, p_perpendicular, direction, normal, traction, status)
        expected_traction = matmul(expected, normal)
        call check_condition(status%code == 0 .and. &
            maxval(abs(traction - expected_traction)) < 2.0e-14_dp, &
            "stress traction is the pressure tensor contraction")

        call evaluate_cgl_pressure_work( &
            p_parallel, p_perpendicular, direction, gradient, work, status)
        expected_work = sum(expected*gradient)
        call check_condition(status%code == 0 .and. &
            abs(work - expected_work) < 2.0e-14_dp, &
            "tensor-valued stress work matches the independent contraction")

        call evaluate_cgl_pressure_work_jvp( &
            p_parallel, p_perpendicular, direction, gradient, p_parallel_dot, &
            p_perpendicular_dot, direction_dot, gradient_dot, work_dot, status)
        expected_work_dot = sum(expected_dot*gradient) + sum(expected*gradient_dot)
        call check_condition(status%code == 0 .and. &
            abs(work_dot - expected_work_dot) < 2.0e-13_dp, &
            "tensor-valued stress work exposes the product-rule JVP")

        call evaluate_cgl_pressure_work_vjp( &
            p_parallel, p_perpendicular, direction, gradient, work_bar, &
            p_parallel_bar, p_perpendicular_bar, direction_bar, gradient_bar, &
            status)
        lhs = work_bar*work_dot
        rhs = p_parallel_bar*p_parallel_dot + &
            p_perpendicular_bar*p_perpendicular_dot + &
            dot_product(direction_bar, direction_dot) + sum(gradient_bar*gradient_dot)
        call check_condition(status%code == 0 .and. abs(lhs - rhs) < 2.0e-12_dp, &
            "tensor-valued stress work satisfies the independent VJP ledger")

        ! Keep the oracle visibly component-wise, including its symmetric part.
        do row = 1, 3
            do column = 1, 3
                if (row == column) then
                    call check_condition(abs(expected(row, column) - &
                        (p_perpendicular + (p_parallel-p_perpendicular)* &
                        direction(row)*direction(column))) < 2.0e-14_dp, &
                        "pressure oracle retains the isotropic diagonal term")
                end if
            end do
        end do
    end subroutine check_pressure_stress

    subroutine check_field_aligned_anisotropy()
        real(dp), parameter :: parallel = 3.0_dp
        real(dp), parameter :: perpendicular = 1.0_dp
        real(dp), parameter :: hall = 0.5_dp
        real(dp), parameter :: parallel_dot = 0.2_dp
        real(dp), parameter :: perpendicular_dot = -0.1_dp
        real(dp), parameter :: hall_dot = -0.3_dp
        real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
        real(dp), parameter :: direction_dot(3) = [0.2_dp, -0.15_dp, 0.0_dp]
        real(dp), parameter :: tensor_bar(3, 3) = reshape([ &
            0.3_dp, -0.1_dp, 0.4_dp, 0.2_dp, 0.5_dp, -0.6_dp, &
            -0.7_dp, 0.8_dp, 0.7_dp], [3, 3])
        real(dp) :: tensor(3, 3), tensor_dot(3, 3), expected(3, 3)
        real(dp) :: expected_dot(3, 3), parallel_bar, perpendicular_bar
        real(dp) :: hall_bar, direction_bar(3), lhs, rhs
        real(dp) :: pressure(3, 3), pressure_dot(3, 3), skew(3, 3)
        type(fortsparse_status_t) :: status

        call pressure_oracle(parallel, perpendicular, direction, pressure)
        call pressure_oracle_jvp( &
            parallel, perpendicular, direction, parallel_dot, &
            perpendicular_dot, direction_dot, pressure_dot)
        call skew_oracle(hall, direction, skew)
        expected = pressure + skew
        call skew_oracle_jvp(hall, direction, hall_dot, direction_dot, skew)
        expected_dot = pressure_dot + skew

        call evaluate_field_aligned_constitutive_tensor( &
            parallel, perpendicular, direction, tensor, status, hall)
        call check_condition(status%code == 0 .and. &
            maxval(abs(tensor - expected)) < 2.0e-14_dp, &
            "field-aligned tensor consumer matches the anisotropic oracle")

        call evaluate_field_aligned_constitutive_tensor_jvp( &
            parallel, perpendicular, direction, parallel_dot, &
            perpendicular_dot, direction_dot, tensor_dot, status, hall, hall_dot)
        call check_condition(status%code == 0 .and. &
            maxval(abs(tensor_dot - expected_dot)) < 2.0e-13_dp, &
            "field-aligned tensor consumer exposes the Hall-aware JVP")

        call evaluate_field_aligned_constitutive_tensor_vjp( &
            parallel, perpendicular, direction, tensor_bar, parallel_bar, &
            perpendicular_bar, direction_bar, status, hall, hall_bar)
        lhs = sum(tensor_bar*tensor_dot)
        rhs = parallel_bar*parallel_dot + perpendicular_bar*perpendicular_dot + &
            hall_bar*hall_dot + dot_product(direction_bar, direction_dot)
        call check_condition(status%code == 0 .and. abs(lhs - rhs) < 2.0e-12_dp, &
            "field-aligned tensor consumer satisfies the independent VJP ledger")
    end subroutine check_field_aligned_anisotropy

    pure subroutine pressure_oracle( &
            p_parallel, p_perpendicular, direction, tensor)
        real(dp), intent(in) :: p_parallel, p_perpendicular, direction(3)
        real(dp), intent(out) :: tensor(3, 3)
        integer :: row, column

        do column = 1, 3
            do row = 1, 3
                tensor(row, column) = (p_parallel-p_perpendicular)* &
                    direction(row)*direction(column)
                if (row == column) tensor(row, column) = &
                    tensor(row, column) + p_perpendicular
            end do
        end do
    end subroutine pressure_oracle

    pure subroutine pressure_oracle_jvp( &
            p_parallel, p_perpendicular, direction, p_parallel_dot, &
            p_perpendicular_dot, direction_dot, tensor_dot)
        real(dp), intent(in) :: p_parallel, p_perpendicular, direction(3)
        real(dp), intent(in) :: p_parallel_dot, p_perpendicular_dot
        real(dp), intent(in) :: direction_dot(3)
        real(dp), intent(out) :: tensor_dot(3, 3)
        integer :: row, column

        do column = 1, 3
            do row = 1, 3
                tensor_dot(row, column) = (p_parallel_dot- &
                    p_perpendicular_dot)*direction(row)*direction(column) + &
                    (p_parallel-p_perpendicular)*(direction_dot(row)* &
                    direction(column) + direction(row)*direction_dot(column))
                if (row == column) tensor_dot(row, column) = &
                    tensor_dot(row, column) + p_perpendicular_dot
            end do
        end do
    end subroutine pressure_oracle_jvp

    pure subroutine skew_oracle(hall, direction, skew)
        real(dp), intent(in) :: hall, direction(3)
        real(dp), intent(out) :: skew(3, 3)

        skew = 0.0_dp
        skew(1, 2) = -hall*direction(3)
        skew(1, 3) = hall*direction(2)
        skew(2, 1) = hall*direction(3)
        skew(2, 3) = -hall*direction(1)
        skew(3, 1) = -hall*direction(2)
        skew(3, 2) = hall*direction(1)
    end subroutine skew_oracle

    pure subroutine skew_oracle_jvp( &
            hall, direction, hall_dot, direction_dot, skew_dot)
        real(dp), intent(in) :: hall, direction(3), hall_dot, direction_dot(3)
        real(dp), intent(out) :: skew_dot(3, 3)

        call skew_oracle(hall_dot, direction, skew_dot)
        skew_dot = skew_dot + skew_from_direction(hall, direction_dot)
    end subroutine skew_oracle_jvp

    pure function skew_from_direction(hall, direction) result(skew)
        real(dp), intent(in) :: hall, direction(3)
        real(dp) :: skew(3, 3)

        call skew_oracle(hall, direction, skew)
    end function skew_from_direction

end program test_api07_tensor_consumer
