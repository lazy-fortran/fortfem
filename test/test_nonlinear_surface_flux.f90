program test_nonlinear_surface_flux
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_nonlinear_surface_flux, &
        assemble_nonlinear_surface_flux_jvp, &
        assemble_nonlinear_surface_flux_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: quadrature_count = 2, dof_count = 3
    integer, parameter :: component_count = 2, tag_count = 2
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: trace_basis(quadrature_count, dof_count, component_count)
    real(dp) :: surface_weights(quadrature_count), surface_normals(3, quadrature_count)
    real(dp) :: trace_state(quadrature_count, component_count)
    integer :: surface_tags(quadrature_count)
    real(dp) :: residual(dof_count, component_count), ledger(tag_count, component_count)
    real(dp) :: residual_expected(dof_count, component_count)
    real(dp) :: ledger_expected(tag_count, component_count)
    real(dp) :: trace_basis_dot(quadrature_count, dof_count, component_count)
    real(dp) :: surface_weights_dot(quadrature_count)
    real(dp) :: surface_normals_dot(3, quadrature_count)
    real(dp) :: trace_state_dot(quadrature_count, component_count)
    real(dp) :: residual_dot(dof_count, component_count)
    real(dp) :: ledger_dot(tag_count, component_count)
    real(dp) :: residual_plus(dof_count, component_count)
    real(dp) :: ledger_plus(tag_count, component_count)
    real(dp) :: residual_minus(dof_count, component_count)
    real(dp) :: ledger_minus(tag_count, component_count)
    real(dp) :: residual_bar(dof_count, component_count)
    real(dp) :: ledger_bar(tag_count, component_count)
    real(dp) :: trace_basis_bar(quadrature_count, dof_count, component_count)
    real(dp) :: surface_weights_bar(quadrature_count)
    real(dp) :: surface_normals_bar(3, quadrature_count)
    real(dp) :: trace_state_bar(quadrature_count, component_count)
    real(dp) :: lhs, rhs
    integer :: status, status_plus, status_minus

    trace_basis = 0.0_dp
    trace_basis(1, 1, :) = [1.0_dp, 0.5_dp]
    trace_basis(1, 2, :) = [0.2_dp, -0.4_dp]
    trace_basis(1, 3, :) = [-0.3_dp, 0.7_dp]
    trace_basis(2, 1, :) = [0.1_dp, -0.2_dp]
    trace_basis(2, 2, :) = [0.8_dp, 0.6_dp]
    trace_basis(2, 3, :) = [0.4_dp, -0.1_dp]
    surface_weights = [2.0_dp, 3.0_dp]
    surface_normals = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp], [3, quadrature_count])
    trace_state(1, :) = [0.7_dp, -0.2_dp]
    trace_state(2, :) = [-0.4_dp, 0.9_dp]
    surface_tags = [1, 2]
    trace_basis_dot = 0.0_dp
    trace_basis_dot(1, 1, :) = [-0.02_dp, 0.04_dp]
    trace_basis_dot(1, 2, :) = [0.03_dp, -0.01_dp]
    trace_basis_dot(1, 3, :) = [0.01_dp, 0.02_dp]
    trace_basis_dot(2, 1, :) = [0.05_dp, -0.03_dp]
    trace_basis_dot(2, 2, :) = [-0.04_dp, 0.02_dp]
    trace_basis_dot(2, 3, :) = [0.02_dp, 0.01_dp]
    surface_weights_dot = [0.1_dp, -0.2_dp]
    surface_normals_dot = reshape([0.03_dp, -0.02_dp, 0.04_dp, &
        -0.01_dp, 0.05_dp, -0.03_dp], [3, quadrature_count])
    trace_state_dot(1, :) = [0.06_dp, -0.05_dp]
    trace_state_dot(2, :) = [-0.02_dp, 0.07_dp]
    residual_bar = reshape([0.2_dp, -0.3_dp, 0.4_dp, 0.5_dp, &
        -0.6_dp, 0.7_dp], [dof_count, component_count])
    ledger_bar = reshape([0.8_dp, -0.9_dp, 0.1_dp, 0.2_dp], &
        [tag_count, component_count])

    call assemble_nonlinear_surface_flux( &
        trace_basis, surface_weights, surface_normals, surface_tags, trace_state, &
        residual, ledger, nonlinear_flux_value, status)
    call independent_value_oracle( &
        trace_basis, surface_weights, surface_normals, surface_tags, trace_state, &
        residual_expected, ledger_expected)
    call check_condition(status == 0, &
        "nonlinear surface flux value accepts tagged traces")
    call check_condition(maxval(abs(residual - residual_expected)) < 1.0e-14_dp .and. &
        maxval(abs(ledger - ledger_expected)) < 1.0e-14_dp, &
        "nonlinear surface flux value matches an independent ledger oracle")

    call assemble_nonlinear_surface_flux_jvp( &
        trace_basis, surface_weights, surface_normals, surface_tags, trace_state, &
        trace_basis_dot, surface_weights_dot, surface_normals_dot, trace_state_dot, &
        residual_dot, ledger_dot, nonlinear_flux_value, nonlinear_flux_jvp, status)
    call assemble_nonlinear_surface_flux( &
        trace_basis + step*trace_basis_dot, &
        surface_weights + step*surface_weights_dot, &
        surface_normals + step*surface_normals_dot, surface_tags, &
        trace_state + step*trace_state_dot, residual_plus, ledger_plus, &
        nonlinear_flux_value, status_plus)
    call assemble_nonlinear_surface_flux( &
        trace_basis - step*trace_basis_dot, &
        surface_weights - step*surface_weights_dot, &
        surface_normals - step*surface_normals_dot, surface_tags, &
        trace_state - step*trace_state_dot, residual_minus, ledger_minus, &
        nonlinear_flux_value, status_minus)
    call check_condition(status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "nonlinear surface flux JVP accepts fixed tags")
    call check_condition(maxval(abs(residual_dot - (residual_plus - residual_minus)/ &
        (2.0_dp*step))) < 5.0e-9_dp .and. maxval(abs(ledger_dot - &
        (ledger_plus - ledger_minus)/(2.0_dp*step))) < 5.0e-9_dp, &
        "nonlinear surface flux JVP matches central differences")

    call assemble_nonlinear_surface_flux_vjp( &
        trace_basis, surface_weights, surface_normals, surface_tags, trace_state, &
        residual_bar, ledger_bar, trace_basis_bar, surface_weights_bar, &
        surface_normals_bar, trace_state_bar, nonlinear_flux_value, &
        nonlinear_flux_vjp, status)
    lhs = sum(residual_bar*residual_dot) + sum(ledger_bar*ledger_dot)
    rhs = sum(trace_basis_bar*trace_basis_dot) + &
        sum(surface_weights_bar*surface_weights_dot) + &
        sum(surface_normals_bar*surface_normals_dot) + &
        sum(trace_state_bar*trace_state_dot)
    call check_condition(status == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "nonlinear surface flux VJP satisfies the full dot-product identity")

    surface_weights(1) = 0.0_dp
    call assemble_nonlinear_surface_flux( &
        trace_basis, surface_weights, surface_normals, surface_tags, trace_state, &
        residual, ledger, nonlinear_flux_value, status)
    call check_condition(status /= 0, &
        "nonlinear surface flux rejects non-positive weights")
    surface_weights(1) = 2.0_dp
    surface_tags(1) = tag_count + 1
    call assemble_nonlinear_surface_flux( &
        trace_basis, surface_weights, surface_normals, surface_tags, trace_state, &
        residual, ledger, nonlinear_flux_value, status)
    call check_condition(status /= 0, "nonlinear surface flux rejects missing tags")
    call check_summary("Nonlinear material-surface flux")

contains

    subroutine nonlinear_flux_value(state, normal, tag, flux, status)
        real(dp), intent(in) :: state(:), normal(3)
        integer, intent(in) :: tag
        real(dp), intent(out) :: flux(:)
        integer, intent(out) :: status

        real(dp), parameter :: coefficient(2, 2) = reshape([ &
            1.2_dp, 0.8_dp, 0.6_dp, 1.5_dp], [2, 2])

        flux = 0.0_dp
        status = 1
        if (size(state) /= component_count .or. size(flux) /= component_count .or. &
            tag < 1 .or. tag > tag_count) return
        flux(1) = coefficient(1, tag)*state(1)**2 + 0.1_dp*normal(1) + &
            0.2_dp*normal(2)
        flux(2) = coefficient(2, tag)*state(2)**2 - 0.3_dp*normal(2) + &
            0.05_dp*normal(3)
        status = 0
    end subroutine nonlinear_flux_value

    subroutine nonlinear_flux_jvp(state, normal, tag, state_dot, normal_dot, &
            flux_dot, status)
        real(dp), intent(in) :: state(:), normal(3), state_dot(:), normal_dot(3)
        integer, intent(in) :: tag
        real(dp), intent(out) :: flux_dot(:)
        integer, intent(out) :: status

        real(dp), parameter :: coefficient(2, 2) = reshape([ &
            1.2_dp, 0.8_dp, 0.6_dp, 1.5_dp], [2, 2])

        flux_dot = 0.0_dp
        status = 1
        if (size(state) /= component_count .or. &
            size(state_dot) /= component_count .or. &
            size(flux_dot) /= component_count .or. tag < 1 .or. tag > tag_count) return
        flux_dot(1) = 2.0_dp*coefficient(1, tag)*state(1)*state_dot(1) + &
            0.1_dp*normal_dot(1) + 0.2_dp*normal_dot(2)
        flux_dot(2) = 2.0_dp*coefficient(2, tag)*state(2)*state_dot(2) - &
            0.3_dp*normal_dot(2) + 0.05_dp*normal_dot(3)
        status = 0
    end subroutine nonlinear_flux_jvp

    subroutine nonlinear_flux_vjp(state, normal, tag, flux_bar, state_bar, &
            normal_bar, status)
        real(dp), intent(in) :: state(:), normal(3), flux_bar(:)
        integer, intent(in) :: tag
        real(dp), intent(out) :: state_bar(:), normal_bar(3)
        integer, intent(out) :: status

        real(dp), parameter :: coefficient(2, 2) = reshape([ &
            1.2_dp, 0.8_dp, 0.6_dp, 1.5_dp], [2, 2])

        state_bar = 0.0_dp
        normal_bar = 0.0_dp
        status = 1
        if (size(state) /= component_count .or. size(flux_bar) /= component_count .or. &
            size(state_bar) /= component_count .or. tag < 1 .or. tag > tag_count) return
        state_bar(1) = flux_bar(1)*2.0_dp*coefficient(1, tag)*state(1)
        state_bar(2) = flux_bar(2)*2.0_dp*coefficient(2, tag)*state(2)
        normal_bar(1) = 0.1_dp*flux_bar(1)
        normal_bar(2) = 0.2_dp*flux_bar(1) - 0.3_dp*flux_bar(2)
        normal_bar(3) = 0.05_dp*flux_bar(2)
        status = 0
    end subroutine nonlinear_flux_vjp

    subroutine independent_value_oracle( &
            basis, weights, normals, tags, state, residual_oracle, ledger_oracle)
        real(dp), intent(in) :: basis(:, :, :), weights(:), normals(:, :), state(:, :)
        integer, intent(in) :: tags(:)
        real(dp), intent(out) :: residual_oracle(:, :), ledger_oracle(:, :)

        real(dp) :: flux(component_count)
        integer :: quadrature, dof, component, local_status

        residual_oracle = 0.0_dp
        ledger_oracle = 0.0_dp
        do quadrature = 1, quadrature_count
            call nonlinear_flux_value( &
                state(quadrature, :), normals(:, quadrature), tags(quadrature), &
                flux, local_status)
            do component = 1, component_count
                ledger_oracle(tags(quadrature), component) = &
                    ledger_oracle(tags(quadrature), component) + &
                    weights(quadrature)*flux(component)
                do dof = 1, dof_count
                    residual_oracle(dof, component) = &
                        residual_oracle(dof, component) + &
                        basis(quadrature, dof, component)*weights(quadrature)* &
                        flux(component)
                end do
            end do
        end do
    end subroutine independent_value_oracle

end program test_nonlinear_surface_flux
