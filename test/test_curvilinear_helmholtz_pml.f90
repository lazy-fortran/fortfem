program test_curvilinear_helmholtz_pml
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        curvilinear_curl_curl_pml_coefficients, &
        curvilinear_curl_curl_pml_coefficients_jvp, &
        curvilinear_curl_curl_pml_coefficients_vjp, &
        curvilinear_scalar_helmholtz_pml_coefficients, &
        curvilinear_scalar_helmholtz_pml_coefficients_jvp, &
        curvilinear_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_kinds, only: dp
    implicit none

    complex(dp) :: stretch(3, 3), stretch_dot(3, 3), singular_stretch(3, 3)
    complex(dp) :: gradient(3, 3), gradient_dot(3, 3)
    complex(dp) :: gradient_minus(3, 3), gradient_plus(3, 3), gradient_bar(3, 3)
    complex(dp) :: curl(3, 3), curl_dot(3, 3), curl_minus(3, 3), curl_plus(3, 3)
    complex(dp) :: curl_bar(3, 3), mass_tensor(3, 3), mass_tensor_dot(3, 3)
    complex(dp) :: mass_tensor_minus(3, 3), mass_tensor_plus(3, 3)
    complex(dp) :: mass_tensor_bar(3, 3), mass, mass_dot, mass_minus, mass_plus
    complex(dp) :: mass_bar, stretch_bar(3, 3)
    complex(dp) :: diagonal(3, 3), diagonal_gradient(3, 3), diagonal_curl(3, 3)
    complex(dp) :: diagonal_mass_tensor(3, 3)
    complex(dp) :: shear(3, 3), shear_gradient(3, 3), shear_curl(3, 3)
    complex(dp) :: shear_mass_tensor(3, 3), shear_gradient_expected(3, 3)
    complex(dp) :: shear_curl_expected(3, 3), shear_mass
    real(dp) :: forward_pairing, reverse_pairing, step
    integer :: status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    stretch = reshape([ &
        cmplx(1.2_dp, 0.1_dp, dp), cmplx(0.1_dp, -0.2_dp, dp), &
        cmplx(-0.2_dp, 0.05_dp, dp), cmplx(0.3_dp, 0.2_dp, dp), &
        cmplx(1.1_dp, 0.15_dp, dp), cmplx(0.08_dp, -0.1_dp, dp), &
        cmplx(0.04_dp, 0.03_dp, dp), cmplx(-0.15_dp, 0.12_dp, dp), &
        cmplx(0.9_dp, 0.25_dp, dp)], [3, 3])
    stretch_dot = reshape([ &
        cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(0.05_dp, 0.01_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.02_dp, -0.03_dp, dp), &
        cmplx(-0.01_dp, 0.03_dp, dp), cmplx(0.04_dp, -0.02_dp, dp), &
        cmplx(0.01_dp, 0.01_dp, dp)], [3, 3])

    call curvilinear_scalar_helmholtz_pml_coefficients( &
        stretch, gradient, mass, status)
    call record_condition(status == 0, "curvilinear scalar PML accepts a full stretch")
    call curvilinear_scalar_helmholtz_pml_coefficients_jvp( &
        stretch, stretch_dot, gradient_dot, mass_dot, status)
    call record_condition(status == 0, &
        "curvilinear scalar PML JVP accepts a full stretch")

    step = 1.0e-6_dp
    call curvilinear_scalar_helmholtz_pml_coefficients( &
        stretch - step*stretch_dot, gradient_minus, mass_minus, status_minus)
    call curvilinear_scalar_helmholtz_pml_coefficients( &
        stretch + step*stretch_dot, gradient_plus, mass_plus, status_plus)
    call record_condition(status_minus == 0 .and. status_plus == 0, &
        "curvilinear scalar PML central-difference states assemble")
    call record_condition(maxval(abs(gradient_dot - &
        (gradient_plus - gradient_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "curvilinear scalar PML gradient JVP matches central differences")
    call record_condition(abs(mass_dot - &
        (mass_plus - mass_minus)/(2.0_dp*step)) < 3.0e-9_dp, &
        "curvilinear scalar PML mass JVP matches central differences")

    gradient_bar = reshape([ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp), &
        cmplx(0.1_dp, 0.4_dp, dp), cmplx(0.4_dp, 0.1_dp, dp), &
        cmplx(-0.2_dp, -0.3_dp, dp), cmplx(0.3_dp, 0.2_dp, dp), &
        cmplx(-0.1_dp, 0.2_dp, dp), cmplx(0.1_dp, -0.4_dp, dp), &
        cmplx(0.2_dp, 0.3_dp, dp)], [3, 3])
    mass_bar = cmplx(-0.3_dp, 0.2_dp, dp)
    call curvilinear_scalar_helmholtz_pml_coefficients_vjp( &
        stretch, gradient_bar, mass_bar, stretch_bar, status)
    call record_condition(status == 0, &
        "curvilinear scalar PML VJP accepts a full stretch")
    forward_pairing = real(sum(conjg(gradient_bar)*gradient_dot) + &
        conjg(mass_bar)*mass_dot, dp)
    reverse_pairing = real(sum(conjg(stretch_bar)*stretch_dot), dp)
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-11_dp, &
        "curvilinear scalar PML JVP and VJP satisfy the real complex adjoint identity")

    call curvilinear_curl_curl_pml_coefficients( &
        stretch, curl, mass_tensor, status)
    call record_condition(status == 0, &
        "curvilinear curl-curl PML accepts a full stretch")
    call curvilinear_curl_curl_pml_coefficients_jvp( &
        stretch, stretch_dot, curl_dot, mass_tensor_dot, status)
    call record_condition(status == 0, &
        "curvilinear curl-curl PML JVP accepts a full stretch")
    call curvilinear_curl_curl_pml_coefficients( &
        stretch - step*stretch_dot, curl_minus, mass_tensor_minus, status_minus)
    call curvilinear_curl_curl_pml_coefficients( &
        stretch + step*stretch_dot, curl_plus, mass_tensor_plus, status_plus)
    call record_condition(maxval(abs(curl_dot - &
        (curl_plus - curl_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "curvilinear curl-curl coefficient JVP matches central differences")
    call record_condition(maxval(abs(mass_tensor_dot - &
        (mass_tensor_plus - mass_tensor_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "curvilinear Maxwell mass JVP matches central differences")

    curl_bar = reshape([ &
        cmplx(-0.2_dp, 0.3_dp, dp), cmplx(0.1_dp, 0.2_dp, dp), &
        cmplx(0.3_dp, -0.1_dp, dp), cmplx(0.2_dp, 0.1_dp, dp), &
        cmplx(-0.1_dp, 0.4_dp, dp), cmplx(0.2_dp, -0.3_dp, dp), &
        cmplx(0.1_dp, -0.2_dp, dp), cmplx(0.3_dp, 0.1_dp, dp), &
        cmplx(-0.2_dp, 0.2_dp, dp)], [3, 3])
    mass_tensor_bar = reshape([ &
        cmplx(0.1_dp, 0.2_dp, dp), cmplx(-0.2_dp, 0.1_dp, dp), &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(0.3_dp, 0.0_dp, dp), &
        cmplx(-0.1_dp, 0.2_dp, dp), cmplx(0.2_dp, 0.1_dp, dp), &
        cmplx(-0.3_dp, 0.1_dp, dp), cmplx(0.2_dp, -0.2_dp, dp), &
        cmplx(0.1_dp, 0.3_dp, dp)], [3, 3])
    call curvilinear_curl_curl_pml_coefficients_vjp( &
        stretch, curl_bar, mass_tensor_bar, stretch_bar, status)
    call record_condition(status == 0, &
        "curvilinear curl-curl PML VJP accepts a full stretch")
    forward_pairing = real(sum(conjg(curl_bar)*curl_dot) + &
        sum(conjg(mass_tensor_bar)*mass_tensor_dot), dp)
    reverse_pairing = real(sum(conjg(stretch_bar)*stretch_dot), dp)
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-11_dp, &
        "curvilinear Maxwell JVP and VJP satisfy the real complex adjoint identity")

    diagonal = cmplx(0.0_dp, 0.0_dp, dp)
    diagonal(1, 1) = cmplx(1.2_dp, 0.1_dp, dp)
    diagonal(2, 2) = cmplx(0.9_dp, 0.2_dp, dp)
    diagonal(3, 3) = cmplx(1.1_dp, 0.3_dp, dp)
    call curvilinear_scalar_helmholtz_pml_coefficients( &
        diagonal, diagonal_gradient, mass, status)
    call record_condition(status == 0, "diagonal stretch uses curvilinear scalar path")
    call curvilinear_curl_curl_pml_coefficients( &
        diagonal, diagonal_curl, diagonal_mass_tensor, status)
    call record_condition(status == 0, "diagonal stretch uses curvilinear Maxwell path")
    call record_condition(maxval(abs(diagonal_gradient - &
        diagonal_mass_tensor)) < 2.0e-14_dp, &
        "scalar and Maxwell mass tensors agree for diagonal stretches")
    call record_condition(abs(diagonal_gradient(1, 1) - &
        diagonal(2, 2)*diagonal(3, 3)/diagonal(1, 1)) < 2.0e-14_dp, &
        "curvilinear scalar tensor reduces to Cartesian diagonal form")

    shear = cmplx(0.0_dp, 0.0_dp, dp)
    shear(1, 1) = cmplx(2.0_dp, 0.0_dp, dp)
    shear(1, 2) = cmplx(1.0_dp, 0.0_dp, dp)
    shear(2, 2) = cmplx(3.0_dp, 0.0_dp, dp)
    shear(3, 3) = cmplx(4.0_dp, 0.0_dp, dp)
    shear_gradient_expected = cmplx(0.0_dp, 0.0_dp, dp)
    shear_gradient_expected(1, 1) = cmplx(20.0_dp/3.0_dp, 0.0_dp, dp)
    shear_gradient_expected(1, 2) = cmplx(-4.0_dp/3.0_dp, 0.0_dp, dp)
    shear_gradient_expected(2, 1) = cmplx(-4.0_dp/3.0_dp, 0.0_dp, dp)
    shear_gradient_expected(2, 2) = cmplx(8.0_dp/3.0_dp, 0.0_dp, dp)
    shear_gradient_expected(3, 3) = cmplx(3.0_dp/2.0_dp, 0.0_dp, dp)
    shear_curl_expected = cmplx(0.0_dp, 0.0_dp, dp)
    shear_curl_expected(1, 1) = cmplx(1.0_dp/6.0_dp, 0.0_dp, dp)
    shear_curl_expected(1, 2) = cmplx(1.0_dp/12.0_dp, 0.0_dp, dp)
    shear_curl_expected(2, 1) = cmplx(1.0_dp/12.0_dp, 0.0_dp, dp)
    shear_curl_expected(2, 2) = cmplx(5.0_dp/12.0_dp, 0.0_dp, dp)
    shear_curl_expected(3, 3) = cmplx(2.0_dp/3.0_dp, 0.0_dp, dp)
    call curvilinear_scalar_helmholtz_pml_coefficients( &
        shear, shear_gradient, shear_mass, status)
    call curvilinear_curl_curl_pml_coefficients( &
        shear, shear_curl, shear_mass_tensor, status)
    call record_condition(abs(shear_mass - cmplx(24.0_dp, 0.0_dp, dp)) < &
        2.0e-14_dp, "curvilinear stretch uses the independent triangular determinant")
    call record_condition(maxval(abs(shear_gradient - &
        shear_gradient_expected)) < 2.0e-14_dp, &
        "curvilinear scalar tensor retains off-diagonal metric coupling")
    call record_condition(maxval(abs(shear_curl - shear_curl_expected)) < &
        2.0e-14_dp, "curvilinear curl tensor matches the independent shear oracle")
    call record_condition(maxval(abs(shear_mass_tensor - &
        shear_gradient_expected)) < 2.0e-14_dp, &
        "curvilinear Maxwell mass tensor matches the shear oracle")

    singular_stretch = diagonal
    singular_stretch(3, :) = singular_stretch(2, :)
    call curvilinear_scalar_helmholtz_pml_coefficients( &
        singular_stretch, gradient, mass, status)
    call record_condition(status /= 0, "curvilinear PML rejects a singular stretch")

    call check_summary("Curvilinear Helmholtz and Maxwell PML")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_curvilinear_helmholtz_pml
