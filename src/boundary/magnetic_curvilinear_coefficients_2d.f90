module fortfem_magnetic_curvilinear_coefficients_2d
    use fortfem_generated_magnetic_curvilinear_coefficients_2d, only: &
        generated_magnetic_curvilinear_coefficients_2d
    use fortfem_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    private

    public :: scalar_reluctivity_curvilinear_fourier_coefficients

contains

    pure subroutine scalar_reluctivity_curvilinear_fourier_coefficients( &
            metric, reluctivity, curl_weight, mass_tensor, status)
        real(dp), intent(in) :: metric(3, 3), reluctivity
        real(dp), intent(out) :: curl_weight, mass_tensor(2, 2)
        integer, intent(out) :: status

        real(dp) :: determinant, scale, tolerance

        curl_weight = 0.0_dp
        mass_tensor = 0.0_dp
        status = 1
        if (.not. all(ieee_is_finite(metric))) return
        if (.not. ieee_is_finite(reluctivity)) return
        if (reluctivity < 0.0_dp) return
        scale = max(1.0_dp, maxval(abs(metric)))
        tolerance = 128.0_dp * epsilon(1.0_dp) * scale
        if (maxval(abs(metric - transpose(metric))) > tolerance) return
        if (max(abs(metric(1, 3)), abs(metric(2, 3))) > tolerance) return
        if (metric(1, 1) <= tolerance) return
        determinant = &
            metric(1, 1) * metric(2, 2) - metric(1, 2)**2
        if (determinant <= tolerance * scale) return
        if (metric(3, 3) <= tolerance) return

        call generated_magnetic_curvilinear_coefficients_2d( &
            metric(1, 1), metric(1, 2), metric(2, 2), metric(3, 3), &
            reluctivity, curl_weight, mass_tensor)
        status = 0
    end subroutine scalar_reluctivity_curvilinear_fourier_coefficients

end module fortfem_magnetic_curvilinear_coefficients_2d
