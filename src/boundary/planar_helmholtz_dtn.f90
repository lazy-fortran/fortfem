module fortfem_planar_helmholtz_dtn
    use fortfem_kinds, only: dp
    use fortfem_generated_planar_helmholtz_evanescent_jvp, only: &
        generated_planar_helmholtz_evanescent_jvp
    use fortfem_generated_planar_helmholtz_evanescent_vjp, only: &
        generated_planar_helmholtz_evanescent_vjp
    use fortfem_generated_planar_helmholtz_propagating_jvp, only: &
        generated_planar_helmholtz_propagating_jvp
    use fortfem_generated_planar_helmholtz_propagating_vjp, only: &
        generated_planar_helmholtz_propagating_vjp
    use fortnum_fft, only: fft_c2c, fft_c2c_vjp
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: apply_planar_helmholtz_dtn
    public :: apply_planar_helmholtz_dtn_jvp
    public :: apply_planar_helmholtz_dtn_vjp
    public :: assemble_planar_helmholtz_dtn_form
    public :: assemble_planar_helmholtz_dtn_form_jvp
    public :: assemble_planar_helmholtz_dtn_form_vjp

contains

    subroutine assemble_planar_helmholtz_dtn_form( &
            sample_count, wavenumber, period, form, status)
        integer, intent(in) :: sample_count
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(out) :: form(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: derivative(:), trace(:)
        complex(dp), allocatable :: strong_operator(:, :)
        real(dp) :: spacing
        integer :: column, next, previous, row

        form = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (sample_count < 3 .or. wavenumber < 0.0_dp .or. &
            period <= 0.0_dp) return
        if (size(form, 1) /= sample_count .or. &
            size(form, 2) /= sample_count) return

        allocate(trace(sample_count), derivative(sample_count))
        allocate(strong_operator(sample_count, sample_count))
        do column = 1, sample_count
            trace = cmplx(0.0_dp, 0.0_dp, dp)
            trace(column) = cmplx(1.0_dp, 0.0_dp, dp)
            call apply_planar_helmholtz_dtn( &
                trace, wavenumber, period, derivative)
            strong_operator(:, column) = derivative
        end do

        spacing = period/real(sample_count, dp)
        do row = 1, sample_count
            previous = modulo(row - 2, sample_count) + 1
            next = modulo(row, sample_count) + 1
            form(row, :) = spacing/6.0_dp*( &
                strong_operator(previous, :) + &
                4.0_dp*strong_operator(row, :) + &
                strong_operator(next, :))
        end do
        status = 0
    end subroutine assemble_planar_helmholtz_dtn_form

    subroutine assemble_planar_helmholtz_dtn_form_jvp( &
            sample_count, wavenumber, period, wavenumber_dot, period_dot, &
            form_dot, status)
        integer, intent(in) :: sample_count
        real(dp), intent(in) :: wavenumber, period
        real(dp), intent(in) :: wavenumber_dot, period_dot
        complex(dp), intent(out) :: form_dot(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: derivative(:), derivative_dot(:)
        complex(dp), allocatable :: strong_operator(:, :)
        complex(dp), allocatable :: strong_operator_dot(:, :), trace(:)
        complex(dp), allocatable :: trace_dot(:)
        real(dp) :: spacing, spacing_dot
        integer :: column, next, previous, row

        form_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (sample_count < 3 .or. wavenumber < 0.0_dp .or. &
            period <= 0.0_dp) return
        if (size(form_dot, 1) /= sample_count .or. &
            size(form_dot, 2) /= sample_count) return

        allocate (trace(sample_count), trace_dot(sample_count))
        allocate (derivative(sample_count), derivative_dot(sample_count))
        allocate (strong_operator(sample_count, sample_count))
        allocate (strong_operator_dot(sample_count, sample_count))
        trace_dot = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, sample_count
            trace = cmplx(0.0_dp, 0.0_dp, dp)
            trace(column) = cmplx(1.0_dp, 0.0_dp, dp)
            call apply_planar_helmholtz_dtn( &
                trace, wavenumber, period, derivative)
            call apply_planar_helmholtz_dtn_jvp( &
                trace, wavenumber, period, trace_dot, wavenumber_dot, &
                period_dot, derivative_dot)
            strong_operator(:, column) = derivative
            strong_operator_dot(:, column) = derivative_dot
        end do

        spacing = period/real(sample_count, dp)
        spacing_dot = period_dot/real(sample_count, dp)
        do row = 1, sample_count
            previous = modulo(row - 2, sample_count) + 1
            next = modulo(row, sample_count) + 1
            form_dot(row, :) = spacing_dot/6.0_dp*( &
                strong_operator(previous, :) + &
                4.0_dp*strong_operator(row, :) + &
                strong_operator(next, :)) + spacing/6.0_dp*( &
                strong_operator_dot(previous, :) + &
                4.0_dp*strong_operator_dot(row, :) + &
                strong_operator_dot(next, :))
        end do
        status = 0
    end subroutine assemble_planar_helmholtz_dtn_form_jvp

    subroutine assemble_planar_helmholtz_dtn_form_vjp( &
            sample_count, wavenumber, period, form_bar, wavenumber_bar, &
            period_bar, status)
        integer, intent(in) :: sample_count
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: form_bar(:, :)
        real(dp), intent(out) :: wavenumber_bar, period_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: derivative(:), strong_bar(:, :)
        complex(dp), allocatable :: strong_operator(:, :), trace(:), trace_bar(:)
        real(dp) :: local_period_bar, local_wavenumber_bar, spacing
        integer :: column, next, previous, row

        wavenumber_bar = 0.0_dp
        period_bar = 0.0_dp
        status = 1
        if (sample_count < 3 .or. wavenumber < 0.0_dp .or. &
            period <= 0.0_dp) return
        if (size(form_bar, 1) /= sample_count .or. &
            size(form_bar, 2) /= sample_count) return

        allocate (trace(sample_count), trace_bar(sample_count))
        allocate (derivative(sample_count))
        allocate (strong_operator(sample_count, sample_count))
        allocate (strong_bar(sample_count, sample_count))
        do column = 1, sample_count
            trace = cmplx(0.0_dp, 0.0_dp, dp)
            trace(column) = cmplx(1.0_dp, 0.0_dp, dp)
            call apply_planar_helmholtz_dtn( &
                trace, wavenumber, period, derivative)
            strong_operator(:, column) = derivative
        end do

        spacing = period/real(sample_count, dp)
        strong_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do row = 1, sample_count
            previous = modulo(row - 2, sample_count) + 1
            next = modulo(row, sample_count) + 1
            strong_bar(previous, :) = strong_bar(previous, :) + &
                spacing/6.0_dp*form_bar(row, :)
            strong_bar(row, :) = strong_bar(row, :) + &
                4.0_dp*spacing/6.0_dp*form_bar(row, :)
            strong_bar(next, :) = strong_bar(next, :) + &
                spacing/6.0_dp*form_bar(row, :)
            period_bar = period_bar + real(sum(conjg(form_bar(row, :))*( &
                strong_operator(previous, :) + &
                4.0_dp*strong_operator(row, :) + &
                strong_operator(next, :)))/ &
                (6.0_dp*real(sample_count, dp)), dp)
        end do

        do column = 1, sample_count
            trace = cmplx(0.0_dp, 0.0_dp, dp)
            trace(column) = cmplx(1.0_dp, 0.0_dp, dp)
            call apply_planar_helmholtz_dtn_vjp( &
                trace, wavenumber, period, strong_bar(:, column), trace_bar, &
                local_wavenumber_bar, local_period_bar)
            wavenumber_bar = wavenumber_bar + local_wavenumber_bar
            period_bar = period_bar + local_period_bar
        end do
        status = 0
    end subroutine assemble_planar_helmholtz_dtn_form_vjp

    subroutine apply_planar_helmholtz_dtn( &
            trace, wavenumber, period, normal_derivative)
        complex(dp), intent(in) :: trace(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(out) :: normal_derivative(:)

        complex(dp), allocatable :: spectrum(:)
        complex(dp) :: dtn_eigenvalue
        real(dp) :: tangential_wavenumber, radicand
        integer :: n, bin, mode

        n = size(trace)
        if (n < 1) error stop "Planar DtN requires at least one trace point"
        if (size(normal_derivative) /= n) then
            error stop "Planar DtN input and output sizes differ"
        end if
        if (wavenumber < 0.0_dp) then
            error stop "Planar DtN wavenumber must be nonnegative"
        end if
        if (period <= 0.0_dp) then
            error stop "Planar DtN period must be positive"
        end if

        allocate(spectrum(n))
        spectrum = trace
        call fft_c2c(spectrum, -1)

        do bin = 0, n - 1
            if (bin <= n / 2) then
                mode = bin
            else
                mode = bin - n
            end if
            tangential_wavenumber = 2.0_dp * pi * real(mode, dp) / period
            radicand = wavenumber**2 - tangential_wavenumber**2
            if (radicand >= 0.0_dp) then
                dtn_eigenvalue = &
                    cmplx(0.0_dp, sqrt(radicand), kind=dp)
            else
                dtn_eigenvalue = &
                    cmplx(-sqrt(-radicand), 0.0_dp, kind=dp)
            end if
            spectrum(bin + 1) = dtn_eigenvalue * spectrum(bin + 1)
        end do

        call fft_c2c(spectrum, 1)
        normal_derivative = spectrum / real(n, dp)
    end subroutine apply_planar_helmholtz_dtn

    subroutine apply_planar_helmholtz_dtn_jvp( &
            trace, wavenumber, period, trace_dot, wavenumber_dot, period_dot, &
            normal_derivative_dot)
        complex(dp), intent(in) :: trace(:), trace_dot(:)
        real(dp), intent(in) :: wavenumber, period
        real(dp), intent(in) :: wavenumber_dot, period_dot
        complex(dp), intent(out) :: normal_derivative_dot(:)

        complex(dp), allocatable :: spectrum(:), spectrum_dot(:)
        complex(dp) :: dtn_eigenvalue, dtn_eigenvalue_dot
        real(dp) :: beta, beta_dot, cutoff_scale, radicand
        real(dp) :: tangential_wavenumber, tangential_wavenumber_dot
        integer :: bin, mode, n

        n = size(trace)
        call validate_derivative_inputs( &
            n, size(trace_dot), size(normal_derivative_dot), wavenumber, period)
        allocate (spectrum(n), spectrum_dot(n))
        spectrum = trace
        spectrum_dot = trace_dot
        call fft_c2c(spectrum, -1)
        call fft_c2c(spectrum_dot, -1)
        do bin = 0, n - 1
            mode = signed_mode(bin, n)
            tangential_wavenumber = 2.0_dp*pi*real(mode, dp)/period
            tangential_wavenumber_dot = &
                -tangential_wavenumber*period_dot/period
            radicand = wavenumber**2 - tangential_wavenumber**2
            cutoff_scale = max( &
                1.0_dp, wavenumber**2, tangential_wavenumber**2)
            if (abs(radicand) <= epsilon(1.0_dp)*cutoff_scale) then
                error stop "Planar DtN parameter derivative is singular at cutoff"
            end if
            if (radicand > 0.0_dp) then
                beta = sqrt(radicand)
                call generated_planar_helmholtz_propagating_jvp( &
                    wavenumber, tangential_wavenumber, wavenumber_dot, &
                    tangential_wavenumber_dot, beta_dot)
                dtn_eigenvalue = cmplx(0.0_dp, beta, dp)
                dtn_eigenvalue_dot = cmplx(0.0_dp, beta_dot, dp)
            else
                beta = sqrt(-radicand)
                call generated_planar_helmholtz_evanescent_jvp( &
                    wavenumber, tangential_wavenumber, wavenumber_dot, &
                    tangential_wavenumber_dot, beta_dot)
                dtn_eigenvalue = cmplx(-beta, 0.0_dp, dp)
                dtn_eigenvalue_dot = cmplx(-beta_dot, 0.0_dp, dp)
            end if
            spectrum_dot(bin + 1) = &
                dtn_eigenvalue*spectrum_dot(bin + 1) + &
                dtn_eigenvalue_dot*spectrum(bin + 1)
        end do
        call fft_c2c(spectrum_dot, 1)
        normal_derivative_dot = spectrum_dot/real(n, dp)
    end subroutine apply_planar_helmholtz_dtn_jvp

    subroutine apply_planar_helmholtz_dtn_vjp( &
            trace, wavenumber, period, normal_derivative_bar, trace_bar, &
            wavenumber_bar, period_bar)
        complex(dp), intent(in) :: trace(:), normal_derivative_bar(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(out) :: trace_bar(:)
        real(dp), intent(out) :: wavenumber_bar, period_bar

        complex(dp), allocatable :: spectrum(:), spectrum_bar(:)
        complex(dp) :: dtn_eigenvalue, modal_derivative
        real(dp) :: beta, beta_bar, cutoff_scale, local_tangential_bar
        real(dp) :: local_wavenumber_bar, radicand, tangential_wavenumber
        integer :: bin, mode, n

        n = size(trace)
        call validate_derivative_inputs( &
            n, size(normal_derivative_bar), size(trace_bar), wavenumber, period)
        allocate (spectrum(n), spectrum_bar(n))
        spectrum = trace
        call fft_c2c(spectrum, -1)
        spectrum_bar = normal_derivative_bar/real(n, dp)
        call fft_c2c_vjp(spectrum_bar, 1)
        wavenumber_bar = 0.0_dp
        period_bar = 0.0_dp
        do bin = 0, n - 1
            mode = signed_mode(bin, n)
            tangential_wavenumber = 2.0_dp*pi*real(mode, dp)/period
            radicand = wavenumber**2 - tangential_wavenumber**2
            cutoff_scale = max( &
                1.0_dp, wavenumber**2, tangential_wavenumber**2)
            if (abs(radicand) <= epsilon(1.0_dp)*cutoff_scale) then
                error stop "Planar DtN parameter derivative is singular at cutoff"
            end if
            if (radicand > 0.0_dp) then
                beta = sqrt(radicand)
                dtn_eigenvalue = cmplx(0.0_dp, beta, dp)
                modal_derivative = &
                    cmplx(0.0_dp, 1.0_dp, dp)*spectrum(bin + 1)
                beta_bar = real( &
                    conjg(spectrum_bar(bin + 1))*modal_derivative, dp)
                call generated_planar_helmholtz_propagating_vjp( &
                    wavenumber, tangential_wavenumber, beta_bar, &
                    local_wavenumber_bar, local_tangential_bar)
            else
                beta = sqrt(-radicand)
                dtn_eigenvalue = cmplx(-beta, 0.0_dp, dp)
                modal_derivative = -spectrum(bin + 1)
                beta_bar = real( &
                    conjg(spectrum_bar(bin + 1))*modal_derivative, dp)
                call generated_planar_helmholtz_evanescent_vjp( &
                    wavenumber, tangential_wavenumber, beta_bar, &
                    local_wavenumber_bar, local_tangential_bar)
            end if
            wavenumber_bar = wavenumber_bar + local_wavenumber_bar
            period_bar = period_bar - &
                local_tangential_bar*tangential_wavenumber/period
            spectrum_bar(bin + 1) = &
                conjg(dtn_eigenvalue)*spectrum_bar(bin + 1)
        end do
        call fft_c2c_vjp(spectrum_bar, -1)
        trace_bar = spectrum_bar
    end subroutine apply_planar_helmholtz_dtn_vjp

    pure integer function signed_mode(bin, sample_count) result(mode)
        integer, intent(in) :: bin, sample_count

        if (bin <= sample_count/2) then
            mode = bin
        else
            mode = bin - sample_count
        end if
    end function signed_mode

    subroutine validate_derivative_inputs( &
            trace_size, direction_size, output_size, wavenumber, period)
        integer, intent(in) :: trace_size, direction_size, output_size
        real(dp), intent(in) :: wavenumber, period

        if (trace_size < 1) then
            error stop "Planar DtN requires at least one trace point"
        end if
        if (direction_size /= trace_size .or. output_size /= trace_size) then
            error stop "Planar DtN derivative sizes differ"
        end if
        if (wavenumber < 0.0_dp) then
            error stop "Planar DtN wavenumber must be nonnegative"
        end if
        if (period <= 0.0_dp) then
            error stop "Planar DtN period must be positive"
        end if
    end subroutine validate_derivative_inputs

end module fortfem_planar_helmholtz_dtn
