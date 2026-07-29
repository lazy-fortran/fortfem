module fortfem_planar_helmholtz_dtn
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: apply_planar_helmholtz_dtn
    public :: assemble_planar_helmholtz_dtn_form

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

end module fortfem_planar_helmholtz_dtn
