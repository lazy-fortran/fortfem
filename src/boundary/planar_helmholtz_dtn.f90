module fortfem_planar_helmholtz_dtn
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: apply_planar_helmholtz_dtn

contains

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
