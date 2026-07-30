module fortfem_planar_maxwell_dtn
    !! Fourier Maxwell capacity operator on a biperiodic planar boundary.
    !!
    !! The modal transparent condition follows Jiang et al.,
    !! arXiv:1811.12449. For tangential wavevector xi and outgoing normal
    !! wavenumber beta=sqrt(k^2-|xi|^2), this routine returns
    !! n x curl(E) = -i [beta I + xi xi^T/beta] E_t.
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c
    implicit none
    private

    public :: apply_planar_maxwell_dtn
    public :: assemble_planar_maxwell_dtn_form

contains

    subroutine assemble_planar_maxwell_dtn_form( &
            nx, ny, wave_number, length_x, length_y, form, status)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: wave_number, length_x, length_y
        complex(dp), allocatable, intent(out) :: form(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: derivative(:, :, :), trace(:, :, :)
        real(dp) :: weight
        integer :: column, component, i, j, quotient

        status = 1
        if (allocated(form)) deallocate(form)
        if (nx < 2 .or. ny < 2) return
        if (modulo(nx, 2) == 0 .or. modulo(ny, 2) == 0) return
        if (wave_number <= 0.0_dp) return
        if (length_x <= 0.0_dp .or. length_y <= 0.0_dp) return
        allocate(form(2*nx*ny, 2*nx*ny))
        allocate(trace(2, nx, ny), derivative(2, nx, ny))
        weight = length_x*length_y/real(nx*ny, dp)
        do column = 1, size(form, 2)
            trace = cmplx(0.0_dp, 0.0_dp, dp)
            component = modulo(column - 1, 2) + 1
            quotient = (column - 1)/2
            i = modulo(quotient, nx) + 1
            j = quotient/nx + 1
            trace(component, i, j) = cmplx(1.0_dp, 0.0_dp, dp)
            call apply_planar_maxwell_dtn( &
                trace, wave_number, length_x, length_y, derivative, status)
            if (status /= 0) return
            form(:, column) = weight*reshape(derivative, [size(derivative)])
        end do
        status = 0
    end subroutine assemble_planar_maxwell_dtn_form

    subroutine apply_planar_maxwell_dtn( &
            tangential_trace, wave_number, length_x, length_y, &
            tangential_derivative, status)
        complex(dp), intent(in) :: tangential_trace(:, :, :)
        real(dp), intent(in) :: wave_number, length_x, length_y
        complex(dp), intent(out) :: tangential_derivative(:, :, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: line(:), spectrum(:, :, :)
        complex(dp) :: beta, coefficient(2), imaginary_unit
        real(dp) :: radicand, xi(2)
        integer :: component, i, j, mode_x, mode_y, nx, ny

        status = 1
        tangential_derivative = cmplx(0.0_dp, 0.0_dp, dp)
        if (size(tangential_trace, 1) /= 2) return
        nx = size(tangential_trace, 2)
        ny = size(tangential_trace, 3)
        if (nx < 2 .or. ny < 2) return
        if (size(tangential_derivative, 1) /= 2) return
        if (size(tangential_derivative, 2) /= nx) return
        if (size(tangential_derivative, 3) /= ny) return
        if (wave_number <= 0.0_dp) return
        if (length_x <= 0.0_dp .or. length_y <= 0.0_dp) return

        allocate(spectrum(2, nx, ny))
        spectrum = tangential_trace
        do component = 1, 2
            allocate(line(nx))
            do j = 1, ny
                line = spectrum(component, :, j)
                call fft_c2c(line, -1)
                spectrum(component, :, j) = line
            end do
            deallocate(line)
            allocate(line(ny))
            do i = 1, nx
                line = spectrum(component, i, :)
                call fft_c2c(line, -1)
                spectrum(component, i, :) = line
            end do
            deallocate(line)
        end do

        imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)
        do j = 1, ny
            mode_y = j - 1
            if (mode_y > ny/2) mode_y = mode_y - ny
            xi(2) = 2.0_dp*acos(-1.0_dp)*real(mode_y, dp)/length_y
            do i = 1, nx
                mode_x = i - 1
                if (mode_x > nx/2) mode_x = mode_x - nx
                xi(1) = 2.0_dp*acos(-1.0_dp)*real(mode_x, dp)/length_x
                radicand = wave_number**2 - sum(xi**2)
                if (radicand >= 0.0_dp) then
                    beta = cmplx(sqrt(radicand), 0.0_dp, dp)
                else
                    beta = cmplx(0.0_dp, sqrt(-radicand), dp)
                end if
                if (abs(beta) <= 64.0_dp*epsilon(1.0_dp)*wave_number) then
                    if (maxval(abs(spectrum(:, i, j))) > &
                        1.0e3_dp*epsilon(1.0_dp)) return
                    spectrum(:, i, j) = cmplx(0.0_dp, 0.0_dp, dp)
                    cycle
                end if
                coefficient = spectrum(:, i, j)
                spectrum(:, i, j) = -imaginary_unit*( &
                    beta*coefficient + xi*sum(xi*coefficient)/beta)
            end do
        end do

        do component = 1, 2
            allocate(line(ny))
            do i = 1, nx
                line = spectrum(component, i, :)
                call fft_c2c(line, 1)
                spectrum(component, i, :) = line
            end do
            deallocate(line)
            allocate(line(nx))
            do j = 1, ny
                line = spectrum(component, :, j)
                call fft_c2c(line, 1)
                spectrum(component, :, j) = line
            end do
            deallocate(line)
        end do
        tangential_derivative = spectrum/real(nx*ny, dp)
        status = 0
    end subroutine apply_planar_maxwell_dtn

end module fortfem_planar_maxwell_dtn
