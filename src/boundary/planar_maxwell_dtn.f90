module fortfem_planar_maxwell_dtn
    !! Fourier Maxwell capacity operator on a biperiodic planar boundary.
    !!
    !! The modal transparent condition follows Jiang et al.,
    !! arXiv:1811.12449. For tangential wavevector xi and outgoing normal
    !! wavenumber beta=sqrt(k^2-|xi|^2), this routine returns
    !! n x curl(E) = -i [beta I + xi xi^T/beta] E_t.
    use fortfem_generated_planar_helmholtz_evanescent_jvp, only: &
        generated_planar_helmholtz_evanescent_jvp
    use fortfem_generated_planar_helmholtz_propagating_jvp, only: &
        generated_planar_helmholtz_propagating_jvp
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c
    implicit none
    private

    public :: apply_planar_maxwell_dtn
    public :: apply_planar_maxwell_dtn_jvp
    public :: apply_planar_maxwell_dtn_vjp
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

    subroutine apply_planar_maxwell_dtn_jvp( &
            tangential_trace, wave_number, length_x, length_y, &
            tangential_trace_dot, wave_number_dot, length_x_dot, length_y_dot, &
            tangential_derivative_dot, status)
        complex(dp), intent(in) :: tangential_trace(:, :, :)
        real(dp), intent(in) :: wave_number, length_x, length_y
        complex(dp), intent(in) :: tangential_trace_dot(:, :, :)
        real(dp), intent(in) :: wave_number_dot, length_x_dot, length_y_dot
        complex(dp), intent(out) :: tangential_derivative_dot(:, :, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:, :, :), spectrum_dot(:, :, :)
        complex(dp) :: matrix(2, 2), matrix_dot(2, 2)
        integer :: i, j, nx, ny

        tangential_derivative_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_product_inputs( &
            tangential_trace, tangential_trace_dot, tangential_derivative_dot, &
            wave_number, length_x, length_y, nx, ny, status)
        if (status /= 0) return
        allocate(spectrum(2, nx, ny), spectrum_dot(2, nx, ny))
        spectrum = tangential_trace
        spectrum_dot = tangential_trace_dot
        call transform_forward(spectrum)
        call transform_forward(spectrum_dot)
        do j = 1, ny
            do i = 1, nx
                call modal_matrix_product( &
                    i, j, nx, ny, wave_number, length_x, length_y, &
                    wave_number_dot, length_x_dot, length_y_dot, matrix, &
                    matrix_dot, status)
                if (status /= 0) return
                spectrum_dot(:, i, j) = &
                    matmul(matrix, spectrum_dot(:, i, j)) + &
                    matmul(matrix_dot, spectrum(:, i, j))
            end do
        end do
        call transform_inverse(spectrum_dot)
        tangential_derivative_dot = spectrum_dot/real(nx*ny, dp)
        status = 0
    end subroutine apply_planar_maxwell_dtn_jvp

    subroutine apply_planar_maxwell_dtn_vjp( &
            tangential_trace, wave_number, length_x, length_y, &
            tangential_derivative_bar, tangential_trace_bar, wave_number_bar, &
            length_x_bar, length_y_bar, status)
        complex(dp), intent(in) :: tangential_trace(:, :, :)
        real(dp), intent(in) :: wave_number, length_x, length_y
        complex(dp), intent(in) :: tangential_derivative_bar(:, :, :)
        complex(dp), intent(out) :: tangential_trace_bar(:, :, :)
        real(dp), intent(out) :: wave_number_bar, length_x_bar, length_y_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:, :, :), spectrum_bar(:, :, :)
        complex(dp) :: matrix(2, 2), matrix_dot(2, 2), product(2)
        real(dp) :: direction(3), sensitivity(3), scale
        integer :: active, i, j, nx, ny

        tangential_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        length_x_bar = 0.0_dp
        length_y_bar = 0.0_dp
        call validate_product_inputs( &
            tangential_trace, tangential_derivative_bar, tangential_trace_bar, &
            wave_number, length_x, length_y, nx, ny, status)
        if (status /= 0) return
        allocate(spectrum(2, nx, ny), spectrum_bar(2, nx, ny))
        spectrum = tangential_trace
        spectrum_bar = tangential_derivative_bar
        call transform_forward(spectrum)
        call transform_forward(spectrum_bar)
        scale = 1.0_dp/real(nx*ny, dp)
        do j = 1, ny
            do i = 1, nx
                call modal_matrix_product( &
                    i, j, nx, ny, wave_number, length_x, length_y, &
                    0.0_dp, 0.0_dp, 0.0_dp, matrix, matrix_dot, status)
                if (status /= 0) return
                sensitivity = 0.0_dp
                do active = 1, 3
                    direction = 0.0_dp
                    direction(active) = 1.0_dp
                    call modal_matrix_product( &
                        i, j, nx, ny, wave_number, length_x, length_y, &
                        direction(1), direction(2), direction(3), matrix, &
                        matrix_dot, status)
                    if (status /= 0) return
                    product = matmul(matrix_dot, spectrum(:, i, j))
                    sensitivity(active) = scale*real( &
                        sum(conjg(spectrum_bar(:, i, j))*product), dp)
                end do
                wave_number_bar = wave_number_bar + sensitivity(1)
                length_x_bar = length_x_bar + sensitivity(2)
                length_y_bar = length_y_bar + sensitivity(3)
                spectrum_bar(:, i, j) = &
                    matmul(conjg(transpose(matrix)), spectrum_bar(:, i, j))
            end do
        end do
        call transform_inverse(spectrum_bar)
        tangential_trace_bar = spectrum_bar/real(nx*ny, dp)
        status = 0
    end subroutine apply_planar_maxwell_dtn_vjp

    subroutine modal_matrix_product( &
            i, j, nx, ny, wave_number, length_x, length_y, wave_number_dot, &
            length_x_dot, length_y_dot, matrix, matrix_dot, status)
        integer, intent(in) :: i, j, nx, ny
        real(dp), intent(in) :: wave_number, length_x, length_y
        real(dp), intent(in) :: wave_number_dot, length_x_dot, length_y_dot
        complex(dp), intent(out) :: matrix(2, 2), matrix_dot(2, 2)
        integer, intent(out) :: status

        complex(dp) :: beta, beta_dot, imaginary_unit
        real(dp) :: beta_real, beta_real_dot, cutoff_scale
        real(dp) :: xi(2), xi_dot(2), q, q_dot, radicand
        integer :: a, b, mode_x, mode_y

        mode_x = i - 1
        if (mode_x > nx/2) mode_x = mode_x - nx
        mode_y = j - 1
        if (mode_y > ny/2) mode_y = mode_y - ny
        xi = 2.0_dp*acos(-1.0_dp)*[ &
            real(mode_x, dp)/length_x, real(mode_y, dp)/length_y]
        xi_dot = -xi*[length_x_dot/length_x, length_y_dot/length_y]
        q = sqrt(sum(xi**2))
        if (q > 0.0_dp) then
            q_dot = sum(xi*xi_dot)/q
        else
            q_dot = 0.0_dp
        end if
        radicand = wave_number**2 - q**2
        cutoff_scale = max(1.0_dp, wave_number**2, q**2)
        status = 1
        if (abs(radicand) <= 64.0_dp*epsilon(1.0_dp)*cutoff_scale) return
        if (radicand > 0.0_dp) then
            call generated_planar_helmholtz_propagating_jvp( &
                wave_number, q, wave_number_dot, q_dot, beta_real_dot)
            beta_real = sqrt(radicand)
            beta = cmplx(beta_real, 0.0_dp, dp)
            beta_dot = cmplx(beta_real_dot, 0.0_dp, dp)
        else
            call generated_planar_helmholtz_evanescent_jvp( &
                wave_number, q, wave_number_dot, q_dot, beta_real_dot)
            beta_real = sqrt(-radicand)
            beta = cmplx(0.0_dp, beta_real, dp)
            beta_dot = cmplx(0.0_dp, beta_real_dot, dp)
        end if
        imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)
        do b = 1, 2
            do a = 1, 2
                matrix(a, b) = -imaginary_unit*xi(a)*xi(b)/beta
                matrix_dot(a, b) = -imaginary_unit*( &
                    (xi_dot(a)*xi(b) + xi(a)*xi_dot(b))/beta - &
                    xi(a)*xi(b)*beta_dot/beta**2)
                if (a == b) then
                    matrix(a, b) = matrix(a, b) - imaginary_unit*beta
                    matrix_dot(a, b) = &
                        matrix_dot(a, b) - imaginary_unit*beta_dot
                end if
            end do
        end do
        status = 0
    end subroutine modal_matrix_product

    subroutine transform_forward(values)
        complex(dp), intent(inout) :: values(:, :, :)

        complex(dp), allocatable :: line(:)
        integer :: component, i, j

        do component = 1, 2
            allocate(line(size(values, 2)))
            do j = 1, size(values, 3)
                line = values(component, :, j)
                call fft_c2c(line, -1)
                values(component, :, j) = line
            end do
            deallocate(line)
            allocate(line(size(values, 3)))
            do i = 1, size(values, 2)
                line = values(component, i, :)
                call fft_c2c(line, -1)
                values(component, i, :) = line
            end do
            deallocate(line)
        end do
    end subroutine transform_forward

    subroutine transform_inverse(values)
        complex(dp), intent(inout) :: values(:, :, :)

        complex(dp), allocatable :: line(:)
        integer :: component, i, j

        do component = 1, 2
            allocate(line(size(values, 3)))
            do i = 1, size(values, 2)
                line = values(component, i, :)
                call fft_c2c(line, 1)
                values(component, i, :) = line
            end do
            deallocate(line)
            allocate(line(size(values, 2)))
            do j = 1, size(values, 3)
                line = values(component, :, j)
                call fft_c2c(line, 1)
                values(component, :, j) = line
            end do
            deallocate(line)
        end do
    end subroutine transform_inverse

    subroutine validate_product_inputs( &
            trace, seed, result, wave_number, length_x, length_y, nx, ny, status)
        complex(dp), intent(in) :: trace(:, :, :), seed(:, :, :)
        complex(dp), intent(out) :: result(:, :, :)
        real(dp), intent(in) :: wave_number, length_x, length_y
        integer, intent(out) :: nx, ny, status

        status = 1
        nx = size(trace, 2)
        ny = size(trace, 3)
        if (size(trace, 1) /= 2 .or. nx < 2 .or. ny < 2) return
        if (any(shape(seed) /= shape(trace))) return
        if (any(shape(result) /= shape(trace))) return
        if (wave_number <= 0.0_dp) return
        if (length_x <= 0.0_dp .or. length_y <= 0.0_dp) return
        status = 0
    end subroutine validate_product_inputs

end module fortfem_planar_maxwell_dtn
