module fortfem_circular_dtn_2d
    use fortfem_kinds, only: dp
    use fortnum_special_complex_bessel, only: hankel_h1_real
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: circular_helmholtz_dtn_eigenvalue

contains

    subroutine circular_helmholtz_dtn_eigenvalue( &
            mode, wavenumber, radius, eigenvalue, status)
        integer, intent(in) :: mode
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(out) :: eigenvalue
        integer, intent(out) :: status

        complex(dp) :: derivative, hankel, hankel_minus, hankel_plus
        type(fortnum_status_t) :: special_status
        real(dp) :: argument
        integer :: order

        eigenvalue = (0.0_dp, 0.0_dp)
        status = -1
        if (wavenumber <= 0.0_dp .or. radius <= 0.0_dp) return

        order = abs(mode)
        argument = wavenumber * radius
        call hankel_h1_real(order, argument, hankel, special_status)
        if (.not. status_ok(special_status)) then
            status = -2
            return
        end if
        if (order == 0) then
            call hankel_h1_real(1, argument, hankel_plus, special_status)
            derivative = -hankel_plus
        else
            call hankel_h1_real( &
                order - 1, argument, hankel_minus, special_status)
            call hankel_h1_real( &
                order + 1, argument, hankel_plus, special_status)
            derivative = 0.5_dp * (hankel_minus - hankel_plus)
        end if
        if (.not. status_ok(special_status)) then
            status = -2
            return
        end if

        eigenvalue = wavenumber * derivative / hankel
        status = 0
    end subroutine circular_helmholtz_dtn_eigenvalue

end module fortfem_circular_dtn_2d
