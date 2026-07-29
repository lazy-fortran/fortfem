program test_tetra_duffy_quadrature
    use check, only: check_condition, check_summary
    use fortfem_api, only: tetra_duffy_quadrature
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: weights(:), x(:), y(:), z(:)
    real(dp) :: exact, numerical
    integer :: degree, status, x_degree, y_degree, z_degree
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 7
        call tetra_duffy_quadrature(degree, x, y, z, weights, status)
        call record_condition(status == 0, &
            "Tetrahedron Duffy rule accepts a nonnegative degree")
        call record_condition(abs(sum(weights) - 1.0_dp / 6.0_dp) < &
            2.0e-14_dp, "Tetrahedron Duffy weights have exact volume")
        do z_degree = 0, degree
            do y_degree = 0, degree - z_degree
                do x_degree = 0, degree - y_degree - z_degree
                    numerical = sum(weights * x**x_degree * &
                        y**y_degree * z**z_degree)
                    exact = monomial_integral( &
                        x_degree, y_degree, z_degree)
                    call record_condition(abs(numerical - exact) < &
                        4.0e-14_dp, &
                        "Duffy rule integrates a tetrahedron monomial exactly")
                end do
            end do
        end do
        deallocate(x, y, z, weights)
    end do

    call tetra_duffy_quadrature(-1, x, y, z, weights, status)
    call record_condition(status /= 0, &
        "Tetrahedron Duffy rule rejects a negative degree")

    call check_summary("Tetrahedron Duffy quadrature")
    if (.not. all_passed) error stop 1

contains

    pure function monomial_integral( &
            x_degree, y_degree, z_degree) result(value)
        integer, intent(in) :: x_degree, y_degree, z_degree
        real(dp) :: value

        value = factorial(x_degree) * factorial(y_degree) * &
            factorial(z_degree) / &
            factorial(x_degree + y_degree + z_degree + 3)
    end function monomial_integral

    pure function factorial(argument) result(value)
        integer, intent(in) :: argument
        real(dp) :: value

        integer :: factor

        value = 1.0_dp
        do factor = 2, argument
            value = value * real(factor, dp)
        end do
    end function factorial

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_duffy_quadrature
