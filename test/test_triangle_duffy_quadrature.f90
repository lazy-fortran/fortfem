program test_triangle_duffy_quadrature
    use check, only: check_condition, check_summary
    use fortfem_feec, only: triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: weights(:), xi(:), eta(:)
    real(dp) :: exact, numerical
    integer :: degree, status, x_degree, y_degree
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 12
        call triangle_duffy_quadrature(degree, xi, eta, weights, status)
        call record_condition(status == 0, &
            "Arbitrary-degree triangle quadrature constructs its rule")
        do x_degree = 0, degree
            do y_degree = 0, degree - x_degree
                numerical = sum( &
                    weights * xi**x_degree * eta**y_degree)
                exact = factorial(x_degree) * factorial(y_degree) / &
                    factorial(x_degree + y_degree + 2)
                call record_condition(abs(numerical - exact) < 3.0e-14_dp, &
                    "Duffy rule integrates a triangle monomial exactly")
            end do
        end do
    end do

    call triangle_duffy_quadrature(-1, xi, eta, weights, status)
    call record_condition(status /= 0, &
        "Triangle quadrature rejects a negative polynomial degree")

    call check_summary("Arbitrary-degree triangle Duffy quadrature")
    if (.not. all_passed) error stop 1

contains

    pure function factorial(n) result(value)
        integer, intent(in) :: n
        real(dp) :: value

        integer :: factor

        value = 1.0_dp
        do factor = 2, n
            value = value * real(factor, dp)
        end do
    end function factorial

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_duffy_quadrature
