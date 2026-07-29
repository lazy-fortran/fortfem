program test_triangle_quadrature_degree7
    use fortfem_kinds, only: dp
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    use check, only: check_condition, check_summary
    implicit none

    type(gauss_quadrature_triangle_t) :: quadrature
    real(dp) :: numerical, exact
    integer :: x_degree, y_degree
    logical :: all_passed

    all_passed = .true.
    quadrature = get_gauss_quadrature_triangle(7)

    do x_degree = 0, 7
        do y_degree = 0, 7 - x_degree
            numerical = sum(quadrature%weights * &
                quadrature%xi**x_degree * quadrature%eta**y_degree)
            exact = factorial(x_degree) * factorial(y_degree) / &
                factorial(x_degree + y_degree + 2)
            call record_condition(abs(numerical - exact) < 1.0e-12_dp, &
                "Degree-seven triangle rule integrates each monomial exactly")
        end do
    end do

    call check_summary("Degree-seven triangle quadrature")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function factorial(n) result(value)
        integer, intent(in) :: n
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

end program test_triangle_quadrature_degree7
