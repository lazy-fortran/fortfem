module triangle_level_context_cases
    use fortfem_kinds, only: dp
    use fortfem_triangle_level_quadrature, only: triangle_level_function_t
    implicit none
    type, extends(triangle_level_function_t) :: affine_integrand_t
        real(dp) :: scale = 1
    contains
        procedure :: evaluate => affine_evaluate
    end type
contains
    subroutine affine_evaluate(self, point, values, status)
        class(affine_integrand_t), intent(in) :: self
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status
        values = self%scale*(1 + point(1) + point(2))
        status = 0
    end subroutine
    subroutine legacy_affine(point, values, status)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status
        values = 1 + point(1) + point(2)
        status = 0
    end subroutine
end module
program test_triangle_level_context
    use triangle_level_context_cases
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_triangle_level_quadrature, only: triangle_level_workspace_t, &
        initialize_triangle_level_workspace, integrate_triangle_levels
    use check, only: check_condition, check_summary
    implicit none
    type(affine_integrand_t) :: first, second
    type(triangle_level_workspace_t) :: work
    real(dp) :: coefficients(6), levels(0), value(1), error(1), old(1), old_error(1)
    integer :: status
    character(180) :: message
    coefficients = 0
    call initialize_triangle_level_workspace(work, 1, 1, status)
    call check_condition(status == 0, 'context workspace')
    call integrate_triangle_levels(coefficients, levels, first, [1e-12_dp], &
        1e-11_dp, work, value, error, status, message)
    call check_condition(status == 0, 'class integrand converges')
    call check_condition(abs(value(1) - 5._dp/6) < 2e-12_dp, &
        'affine integral matches independent simplex moments')
    call integrate_triangle_levels(coefficients, levels, legacy_affine, [1e-12_dp], &
        1e-11_dp, work, old, old_error, status, message)
    call check_condition(status == 0, 'procedure adapter converges')
    call check_condition(all(transfer(value, [0_int64], 1) == &
        transfer(old, [0_int64], 1)), 'class and procedure values are bitwise equal')
    call check_condition(all(transfer(error, [0_int64], 1) == &
        transfer(old_error, [0_int64], 1)), &
            'class and procedure errors are bitwise equal')
    second%scale = 2
    call integrate_triangle_levels(coefficients, levels, second, [1e-12_dp], &
        1e-11_dp, work, value, error, status, message)
    call check_condition(status == 0 .and. abs(value(1) - 5._dp/3) < 2e-12_dp, &
        'independent contexts preserve their parameters')
    call check_summary('triangle level context')
end program
