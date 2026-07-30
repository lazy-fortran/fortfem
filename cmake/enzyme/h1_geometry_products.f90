module fortfem_enzyme_h1_geometry_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private

    public :: h1_geometry_objective

contains

    pure function h1_geometry_objective(geometry, data) result(value) &
            bind(c, name="fortfem_enzyme_h1_geometry_objective")
        real(c_double), intent(in) :: geometry(4), data(10)
        real(c_double) :: column_physical(2), determinant
        real(c_double) :: contribution, row_physical(2), value

        determinant = geometry(1)*geometry(4) - geometry(3)*geometry(2)
        row_physical(1) = ( &
            geometry(4)*data(1) - geometry(2)*data(2))/determinant
        row_physical(2) = ( &
            -geometry(3)*data(1) + geometry(1)*data(2))/determinant
        column_physical(1) = ( &
            geometry(4)*data(3) - geometry(2)*data(4))/determinant
        column_physical(2) = ( &
            -geometry(3)*data(3) + geometry(1)*data(4))/determinant
        contribution = data(9)*determinant*( &
            data(7)*dot_product(row_physical, column_physical) + &
            data(8)*data(5)*data(6))
        value = data(10)*contribution
    end function h1_geometry_objective

end module fortfem_enzyme_h1_geometry_kernel

program fortfem_enzyme_h1_geometry_products
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use fortfem_enzyme_h1_geometry_kernel, only: h1_geometry_objective
    use fortfem_generated_bspline_h1_geometry_jvp, only: &
        generated_bspline_h1_geometry_jvp
    use fortfem_generated_bspline_h1_geometry_vjp, only: &
        generated_bspline_h1_geometry_vjp
    implicit none

    interface
        function enzyme_fwddiff( &
                function_pointer, geometry, geometry_dot, data, data_dot) &
                result(value_dot) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: function_pointer
            real(c_double), intent(in) :: geometry(*), geometry_dot(*)
            real(c_double), intent(in) :: data(*), data_dot(*)
            real(c_double) :: value_dot
        end function enzyme_fwddiff

        function enzyme_autodiff( &
                function_pointer, geometry, geometry_bar, data, data_bar) &
                result(value) bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: function_pointer
            real(c_double), intent(in) :: geometry(*), data(*)
            real(c_double), intent(inout) :: geometry_bar(*), data_bar(*)
            real(c_double) :: value
        end function enzyme_autodiff
    end interface

    integer, parameter :: benchmark_iterations = 1000000
    real(c_double), parameter :: tolerance = 2.0e-11_c_double
    real(c_double), parameter :: difference_step = 1.0e-6_c_double
    real(c_double) :: analytical_bar(4), analytical_dot
    real(c_double) :: analytical_jvp_seconds, analytical_vjp_seconds
    real(c_double) :: automatic_bar(4), automatic_dot
    real(c_double) :: automatic_jvp_seconds, automatic_vjp_seconds
    real(c_double) :: data(10), data_bar(10), data_dot(10)
    real(c_double) :: finite_difference
    real(c_double) :: geometry(4), geometry_dot(4)
    real(c_double) :: geometry_reference(4)
    real(c_double) :: geometry_minus(4), geometry_plus(4)
    real(c_double) :: data_reference(10)
    real(c_double), volatile :: autodiff_return, sink
    real(c_double) :: start_time
    integer :: i

    geometry = [1.2_c_double, 0.1_c_double, -0.2_c_double, 0.9_c_double]
    geometry_dot = [ &
        0.03_c_double, -0.02_c_double, 0.01_c_double, 0.04_c_double]
    data = [ &
        0.4_c_double, -0.7_c_double, 0.2_c_double, 0.8_c_double, &
        0.6_c_double, 0.3_c_double, 1.1_c_double, 0.5_c_double, &
        0.07_c_double, -0.9_c_double]
    geometry_reference = geometry
    data_reference = data
    data_dot = 0.0_c_double
    call analytical_jvp(analytical_dot)
    geometry_minus = geometry - difference_step*geometry_dot
    geometry_plus = geometry + difference_step*geometry_dot
    finite_difference = ( &
        h1_geometry_objective(geometry_plus, data) - &
        h1_geometry_objective(geometry_minus, data))/(2.0_c_double*difference_step)
    if (abs(analytical_dot - finite_difference) > &
        2.0e-9_c_double*max(1.0_c_double, abs(finite_difference))) then
        error stop "FortSym H1 geometry JVP does not match central difference"
    end if
    automatic_dot = enzyme_fwddiff( &
        c_funloc(h1_geometry_objective), geometry, geometry_dot, data, data_dot)
    if (abs(automatic_dot - analytical_dot) > &
        tolerance*max(1.0_c_double, abs(analytical_dot))) then
        error stop "Enzyme H1 geometry JVP does not match FortSym"
    end if

    call analytical_vjp(analytical_bar)
    if (abs(dot_product(analytical_bar, geometry_dot) - analytical_dot) > &
        tolerance*max(1.0_c_double, abs(analytical_dot))) then
        error stop "FortSym H1 geometry products violate the dot identity"
    end if
    automatic_bar = 0.0_c_double
    data_bar = 0.0_c_double
    autodiff_return = enzyme_autodiff( &
        c_funloc(h1_geometry_objective), geometry, automatic_bar, data, data_bar)
    if (any(geometry /= geometry_reference) .or. any(data /= data_reference)) then
        error stop "Enzyme H1 geometry VJP mutated a primal argument"
    end if
    if (maxval(abs(automatic_bar - analytical_bar)) > &
        tolerance*max(1.0_c_double, maxval(abs(analytical_bar)))) then
        error stop "Enzyme H1 geometry VJP does not match FortSym"
    end if

    sink = 0.0_c_double
    call cpu_time(start_time)
    do i = 1, benchmark_iterations
        geometry(1) = geometry_reference(1) + &
            1.0e-8_c_double*real(mod(i, 251), c_double)
        call analytical_jvp(analytical_dot)
        sink = analytical_dot
    end do
    call cpu_time(analytical_jvp_seconds)
    analytical_jvp_seconds = analytical_jvp_seconds - start_time
    call cpu_time(start_time)
    do i = 1, benchmark_iterations
        geometry(1) = geometry_reference(1) + &
            1.0e-8_c_double*real(mod(i, 251), c_double)
        automatic_dot = enzyme_fwddiff( &
            c_funloc(h1_geometry_objective), geometry, geometry_dot, &
            data, data_dot)
        sink = automatic_dot
    end do
    call cpu_time(automatic_jvp_seconds)
    automatic_jvp_seconds = automatic_jvp_seconds - start_time

    call cpu_time(start_time)
    do i = 1, benchmark_iterations
        geometry(1) = geometry_reference(1) + &
            1.0e-8_c_double*real(mod(i, 251), c_double)
        call analytical_vjp(analytical_bar)
        sink = sum(analytical_bar)
    end do
    call cpu_time(analytical_vjp_seconds)
    analytical_vjp_seconds = analytical_vjp_seconds - start_time
    call cpu_time(start_time)
    do i = 1, benchmark_iterations
        geometry(1) = geometry_reference(1) + &
            1.0e-8_c_double*real(mod(i, 251), c_double)
        automatic_bar = 0.0_c_double
        data_bar = 0.0_c_double
        autodiff_return = enzyme_autodiff( &
            c_funloc(h1_geometry_objective), geometry, automatic_bar, &
            data, data_bar)
        sink = sum(automatic_bar)
    end do
    call cpu_time(automatic_vjp_seconds)
    automatic_vjp_seconds = automatic_vjp_seconds - start_time

    write (*, "(a)") "PASS"
    write (*, "(a,es12.4)") "analytical JVP seconds: ", analytical_jvp_seconds
    write (*, "(a,es12.4)") "Enzyme JVP seconds: ", automatic_jvp_seconds
    write (*, "(a,es12.4)") "analytical VJP seconds: ", analytical_vjp_seconds
    write (*, "(a,es12.4)") "Enzyme VJP seconds: ", automatic_vjp_seconds
    write (*, "(a,f8.3)") "Enzyme/analytical JVP ratio: ", &
        automatic_jvp_seconds/analytical_jvp_seconds
    write (*, "(a,f8.3)") "Enzyme/analytical VJP ratio: ", &
        automatic_vjp_seconds/analytical_vjp_seconds
    if (sink /= sink) &
        error stop "invalid differentiation benchmark sink"

contains

    subroutine analytical_jvp(value_dot)
        real(c_double), intent(out) :: value_dot

        call generated_bspline_h1_geometry_jvp( &
            geometry(1), geometry(2), geometry(3), geometry(4), &
            data(1), data(2), data(3), data(4), data(5), data(6), &
            data(7), data(8), data(9), geometry_dot(1), geometry_dot(2), &
            geometry_dot(3), geometry_dot(4), value_dot)
        value_dot = data(10)*value_dot
    end subroutine analytical_jvp

    subroutine analytical_vjp(value_bar)
        real(c_double), intent(out) :: value_bar(4)

        call generated_bspline_h1_geometry_vjp( &
            geometry(1), geometry(2), geometry(3), geometry(4), &
            data(1), data(2), data(3), data(4), data(5), data(6), &
            data(7), data(8), data(9), data(10), value_bar(1), &
            value_bar(2), value_bar(3), value_bar(4))
    end subroutine analytical_vjp

end program fortfem_enzyme_h1_geometry_products
