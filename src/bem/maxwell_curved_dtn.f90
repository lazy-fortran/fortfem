module fortfem_maxwell_curved_dtn
    !! Galerkin Maxwell trace-to-flux maps on curved surfaces.
    !!
    !! The matrices have an explicit mixed-trace interpretation.  `electric_form`
    !! maps a surface-current coefficient vector to weak tangential-electric
    !! data, `flux_form` maps the same current to weak tangential magnetic-flux
    !! data, and `trace_mass` maps primal tangential trace coefficients to the
    !! electric test dual.  The resulting discrete DtN map is
    !!
    !!   D = flux_form * electric_form^{-1} * trace_mass.
    !!
    !! This is the finite-dimensional Calderon construction used by the torus
    !! wrapper.  Its output may use a Buffa--Christiansen/RBC dual test space;
    !! it is therefore a weak trace map, not a pointwise interpolation.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_torus_curved_rwg, only: &
        assemble_maxwell_torus_curved_efie_rwg_3d, &
        assemble_maxwell_torus_curved_mfie_rwg_rbc_3d, &
        assemble_maxwell_torus_curved_rwg_mass_matrix
    use fortnum_linalg, only: dense_solve
    implicit none

    private

    public :: apply_maxwell_trace_to_flux
    public :: apply_maxwell_trace_to_flux_jvp
    public :: apply_maxwell_trace_to_flux_map
    public :: apply_maxwell_trace_to_flux_vjp
    public :: assemble_maxwell_trace_to_flux_map
    public :: assemble_maxwell_torus_curved_dtn_rwg_3d

contains

    subroutine assemble_maxwell_trace_to_flux_map( &
            electric_form, flux_form, trace_mass, map, status)
        complex(dp), intent(in) :: electric_form(:, :), flux_form(:, :)
        complex(dp), intent(in) :: trace_mass(:, :)
        complex(dp), allocatable, intent(out) :: map(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: current_map(:, :)

        status = 1
        if (allocated(map)) deallocate(map)
        if (.not. valid_square_forms( &
            electric_form, flux_form, trace_mass)) return
        call solve_left_matrix(electric_form, trace_mass, current_map, status)
        if (status /= 0) return
        allocate(map(size(flux_form, 1), size(trace_mass, 2)))
        map = matmul(flux_form, current_map)
        status = 0
    end subroutine assemble_maxwell_trace_to_flux_map

    subroutine apply_maxwell_trace_to_flux( &
            electric_form, flux_form, trace_mass, trace, flux, status)
        complex(dp), intent(in) :: electric_form(:, :), flux_form(:, :)
        complex(dp), intent(in) :: trace_mass(:, :), trace(:)
        complex(dp), intent(out) :: flux(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: current(:), electric_trace(:)

        flux = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_square_forms( &
            electric_form, flux_form, trace_mass)) return
        if (size(trace) /= size(trace_mass, 2)) return
        if (size(flux) /= size(flux_form, 1)) return
        allocate(electric_trace(size(trace_mass, 1)))
        electric_trace = matmul(trace_mass, trace)
        call solve_left_vector(electric_form, electric_trace, current, status)
        if (status /= 0) return
        flux = matmul(flux_form, current)
        status = 0
    end subroutine apply_maxwell_trace_to_flux

    subroutine apply_maxwell_trace_to_flux_map(map, trace, flux, status)
        complex(dp), intent(in) :: map(:, :), trace(:)
        complex(dp), intent(out) :: flux(:)
        integer, intent(out) :: status

        flux = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(map, 2) /= size(trace)) return
        if (size(map, 1) /= size(flux)) return
        flux = matmul(map, trace)
        status = 0
    end subroutine apply_maxwell_trace_to_flux_map

    subroutine apply_maxwell_trace_to_flux_jvp( &
            electric_form, flux_form, trace_mass, trace, electric_form_dot, &
            flux_form_dot, trace_mass_dot, trace_dot, flux_dot, status)
        complex(dp), intent(in) :: electric_form(:, :), flux_form(:, :)
        complex(dp), intent(in) :: trace_mass(:, :), trace(:)
        complex(dp), intent(in) :: electric_form_dot(:, :), flux_form_dot(:, :)
        complex(dp), intent(in) :: trace_mass_dot(:, :), trace_dot(:)
        complex(dp), intent(out) :: flux_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: current(:), current_dot(:)
        complex(dp), allocatable :: electric_trace(:), electric_trace_dot(:)

        flux_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_square_forms( &
            electric_form, flux_form, trace_mass)) return
        if (.not. same_shape(electric_form, electric_form_dot)) return
        if (.not. same_shape(flux_form, flux_form_dot)) return
        if (.not. same_shape(trace_mass, trace_mass_dot)) return
        if (size(trace) /= size(trace_dot)) return
        if (size(trace) /= size(trace_mass, 2)) return
        if (size(flux_dot) /= size(flux_form, 1)) return
        allocate( &
            electric_trace(size(trace_mass, 1)), &
            electric_trace_dot(size(trace_mass, 1)))
        electric_trace = matmul(trace_mass, trace)
        electric_trace_dot = matmul(trace_mass_dot, trace) + &
            matmul(trace_mass, trace_dot)
        call solve_left_vector(electric_form, electric_trace, current, status)
        if (status /= 0) return
        electric_trace_dot = electric_trace_dot - &
            matmul(electric_form_dot, current)
        call solve_left_vector( &
            electric_form, electric_trace_dot, current_dot, status)
        if (status /= 0) return
        flux_dot = matmul(flux_form_dot, current) + &
            matmul(flux_form, current_dot)
        status = 0
    end subroutine apply_maxwell_trace_to_flux_jvp

    subroutine apply_maxwell_trace_to_flux_vjp( &
            electric_form, flux_form, trace_mass, trace, flux_bar, trace_bar, &
            electric_form_bar, flux_form_bar, trace_mass_bar, status)
        complex(dp), intent(in) :: electric_form(:, :), flux_form(:, :)
        complex(dp), intent(in) :: trace_mass(:, :), trace(:), flux_bar(:)
        complex(dp), intent(out) :: trace_bar(:)
        complex(dp), intent(out) :: electric_form_bar(:, :), flux_form_bar(:, :)
        complex(dp), intent(out) :: trace_mass_bar(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: current(:), electric_trace(:)
        complex(dp), allocatable :: current_bar(:), lambda(:)

        trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        electric_form_bar = cmplx(0.0_dp, 0.0_dp, dp)
        flux_form_bar = cmplx(0.0_dp, 0.0_dp, dp)
        trace_mass_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_square_forms( &
            electric_form, flux_form, trace_mass)) return
        if (size(trace) /= size(trace_mass, 2)) return
        if (size(flux_bar) /= size(flux_form, 1)) return
        if (size(trace_bar) /= size(trace_mass, 2)) return
        if (.not. same_shape(electric_form, electric_form_bar)) return
        if (.not. same_shape(flux_form, flux_form_bar)) return
        if (.not. same_shape(trace_mass, trace_mass_bar)) return
        allocate(electric_trace(size(trace_mass, 1)))
        electric_trace = matmul(trace_mass, trace)
        call solve_left_vector(electric_form, electric_trace, current, status)
        if (status /= 0) return
        allocate(current_bar(size(current)), lambda(size(current)))
        current_bar = matmul(conjg(transpose(flux_form)), flux_bar)
        call solve_left_vector( &
            conjg(transpose(electric_form)), current_bar, lambda, status)
        if (status /= 0) return
        trace_bar = matmul(conjg(transpose(trace_mass)), lambda)
        call rank_one_product(flux_form_bar, flux_bar, current)
        call rank_one_product(trace_mass_bar, lambda, trace)
        call rank_one_product(electric_form_bar, lambda, current)
        electric_form_bar = -electric_form_bar
        status = 0
    end subroutine apply_maxwell_trace_to_flux_vjp

    subroutine assemble_maxwell_torus_curved_dtn_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            mfie_offset, map, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: impedance, tolerance, mfie_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: map(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: electric_form(:, :), flux_form(:, :)
        real(dp), allocatable :: real_mass(:, :)
        complex(dp), allocatable :: trace_mass(:, :)

        status = 1
        if (allocated(map)) deallocate(map)
        call assemble_maxwell_torus_curved_efie_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            electric_form, status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_mfie_rwg_rbc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, mfie_offset, flux_form, status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_rwg_mass_matrix( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, real_mass, status)
        if (status /= 0) return
        allocate(trace_mass(size(real_mass, 1), size(real_mass, 2)))
        trace_mass = cmplx(real_mass, 0.0_dp, dp)
        call assemble_maxwell_trace_to_flux_map( &
            electric_form, flux_form, trace_mass, map, status)
    end subroutine assemble_maxwell_torus_curved_dtn_rwg_3d

    logical pure function valid_square_forms( &
            electric_form, flux_form, trace_mass) result(valid)
        complex(dp), intent(in) :: electric_form(:, :), flux_form(:, :)
        complex(dp), intent(in) :: trace_mass(:, :)

        valid = .false.
        if (size(electric_form, 1) /= size(electric_form, 2)) return
        if (size(flux_form, 2) /= size(electric_form, 2)) return
        if (size(trace_mass, 1) /= size(electric_form, 1)) return
        if (size(trace_mass, 2) /= size(electric_form, 2)) return
        if (size(electric_form, 1) == 0) return
        valid = .true.
    end function valid_square_forms

    logical pure function same_shape(first, second) result(same)
        complex(dp), intent(in) :: first(:, :), second(:, :)

        same = all(shape(first) == shape(second))
    end function same_shape

    subroutine solve_left_matrix(matrix, rhs, solution, status)
        complex(dp), intent(in) :: matrix(:, :), rhs(:, :)
        complex(dp), allocatable, intent(out) :: solution(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix_work(:, :), rhs_work(:, :)
        integer :: info

        status = 1
        if (size(matrix, 1) /= size(matrix, 2)) return
        if (size(rhs, 1) /= size(matrix, 1)) return
        allocate( &
            matrix_work(size(matrix, 1), size(matrix, 2)), &
            rhs_work(size(rhs, 1), size(rhs, 2)), &
            solution(size(rhs, 1), size(rhs, 2)))
        matrix_work = matrix
        rhs_work = rhs
        call dense_solve(matrix_work, rhs_work, solution, info)
        if (info /= 0) return
        status = 0
    end subroutine solve_left_matrix

    subroutine solve_left_vector(matrix, rhs, solution, status)
        complex(dp), intent(in) :: matrix(:, :), rhs(:)
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix_work(:, :), rhs_work(:)
        integer :: info

        status = 1
        if (size(matrix, 1) /= size(matrix, 2)) return
        if (size(rhs) /= size(matrix, 1)) return
        allocate( &
            matrix_work(size(matrix, 1), size(matrix, 2)), &
            rhs_work(size(rhs)), solution(size(rhs)))
        matrix_work = matrix
        rhs_work = rhs
        call dense_solve(matrix_work, rhs_work, solution, info)
        if (info /= 0) return
        status = 0
    end subroutine solve_left_vector

    subroutine rank_one_product(target, left, right)
        complex(dp), intent(out) :: target(:, :)
        complex(dp), intent(in) :: left(:), right(:)
        integer :: i, j

        do j = 1, size(target, 2)
            do i = 1, size(target, 1)
                target(i, j) = left(i)*conjg(right(j))
            end do
        end do
    end subroutine rank_one_product

end module fortfem_maxwell_curved_dtn
