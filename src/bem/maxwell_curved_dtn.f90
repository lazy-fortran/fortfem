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
    public :: apply_maxwell_weak_trace_reconstruction
    public :: apply_maxwell_weak_trace_reconstruction_jvp
    public :: apply_maxwell_weak_trace_reconstruction_vjp
    public :: assemble_maxwell_trace_to_flux_map
    public :: assemble_maxwell_weak_trace_reconstruction
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

    subroutine assemble_maxwell_weak_trace_reconstruction( &
            weak_mass, point_basis, reconstruction, status)
        !! Build point values from weak dual coefficients.
        !!
        !! If `A` is the caller-owned weak test/primal mass and `B` evaluates
        !! the primal trace basis at points, the returned map is `B A^{-1}`.
        !! The point basis and mass remain caller-owned so this composition is
        !! valid for RWG/RBC, scalar traces, IGA patches, and cut spaces.
        complex(dp), intent(in) :: weak_mass(:, :), point_basis(:, :)
        complex(dp), allocatable, intent(out) :: reconstruction(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: transposed_solution(:, :)

        status = 1
        if (allocated(reconstruction)) deallocate(reconstruction)
        if (.not. valid_reconstruction_inputs(weak_mass, point_basis)) return
        call solve_left_matrix( &
            transpose(weak_mass), transpose(point_basis), &
            transposed_solution, status)
        if (status /= 0) return
        allocate(reconstruction(size(point_basis, 1), size(point_basis, 2)))
        reconstruction = transpose(transposed_solution)
        status = 0
    end subroutine assemble_maxwell_weak_trace_reconstruction

    subroutine apply_maxwell_weak_trace_reconstruction( &
            weak_mass, point_basis, weak_trace, point_values, status)
        complex(dp), intent(in) :: weak_mass(:, :), point_basis(:, :)
        complex(dp), intent(in) :: weak_trace(:)
        complex(dp), intent(out) :: point_values(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: coefficients(:)

        point_values = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_reconstruction_inputs(weak_mass, point_basis)) return
        if (size(weak_trace) /= size(weak_mass, 1)) return
        if (size(point_values) /= size(point_basis, 1)) return
        call solve_left_vector(weak_mass, weak_trace, coefficients, status)
        if (status /= 0) return
        point_values = matmul(point_basis, coefficients)
        status = 0
    end subroutine apply_maxwell_weak_trace_reconstruction

    subroutine apply_maxwell_weak_trace_reconstruction_jvp( &
            weak_mass, point_basis, weak_trace, weak_mass_dot, point_basis_dot, &
            weak_trace_dot, point_values_dot, status)
        complex(dp), intent(in) :: weak_mass(:, :), point_basis(:, :)
        complex(dp), intent(in) :: weak_trace(:)
        complex(dp), intent(in) :: weak_mass_dot(:, :), point_basis_dot(:, :)
        complex(dp), intent(in) :: weak_trace_dot(:)
        complex(dp), intent(out) :: point_values_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: coefficients(:), coefficients_dot(:)
        complex(dp), allocatable :: right_hand_side(:)

        point_values_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_reconstruction_inputs(weak_mass, point_basis)) return
        if (.not. same_shape(weak_mass, weak_mass_dot)) return
        if (.not. same_shape(point_basis, point_basis_dot)) return
        if (size(weak_trace) /= size(weak_mass, 1) .or. &
            size(weak_trace_dot) /= size(weak_trace)) return
        if (size(point_values_dot) /= size(point_basis, 1)) return
        call solve_left_vector(weak_mass, weak_trace, coefficients, status)
        if (status /= 0) return
        allocate(right_hand_side(size(weak_trace)), coefficients_dot(size(weak_trace)))
        right_hand_side = weak_trace_dot - matmul(weak_mass_dot, coefficients)
        call solve_left_vector( &
            weak_mass, right_hand_side, coefficients_dot, status)
        if (status /= 0) return
        point_values_dot = matmul(point_basis_dot, coefficients) + &
            matmul(point_basis, coefficients_dot)
        status = 0
    end subroutine apply_maxwell_weak_trace_reconstruction_jvp

    subroutine apply_maxwell_weak_trace_reconstruction_vjp( &
            weak_mass, point_basis, weak_trace, point_values_bar, &
            weak_trace_bar, weak_mass_bar, point_basis_bar, status)
        complex(dp), intent(in) :: weak_mass(:, :), point_basis(:, :)
        complex(dp), intent(in) :: weak_trace(:), point_values_bar(:)
        complex(dp), intent(out) :: weak_trace_bar(:)
        complex(dp), intent(out) :: weak_mass_bar(:, :), point_basis_bar(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: coefficients(:), coefficients_bar(:)
        complex(dp), allocatable :: lambda(:)

        weak_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        weak_mass_bar = cmplx(0.0_dp, 0.0_dp, dp)
        point_basis_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_reconstruction_inputs(weak_mass, point_basis)) return
        if (size(weak_trace) /= size(weak_mass, 1)) return
        if (size(point_values_bar) /= size(point_basis, 1)) return
        if (size(weak_trace_bar) /= size(weak_trace)) return
        if (.not. same_shape(weak_mass, weak_mass_bar)) return
        if (.not. same_shape(point_basis, point_basis_bar)) return
        call solve_left_vector(weak_mass, weak_trace, coefficients, status)
        if (status /= 0) return
        allocate(coefficients_bar(size(coefficients)), lambda(size(coefficients)))
        coefficients_bar = matmul( &
            conjg(transpose(point_basis)), point_values_bar)
        call solve_left_vector( &
            conjg(transpose(weak_mass)), coefficients_bar, lambda, status)
        if (status /= 0) return
        weak_trace_bar = lambda
        call rank_one_product(point_basis_bar, point_values_bar, coefficients)
        call rank_one_product(weak_mass_bar, lambda, coefficients)
        weak_mass_bar = -weak_mass_bar
        status = 0
    end subroutine apply_maxwell_weak_trace_reconstruction_vjp

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

    logical pure function valid_reconstruction_inputs( &
            weak_mass, point_basis) result(valid)
        complex(dp), intent(in) :: weak_mass(:, :), point_basis(:, :)

        valid = .false.
        if (size(weak_mass, 1) /= size(weak_mass, 2)) return
        if (size(point_basis, 2) /= size(weak_mass, 1)) return
        if (size(weak_mass, 1) == 0 .or. size(point_basis, 1) == 0) return
        valid = .true.
    end function valid_reconstruction_inputs

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
