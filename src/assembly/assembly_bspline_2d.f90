module fortfem_assembly_bspline_2d
    !! Sparse scalar isogeometric assembly on one rational tensor patch.
    use fortfem_bspline_feec, only: &
        build_bspline_derivative_matrix, evaluate_bspline_basis, &
        evaluate_nurbs_surface_geometry, &
        evaluate_nurbs_surface_geometry_jvp, &
        evaluate_nurbs_surface_geometry_vjp
    use fortfem_generated_bspline_h1_geometry_jvp, only: &
        generated_bspline_h1_geometry_jvp
    use fortfem_generated_bspline_h1_geometry_vjp, only: &
        generated_bspline_h1_geometry_vjp
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: &
        csc_from_triplet, csc_is_valid, csc_matmul, csc_matvec, csc_t, &
        csc_transpose, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, sparse_factor, sparse_free, sparse_solve, &
        sparse_solver_t, status_set
    implicit none
    private

    public :: assemble_bspline_h1_operator_csc
    public :: assemble_bspline_h1_operator_csc_jvp
    public :: assemble_bspline_h1_operator_csc_vjp
    public :: assemble_bspline_h1_weighted_mass_csc
    public :: assemble_bspline_hcurl_operator_csc
    public :: assemble_bspline_hdiv_operator_csc
    public :: assemble_bspline_h1_hcurl_gradient_csc
    public :: assemble_bspline_hcurl_h1_adjoint_gradient_csc
    public :: assemble_bspline_l2_mass_csc
    public :: assemble_bspline_hcurl_l2_curl_csc
    public :: assemble_bspline_l2_hcurl_adjoint_curl_csc
    public :: assemble_bspline_grad_shafranov_csc
    public :: assemble_bspline_toroidal_fourier_laplacian_csc
    public :: assemble_bspline_poloidal_bracket_csc
    public :: apply_bspline_toroidal_poloidal_bracket
    public :: apply_bspline_toroidal_poloidal_bracket_jvp
    public :: apply_bspline_jorek_flux_rhs
    public :: apply_bspline_jorek_flux_jvp
    public :: apply_bspline_jorek_thermodynamic_rhs
    public :: apply_bspline_jorek_thermodynamic_jvp
    public :: apply_bspline_jorek_density_rhs
    public :: apply_bspline_jorek_density_jvp
    public :: project_bspline_toroidal_product
    public :: advance_bspline_jorek_poloidal_flux_midpoint
    public :: advance_bspline_jorek_poloidal_flux_midpoint_steps
    public :: apply_toroidal_fourier_derivative
    public :: build_bspline_feec_2d_operators_csc
    public :: scalar_weight_2d
    public :: tensor_weight_2d

    abstract interface
        pure subroutine scalar_weight_2d(point, value)
            import :: dp
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value
        end subroutine scalar_weight_2d

        pure subroutine tensor_weight_2d(point, value)
            import :: dp
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value(2, 2)
        end subroutine tensor_weight_2d
    end interface

contains

    subroutine advance_bspline_jorek_poloidal_flux_midpoint( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            electric_potential, quadrature_order, time_step, magnetic_flux, &
            status)
        !! Cayley/implicit-midpoint subflow for dt(psi)=R[psi,u], fixed u.
        !!
        !! The spatial bracket is skew, so this map exactly preserves the
        !! spline mass norm and is time reversible. It is suitable as one
        !! Poisson propagator in a symmetric Hamiltonian splitting.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: electric_potential(:), time_step
        real(dp), intent(inout) :: magnetic_flux(:)
        type(fortsparse_status_t), intent(out) :: status

        call advance_bspline_jorek_poloidal_flux_midpoint_steps( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            electric_potential, quadrature_order, time_step, 1, magnetic_flux, &
            status)
    end subroutine advance_bspline_jorek_poloidal_flux_midpoint

    subroutine advance_bspline_jorek_poloidal_flux_midpoint_steps( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            electric_potential, quadrature_order, time_step, step_count, &
            magnetic_flux, status, magnetic_flux_history)
        !! Repeated midpoint subflow retaining one sparse factorization.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: electric_potential(:), time_step
        integer, intent(in) :: step_count
        real(dp), intent(inout) :: magnetic_flux(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable, intent(out), optional :: &
            magnetic_flux_history(:, :)

        real(dp), allocatable :: new_flux(:), right_hand_side(:)
        type(csc_t) :: left_matrix, right_matrix
        type(sparse_solver_t) :: solver
        integer :: step

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "JOREK midpoint propagator requires a positive step count")
        if (step_count < 1) return
        if (present(magnetic_flux_history)) then
            allocate(magnetic_flux_history(size(magnetic_flux), step_count + 1))
            magnetic_flux_history(:, 1) = magnetic_flux
        end if

        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, left_matrix, status, &
            stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp, &
            advecting_coefficients=electric_potential, &
            advection_coefficient=0.5_dp*time_step, &
            advection_weight_function=radial_weight)
        if (status%code /= 0) return
        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, right_matrix, status, &
            stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp, &
            advecting_coefficients=electric_potential, &
            advection_coefficient=-0.5_dp*time_step, &
            advection_weight_function=radial_weight)
        if (status%code /= 0) return
        allocate(new_flux(size(magnetic_flux)))
        call sparse_factor(solver, left_matrix, status)
        if (status%code /= 0) then
            call sparse_free(solver)
            return
        end if
        do step = 1, step_count
            right_hand_side = csc_matvec(right_matrix, magnetic_flux)
            call sparse_solve(solver, right_hand_side, new_flux, status)
            if (status%code /= 0) then
                call sparse_free(solver)
                return
            end if
            magnetic_flux = new_flux
            if (present(magnetic_flux_history)) then
                magnetic_flux_history(:, step + 1) = magnetic_flux
            end if
        end do
        call sparse_free(solver)

    contains

        pure subroutine radial_weight(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = point(1)
        end subroutine radial_weight

    end subroutine advance_bspline_jorek_poloidal_flux_midpoint_steps

    subroutine assemble_bspline_h1_weighted_mass_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            field_coefficients, quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: field_coefficients(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call assemble_bspline_h1_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient=0.0_dp, &
            mass_coefficient=1.0_dp, &
            mass_field_coefficients=field_coefficients)
    end subroutine assemble_bspline_h1_weighted_mass_csc

    subroutine apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, advecting_coefficients, transported_coefficients, &
            quadrature_order, residual, status, bracket_weight_function)
        !! Energy-skew poloidal bracket with exact Fourier-mode convolution.
        !!
        !! This is the nonlinear mode coupling used by reduced-MHD models:
        !! mode n receives every retained pair p+q=n. The formulation follows
        !! the JOREK reduced-MHD weak structure in Franck et al. (2015),
        !! arXiv:1408.2099, while retaining the skew Galerkin spatial bracket.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: advecting_coefficients(:, :)
        complex(dp), intent(in) :: transported_coefficients(:, :)
        complex(dp), allocatable, intent(out) :: residual(:, :)
        type(fortsparse_status_t), intent(out) :: status
        procedure(scalar_weight_2d), optional :: bracket_weight_function

        real(dp), allocatable :: action(:), advecting_part(:)
        type(csc_t) :: bracket_imaginary, bracket_real
        integer :: mode_advecting, mode_output, mode_transported
        logical :: has_imaginary, has_real

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Toroidal isogeometric poloidal-bracket convolution failed")
        if (size(advecting_coefficients, 1) /= &
            size(transported_coefficients, 1)) return
        if (size(advecting_coefficients, 2) /= size(modes) .or. &
            size(transported_coefficients, 2) /= size(modes)) return
        if (has_duplicate_modes(modes)) return
        allocate( &
            residual(size(transported_coefficients, 1), size(modes)), &
            advecting_part(size(advecting_coefficients, 1)))
        residual = cmplx(0.0_dp, 0.0_dp, dp)
        do mode_advecting = 1, size(modes)
            has_real = any(real( &
                advecting_coefficients(:, mode_advecting), dp) /= 0.0_dp)
            has_imaginary = any(aimag( &
                advecting_coefficients(:, mode_advecting)) /= 0.0_dp)
            if (has_real) then
                advecting_part = &
                    real(advecting_coefficients(:, mode_advecting), dp)
                call assemble_bspline_poloidal_bracket_csc( &
                    knots_r, knots_z, degree_r, degree_z, control_points, &
                    weights, advecting_part, quadrature_order, bracket_real, &
                    status, bracket_weight_function)
                if (status%code /= 0) return
            end if
            if (has_imaginary) then
                advecting_part = &
                    aimag(advecting_coefficients(:, mode_advecting))
                call assemble_bspline_poloidal_bracket_csc( &
                    knots_r, knots_z, degree_r, degree_z, control_points, &
                    weights, advecting_part, quadrature_order, &
                    bracket_imaginary, status, bracket_weight_function)
                if (status%code /= 0) return
            end if
            do mode_transported = 1, size(modes)
                mode_output = find_mode( &
                    modes, modes(mode_advecting) + modes(mode_transported))
                if (mode_output == 0) cycle
                if (has_real) then
                    action = csc_matvec( &
                        bracket_real, real( &
                        transported_coefficients(:, mode_transported), dp))
                    residual(:, mode_output) = residual(:, mode_output) + action
                    action = csc_matvec( &
                        bracket_real, aimag( &
                        transported_coefficients(:, mode_transported)))
                    residual(:, mode_output) = residual(:, mode_output) + &
                        cmplx(0.0_dp, 1.0_dp, dp)*action
                end if
                if (has_imaginary) then
                    action = csc_matvec( &
                        bracket_imaginary, real( &
                        transported_coefficients(:, mode_transported), dp))
                    residual(:, mode_output) = residual(:, mode_output) + &
                        cmplx(0.0_dp, 1.0_dp, dp)*action
                    action = csc_matvec( &
                        bracket_imaginary, aimag( &
                        transported_coefficients(:, mode_transported)))
                    residual(:, mode_output) = residual(:, mode_output) - action
                end if
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_bspline_toroidal_poloidal_bracket

    subroutine apply_bspline_toroidal_poloidal_bracket_jvp( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, advecting_coefficients, transported_coefficients, &
            advecting_direction, transported_direction, quadrature_order, &
            jvp, status, bracket_weight_function)
        !! Exact directional derivative of the bilinear Fourier bracket.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: advecting_coefficients(:, :)
        complex(dp), intent(in) :: transported_coefficients(:, :)
        complex(dp), intent(in) :: advecting_direction(:, :)
        complex(dp), intent(in) :: transported_direction(:, :)
        complex(dp), allocatable, intent(out) :: jvp(:, :)
        type(fortsparse_status_t), intent(out) :: status
        procedure(scalar_weight_2d), optional :: bracket_weight_function

        complex(dp), allocatable :: advecting_action(:, :)
        complex(dp), allocatable :: transported_action(:, :)

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Toroidal isogeometric bracket JVP failed")
        if (any(shape(advecting_coefficients) /= &
            shape(transported_coefficients))) return
        if (any(shape(advecting_coefficients) /= &
            shape(advecting_direction))) return
        if (any(shape(advecting_coefficients) /= &
            shape(transported_direction))) return
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, advecting_direction, transported_coefficients, &
            quadrature_order, advecting_action, status, &
            bracket_weight_function)
        if (status%code /= 0) return
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, advecting_coefficients, transported_direction, &
            quadrature_order, transported_action, status, &
            bracket_weight_function)
        if (status%code /= 0) return
        jvp = advecting_action + transported_action
    end subroutine apply_bspline_toroidal_poloidal_bracket_jvp

    subroutine apply_bspline_jorek_thermodynamic_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, field, electric_potential, adiabatic_index, &
            quadrature_order, rhs, status)
        !! No-parallel-flow part of Franck et al. equations (26)--(27):
        !! dt(q) = R[q,u] + 2*gamma_q*q*dZ(u).
        !!
        !! Use gamma_q=1 for density and the physical adiabatic index for
        !! pressure. Both terms are assembled as weak Galerkin loads.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: field(:, :), electric_potential(:, :)
        real(dp), intent(in) :: adiabatic_index
        complex(dp), allocatable, intent(out) :: rhs(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: bracket(:, :), compression(:, :)

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "JOREK isogeometric thermodynamic residual failed")
        if (adiabatic_index <= 0.0_dp) return
        if (any(shape(field) /= shape(electric_potential))) return
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, field, electric_potential, quadrature_order, bracket, &
            status, radial_weight)
        if (status%code /= 0) return
        call apply_bspline_toroidal_product_derivative( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, field, electric_potential, [0.0_dp, 1.0_dp], &
            quadrature_order, compression, status)
        if (status%code /= 0) return
        rhs = bracket + 2.0_dp*adiabatic_index*compression

    contains

        pure subroutine radial_weight(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = point(1)
        end subroutine radial_weight

    end subroutine apply_bspline_jorek_thermodynamic_rhs

    subroutine apply_bspline_jorek_thermodynamic_jvp( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, field, electric_potential, field_direction, &
            potential_direction, adiabatic_index, quadrature_order, jvp, &
            status)
        !! Exact product-rule JVP of the no-parallel thermodynamic residual.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: field(:, :), electric_potential(:, :)
        complex(dp), intent(in) :: field_direction(:, :)
        complex(dp), intent(in) :: potential_direction(:, :)
        real(dp), intent(in) :: adiabatic_index
        complex(dp), allocatable, intent(out) :: jvp(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: field_action(:, :), potential_action(:, :)

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "JOREK isogeometric thermodynamic JVP failed")
        if (any(shape(field) /= shape(electric_potential))) return
        if (any(shape(field) /= shape(field_direction))) return
        if (any(shape(field) /= shape(potential_direction))) return
        call apply_bspline_jorek_thermodynamic_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, field_direction, electric_potential, adiabatic_index, &
            quadrature_order, field_action, status)
        if (status%code /= 0) return
        call apply_bspline_jorek_thermodynamic_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, field, potential_direction, adiabatic_index, &
            quadrature_order, potential_action, status)
        if (status%code /= 0) return
        jvp = field_action + potential_action
    end subroutine apply_bspline_jorek_thermodynamic_jvp

    subroutine project_bspline_toroidal_product( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            modes, first_coefficients, second_coefficients, quadrature_order, &
            product_coefficients, status)
        !! L2 Galerkin projection of a retained-mode spline product.
        !!
        !! The convolution is formed before truncation: output n receives all
        !! represented pairs p+q=n. One mass factorization serves every mode.
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: first_coefficients(:, :)
        complex(dp), intent(in) :: second_coefficients(:, :)
        complex(dp), allocatable, intent(out) :: product_coefficients(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: load(:, :)
        real(dp), allocatable :: action(:), field_part(:)
        real(dp), allocatable :: imaginary_solution(:), real_solution(:)
        type(csc_t) :: mass, weighted_imaginary, weighted_real
        type(sparse_solver_t) :: solver
        integer :: first_mode, output_mode, second_mode
        logical :: has_imaginary, has_real

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Toroidal spline product projection failed")
        if (any(shape(first_coefficients) /= shape(second_coefficients))) return
        if (size(first_coefficients, 2) /= size(modes)) return
        if (has_duplicate_modes(modes)) return
        allocate( &
            load(size(first_coefficients, 1), size(modes)), &
            field_part(size(first_coefficients, 1)))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        do first_mode = 1, size(modes)
            has_real = any(real( &
                first_coefficients(:, first_mode), dp) /= 0.0_dp)
            has_imaginary = any(aimag( &
                first_coefficients(:, first_mode)) /= 0.0_dp)
            if (has_real) then
                field_part = real(first_coefficients(:, first_mode), dp)
                call assemble_bspline_h1_weighted_mass_csc( &
                    knots_x, knots_y, degree_x, degree_y, control_points, &
                    weights, field_part, quadrature_order, weighted_real, &
                    status)
                if (status%code /= 0) return
            end if
            if (has_imaginary) then
                field_part = aimag(first_coefficients(:, first_mode))
                call assemble_bspline_h1_weighted_mass_csc( &
                    knots_x, knots_y, degree_x, degree_y, control_points, &
                    weights, field_part, quadrature_order, weighted_imaginary, &
                    status)
                if (status%code /= 0) return
            end if
            do second_mode = 1, size(modes)
                output_mode = find_mode( &
                    modes, modes(first_mode) + modes(second_mode))
                if (output_mode == 0) cycle
                if (has_real) then
                    action = csc_matvec(weighted_real, &
                        real(second_coefficients(:, second_mode), dp))
                    load(:, output_mode) = load(:, output_mode) + action
                    action = csc_matvec( &
                        weighted_real, aimag(second_coefficients(:, second_mode)))
                    load(:, output_mode) = load(:, output_mode) + &
                        cmplx(0.0_dp, 1.0_dp, dp)*action
                end if
                if (has_imaginary) then
                    action = csc_matvec(weighted_imaginary, &
                        real(second_coefficients(:, second_mode), dp))
                    load(:, output_mode) = load(:, output_mode) + &
                        cmplx(0.0_dp, 1.0_dp, dp)*action
                    action = csc_matvec(weighted_imaginary, &
                        aimag(second_coefficients(:, second_mode)))
                    load(:, output_mode) = load(:, output_mode) - action
                end if
            end do
        end do
        call assemble_bspline_h1_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, mass, status, stiffness_coefficient=0.0_dp, &
            mass_coefficient=1.0_dp)
        if (status%code /= 0) return
        call sparse_factor(solver, mass, status)
        if (status%code /= 0) then
            call sparse_free(solver)
            return
        end if
        allocate( &
            product_coefficients(size(load, 1), size(load, 2)), &
            real_solution(size(load, 1)), imaginary_solution(size(load, 1)))
        do output_mode = 1, size(modes)
            call sparse_solve( &
                solver, real(load(:, output_mode), dp), real_solution, status)
            if (status%code /= 0) then
                call sparse_free(solver)
                return
            end if
            call sparse_solve( &
                solver, aimag(load(:, output_mode)), imaginary_solution, status)
            if (status%code /= 0) then
                call sparse_free(solver)
                return
            end if
            product_coefficients(:, output_mode) = &
                cmplx(real_solution, imaginary_solution, dp)
        end do
        call sparse_free(solver)
    end subroutine project_bspline_toroidal_product

    subroutine apply_bspline_jorek_density_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density, electric_potential, magnetic_flux, &
            parallel_velocity, toroidal_field, quadrature_order, rhs, status)
        !! Complete weak density equation (26) of Franck et al.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: density(:, :), electric_potential(:, :)
        complex(dp), intent(in) :: magnetic_flux(:, :)
        complex(dp), intent(in) :: parallel_velocity(:, :)
        real(dp), intent(in) :: toroidal_field
        complex(dp), allocatable, intent(out) :: rhs(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: parallel_product(:, :)
        complex(dp), allocatable :: parallel_rhs(:, :)

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "JOREK isogeometric density residual failed")
        if (any(shape(density) /= shape(electric_potential))) return
        if (any(shape(density) /= shape(magnetic_flux))) return
        if (any(shape(density) /= shape(parallel_velocity))) return
        call apply_bspline_jorek_thermodynamic_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density, electric_potential, 1.0_dp, quadrature_order, &
            rhs, status)
        if (status%code /= 0) return
        call project_bspline_toroidal_product( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density, parallel_velocity, quadrature_order, &
            parallel_product, status)
        if (status%code /= 0) return
        call apply_bspline_jorek_density_parallel_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, parallel_product, magnetic_flux, toroidal_field, &
            quadrature_order, parallel_rhs, status)
        if (status%code /= 0) return
        rhs = rhs + parallel_rhs
    end subroutine apply_bspline_jorek_density_rhs

    subroutine apply_bspline_jorek_density_jvp( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density, electric_potential, magnetic_flux, &
            parallel_velocity, density_direction, potential_direction, &
            flux_direction, velocity_direction, toroidal_field, &
            quadrature_order, jvp, status)
        !! Exact product-rule JVP of the complete density equation.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: density(:, :), electric_potential(:, :)
        complex(dp), intent(in) :: magnetic_flux(:, :)
        complex(dp), intent(in) :: parallel_velocity(:, :)
        complex(dp), intent(in) :: density_direction(:, :)
        complex(dp), intent(in) :: potential_direction(:, :)
        complex(dp), intent(in) :: flux_direction(:, :)
        complex(dp), intent(in) :: velocity_direction(:, :)
        real(dp), intent(in) :: toroidal_field
        complex(dp), allocatable, intent(out) :: jvp(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: base_product(:, :), product_direction(:, :)
        complex(dp), allocatable :: product_part(:, :), work(:, :)

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "JOREK isogeometric density JVP failed")
        if (any(shape(density) /= shape(density_direction))) return
        if (any(shape(density) /= shape(potential_direction))) return
        if (any(shape(density) /= shape(flux_direction))) return
        if (any(shape(density) /= shape(velocity_direction))) return
        call apply_bspline_jorek_thermodynamic_jvp( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density, electric_potential, density_direction, &
            potential_direction, 1.0_dp, quadrature_order, jvp, status)
        if (status%code /= 0) return
        call project_bspline_toroidal_product( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density_direction, parallel_velocity, quadrature_order, &
            product_direction, status)
        if (status%code /= 0) return
        call project_bspline_toroidal_product( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density, velocity_direction, quadrature_order, work, status)
        if (status%code /= 0) return
        product_direction = product_direction + work
        call project_bspline_toroidal_product( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density, parallel_velocity, quadrature_order, base_product, &
            status)
        if (status%code /= 0) return
        call apply_bspline_jorek_density_parallel_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, product_direction, magnetic_flux, toroidal_field, &
            quadrature_order, product_part, status)
        if (status%code /= 0) return
        call apply_bspline_jorek_density_parallel_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, base_product, flux_direction, 0.0_dp, quadrature_order, &
            work, status)
        if (status%code /= 0) return
        jvp = jvp + product_part + work
    end subroutine apply_bspline_jorek_density_jvp

    subroutine apply_bspline_jorek_density_parallel_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density_velocity, magnetic_flux, toroidal_field, &
            quadrature_order, rhs, status)
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: density_velocity(:, :)
        complex(dp), intent(in) :: magnetic_flux(:, :)
        real(dp), intent(in) :: toroidal_field
        complex(dp), allocatable, intent(out) :: rhs(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: bracket(:, :), derivative(:, :)
        type(csc_t) :: inverse_radius_squared_mass
        integer :: local_status, mode

        call apply_toroidal_fourier_derivative( &
            modes, density_velocity, derivative, local_status)
        if (local_status /= 0) return
        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, inverse_radius_squared_mass, status, &
            stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp, &
            mass_weight_function=inverse_radius_squared)
        if (status%code /= 0) return
        allocate(rhs(size(density_velocity, 1), size(modes)))
        rhs = cmplx(0.0_dp, 0.0_dp, dp)
        do mode = 1, size(modes)
            call add_complex_matrix_action( &
                inverse_radius_squared_mass, derivative(:, mode), &
                -toroidal_field, rhs(:, mode))
        end do
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, density_velocity, magnetic_flux, quadrature_order, bracket, &
            status, inverse_radius)
        if (status%code /= 0) return
        rhs = rhs - bracket

    contains

        pure subroutine inverse_radius_squared(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = 1.0_dp/point(1)**2
        end subroutine inverse_radius_squared

        pure subroutine inverse_radius(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = 1.0_dp/point(1)
        end subroutine inverse_radius

    end subroutine apply_bspline_jorek_density_parallel_rhs

    subroutine apply_bspline_toroidal_product_derivative( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            modes, multiplying_coefficients, differentiated_coefficients, &
            direction, quadrature_order, residual, status)
        !! Weak Fourier convolution of a_h * direction.grad(b_h).
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: multiplying_coefficients(:, :)
        complex(dp), intent(in) :: differentiated_coefficients(:, :)
        real(dp), intent(in) :: direction(2)
        complex(dp), allocatable, intent(out) :: residual(:, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: action(:), multiplying_part(:)
        type(csc_t) :: matrix_imaginary, matrix_real
        integer :: differentiated_mode, multiplying_mode, output_mode
        logical :: has_imaginary, has_real

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Toroidal spline product-derivative convolution failed")
        if (any(shape(multiplying_coefficients) /= &
            shape(differentiated_coefficients))) return
        if (size(multiplying_coefficients, 2) /= size(modes)) return
        if (has_duplicate_modes(modes)) return
        allocate( &
            residual(size(multiplying_coefficients, 1), size(modes)), &
            multiplying_part(size(multiplying_coefficients, 1)))
        residual = cmplx(0.0_dp, 0.0_dp, dp)
        do multiplying_mode = 1, size(modes)
            has_real = any(real( &
                multiplying_coefficients(:, multiplying_mode), dp) /= 0.0_dp)
            has_imaginary = any(aimag( &
                multiplying_coefficients(:, multiplying_mode)) /= 0.0_dp)
            if (has_real) then
                multiplying_part = real( &
                    multiplying_coefficients(:, multiplying_mode), dp)
                call assemble_bspline_product_derivative_csc( &
                    knots_x, knots_y, degree_x, degree_y, control_points, &
                    weights, multiplying_part, direction, quadrature_order, &
                    matrix_real, status)
                if (status%code /= 0) return
            end if
            if (has_imaginary) then
                multiplying_part = aimag( &
                    multiplying_coefficients(:, multiplying_mode))
                call assemble_bspline_product_derivative_csc( &
                    knots_x, knots_y, degree_x, degree_y, control_points, &
                    weights, multiplying_part, direction, quadrature_order, &
                    matrix_imaginary, status)
                if (status%code /= 0) return
            end if
            do differentiated_mode = 1, size(modes)
                output_mode = find_mode( &
                    modes, modes(multiplying_mode) + &
                    modes(differentiated_mode))
                if (output_mode == 0) cycle
                if (has_real) then
                    action = csc_matvec(matrix_real, real( &
                        differentiated_coefficients(:, differentiated_mode), &
                        dp))
                    residual(:, output_mode) = &
                        residual(:, output_mode) + action
                    action = csc_matvec(matrix_real, aimag( &
                        differentiated_coefficients(:, differentiated_mode)))
                    residual(:, output_mode) = residual(:, output_mode) + &
                        cmplx(0.0_dp, 1.0_dp, dp)*action
                end if
                if (has_imaginary) then
                    action = csc_matvec(matrix_imaginary, real( &
                        differentiated_coefficients(:, differentiated_mode), &
                        dp))
                    residual(:, output_mode) = residual(:, output_mode) + &
                        cmplx(0.0_dp, 1.0_dp, dp)*action
                    action = csc_matvec(matrix_imaginary, aimag( &
                        differentiated_coefficients(:, differentiated_mode)))
                    residual(:, output_mode) = &
                        residual(:, output_mode) - action
                end if
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_bspline_toroidal_product_derivative

    subroutine assemble_bspline_product_derivative_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            multiplying_coefficients, direction, quadrature_order, matrix, &
            status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: multiplying_coefficients(:), direction(2)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call assemble_bspline_h1_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient=0.0_dp, &
            mass_coefficient=0.0_dp, &
            reaction_coefficients=multiplying_coefficients, &
            reaction_direction=direction, reaction_coefficient=1.0_dp)
    end subroutine assemble_bspline_product_derivative_csc

    subroutine apply_bspline_jorek_flux_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, magnetic_flux, electric_potential, toroidal_current, &
            resistivity, toroidal_field, quadrature_order, rhs, status)
        !! Weak right-hand side of Franck et al. (2015), equation (10):
        !! dt psi = R[psi,u] + eta*j - F0*dphi(u)
        !!          + eta/R^2*dphi_phi(psi).
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: magnetic_flux(:, :)
        complex(dp), intent(in) :: electric_potential(:, :)
        complex(dp), intent(in) :: toroidal_current(:, :)
        real(dp), intent(in) :: resistivity, toroidal_field
        complex(dp), allocatable, intent(out) :: rhs(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: derivative(:, :), flux_derivative(:, :)
        complex(dp), allocatable :: second_derivative(:, :)
        complex(dp), allocatable :: nonlinear(:, :)
        type(csc_t) :: inverse_radius_squared_mass, mass
        integer :: local_status, mode

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "JOREK isogeometric magnetic-flux residual failed")
        if (any(shape(magnetic_flux) /= shape(electric_potential)) .or. &
            any(shape(magnetic_flux) /= shape(toroidal_current))) return
        if (size(magnetic_flux, 2) /= size(modes)) return
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, magnetic_flux, electric_potential, quadrature_order, &
            nonlinear, status, radial_weight)
        if (status%code /= 0) return
        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, mass, status, stiffness_coefficient=0.0_dp, &
            mass_coefficient=1.0_dp)
        if (status%code /= 0) return
        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, inverse_radius_squared_mass, status, &
            stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp, &
            mass_weight_function=inverse_radius_squared)
        if (status%code /= 0) return
        call apply_toroidal_fourier_derivative( &
            modes, electric_potential, derivative, local_status)
        if (local_status /= 0) return
        call apply_toroidal_fourier_derivative( &
            modes, magnetic_flux, flux_derivative, local_status)
        if (local_status /= 0) return
        call apply_toroidal_fourier_derivative( &
            modes, flux_derivative, second_derivative, local_status)
        if (local_status /= 0) return
        rhs = nonlinear
        do mode = 1, size(modes)
            call add_complex_matrix_action( &
                mass, toroidal_current(:, mode), resistivity, rhs(:, mode))
            call add_complex_matrix_action( &
                mass, derivative(:, mode), -toroidal_field, rhs(:, mode))
            call add_complex_matrix_action( &
                inverse_radius_squared_mass, second_derivative(:, mode), &
                resistivity, rhs(:, mode))
        end do
        call status_set(status, FORTSPARSE_OK, "")

    contains

        pure subroutine radial_weight(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = point(1)
        end subroutine radial_weight

        pure subroutine inverse_radius_squared(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = 1.0_dp/point(1)**2
        end subroutine inverse_radius_squared

    end subroutine apply_bspline_jorek_flux_rhs

    subroutine apply_bspline_jorek_flux_jvp( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, magnetic_flux, electric_potential, magnetic_flux_direction, &
            electric_potential_direction, current_direction, resistivity, &
            toroidal_field, quadrature_order, jvp, status)
        !! Analytical Jacobian-vector product of apply_bspline_jorek_flux_rhs.
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: magnetic_flux(:, :)
        complex(dp), intent(in) :: electric_potential(:, :)
        complex(dp), intent(in) :: magnetic_flux_direction(:, :)
        complex(dp), intent(in) :: electric_potential_direction(:, :)
        complex(dp), intent(in) :: current_direction(:, :)
        real(dp), intent(in) :: resistivity, toroidal_field
        complex(dp), allocatable, intent(out) :: jvp(:, :)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: bracket_cross(:, :)
        complex(dp), allocatable :: bracket_direction(:, :)

        call apply_bspline_jorek_flux_rhs( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, magnetic_flux_direction, electric_potential_direction, &
            current_direction, resistivity, toroidal_field, quadrature_order, &
            jvp, status)
        if (status%code /= 0) return
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, magnetic_flux_direction, electric_potential, &
            quadrature_order, bracket_cross, status, radial_weight)
        if (status%code /= 0) return
        jvp = jvp + bracket_cross
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, magnetic_flux, electric_potential_direction, &
            quadrature_order, bracket_cross, status, radial_weight)
        if (status%code /= 0) return
        jvp = jvp + bracket_cross
        call apply_bspline_toroidal_poloidal_bracket( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            modes, magnetic_flux_direction, electric_potential_direction, &
            quadrature_order, bracket_direction, status, radial_weight)
        if (status%code /= 0) return
        jvp = jvp - bracket_direction
        call status_set(status, FORTSPARSE_OK, "")

    contains

        pure subroutine radial_weight(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = point(1)
        end subroutine radial_weight

    end subroutine apply_bspline_jorek_flux_jvp

    subroutine add_complex_matrix_action(matrix, vector, scale, result)
        type(csc_t), intent(in) :: matrix
        complex(dp), intent(in) :: vector(:)
        real(dp), intent(in) :: scale
        complex(dp), intent(inout) :: result(:)

        result = result + scale*csc_matvec(matrix, real(vector, dp))
        result = result + cmplx(0.0_dp, scale, dp)* &
            csc_matvec(matrix, aimag(vector))
    end subroutine add_complex_matrix_action

    pure subroutine apply_toroidal_fourier_derivative( &
            modes, coefficients, derivative, status)
        !! Exact derivative of exp(i*n*phi) Fourier coefficients.
        integer, intent(in) :: modes(:)
        complex(dp), intent(in) :: coefficients(:, :)
        complex(dp), allocatable, intent(out) :: derivative(:, :)
        integer, intent(out) :: status
        integer :: mode

        status = 1
        if (size(coefficients, 2) /= size(modes)) return
        if (has_duplicate_modes(modes)) return
        allocate(derivative(size(coefficients, 1), size(coefficients, 2)))
        do mode = 1, size(modes)
            derivative(:, mode) = &
                cmplx(0.0_dp, real(modes(mode), dp), dp)*coefficients(:, mode)
        end do
        status = 0
    end subroutine apply_toroidal_fourier_derivative

    pure integer function find_mode(modes, requested) result(location)
        integer, intent(in) :: modes(:), requested
        integer :: mode

        location = 0
        do mode = 1, size(modes)
            if (modes(mode) /= requested) cycle
            location = mode
            return
        end do
    end function find_mode

    pure logical function has_duplicate_modes(modes) result(duplicate)
        integer, intent(in) :: modes(:)
        integer :: first, second

        duplicate = .false.
        do second = 2, size(modes)
            do first = 1, second - 1
                if (modes(first) /= modes(second)) cycle
                duplicate = .true.
                return
            end do
        end do
    end function has_duplicate_modes

    subroutine assemble_bspline_poloidal_bracket_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            advecting_coefficients, quadrature_order, matrix, status, &
            bracket_weight_function)
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: advecting_coefficients(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        procedure(scalar_weight_2d), optional :: bracket_weight_function

        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient=0.0_dp, &
            mass_coefficient=0.0_dp, &
            advecting_coefficients=advecting_coefficients, &
            advection_coefficient=1.0_dp, &
            advection_weight_function=bracket_weight_function)
    end subroutine assemble_bspline_poloidal_bracket_csc

    subroutine assemble_bspline_toroidal_fourier_laplacian_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, toroidal_mode, matrix, status)
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        integer, intent(in) :: toroidal_mode
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient=1.0_dp, &
            mass_coefficient=1.0_dp, &
            stiffness_weight_function=radial_weight, &
            mass_weight_function=fourier_mass_weight)

    contains

        pure subroutine radial_weight(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = point(1)
        end subroutine radial_weight

        pure subroutine fourier_mass_weight(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = real(toroidal_mode**2, dp)/point(1)
        end subroutine fourier_mass_weight

    end subroutine assemble_bspline_toroidal_fourier_laplacian_csc

    subroutine assemble_bspline_grad_shafranov_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_r(:), knots_z(:)
        integer, intent(in) :: degree_r, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call assemble_bspline_h1_operator_csc( &
            knots_r, knots_z, degree_r, degree_z, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient=1.0_dp, &
            mass_coefficient=0.0_dp, &
            stiffness_weight_function=inverse_radial_weight)

    contains

        pure subroutine inverse_radial_weight(point, value)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value

            value = 1.0_dp/point(1)
        end subroutine inverse_radial_weight

    end subroutine assemble_bspline_grad_shafranov_csc

    subroutine assemble_bspline_l2_hcurl_adjoint_curl_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: weak_curl

        call assemble_bspline_hcurl_l2_curl_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, weak_curl, status)
        if (status%code /= 0) return
        call csc_transpose(weak_curl, matrix, status)
    end subroutine assemble_bspline_l2_hcurl_adjoint_curl_csc

    subroutine assemble_bspline_hcurl_l2_curl_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: curl_incidence, gradient_incidence, l2_mass

        call build_bspline_feec_2d_operators_csc( &
            knots_x, knots_y, degree_x, degree_y, gradient_incidence, &
            curl_incidence, status)
        if (status%code /= 0) return
        call assemble_bspline_l2_mass_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, l2_mass, status)
        if (status%code /= 0) return
        call csc_matmul(l2_mass, curl_incidence, matrix, status)
    end subroutine assemble_bspline_hcurl_l2_curl_csc

    subroutine assemble_bspline_l2_mass_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: qw_x(:), qw_y(:), triplet_values(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: determinant, geometry_jacobian(2, 2), geometry_point(2)
        real(dp) :: inverse(2, 2), physical_weight
        integer :: basis_x, basis_y, entry, local_column, local_count
        integer :: local_row, local_status, max_entries, nx, ny
        integer :: point_x, point_y, span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric L2 mass assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        local_count = degree_x*degree_y
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            qw_x(quadrature_order), qw_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, qw_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, qw_x)
                call build_l2_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y(2:size(knots_y) - 1), degree_y - 1, &
                        nodes_y(point_y), value_y, derivative_y, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x(2:size(knots_x) - 1), degree_x - 1, &
                            nodes_x(point_x), value_x, derivative_x, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = qw_x(point_x)*qw_y(point_y)/determinant
                        do local_column = 1, local_count
                            basis_x = span_x - degree_x + &
                                modulo(local_column - 1, degree_x)
                            basis_y = span_y - degree_y + &
                                (local_column - 1)/degree_x
                            do local_row = 1, local_count
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*value_x(basis_x)* &
                                    value_y(basis_y)*l2_local_value( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, value_y)
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            (nx - 1)*(ny - 1), (nx - 1)*(ny - 1), rows(:entry), &
            columns(:entry), triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_l2_mass_csc

    pure subroutine build_l2_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof

        local_dof = 0
        do basis_y = span_y - degree_y, span_y - 1
            do basis_x = span_x - degree_x, span_x - 1
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*(nx - 1)
            end do
        end do
    end subroutine build_l2_local_dofs

    pure function l2_local_value( &
            local_dof, span_x, span_y, degree_x, degree_y, value_x, value_y) &
            result(value)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: value_x(:), value_y(:)
        real(dp) :: value

        integer :: basis_x, basis_y

        basis_x = span_x - degree_x + modulo(local_dof - 1, degree_x)
        basis_y = span_y - degree_y + (local_dof - 1)/degree_x
        value = value_x(basis_x)*value_y(basis_y)
    end function l2_local_value

    subroutine assemble_bspline_hcurl_h1_adjoint_gradient_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: weak_gradient

        call assemble_bspline_h1_hcurl_gradient_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, weak_gradient, status)
        if (status%code /= 0) return
        call csc_transpose(weak_gradient, matrix, status)
    end subroutine assemble_bspline_hcurl_h1_adjoint_gradient_csc

    subroutine assemble_bspline_h1_hcurl_gradient_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: curl_incidence, gradient_incidence, hcurl_mass

        call build_bspline_feec_2d_operators_csc( &
            knots_x, knots_y, degree_x, degree_y, gradient_incidence, &
            curl_incidence, status)
        if (status%code /= 0) return
        call assemble_bspline_hcurl_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, hcurl_mass, status, curl_coefficient=0.0_dp, &
            mass_coefficient=1.0_dp)
        if (status%code /= 0) return
        call csc_matmul(hcurl_mass, gradient_incidence, matrix, status)
    end subroutine assemble_bspline_h1_hcurl_gradient_csc

    subroutine assemble_bspline_hdiv_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, divergence_coefficient, &
            mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: dx(:), dx_reduced(:), dy(:), dy_reduced(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: qw_x(:), qw_y(:), triplet_values(:)
        real(dp), allocatable :: vx(:), vx_reduced(:), vy(:), vy_reduced(:)
        real(dp) :: determinant, divergence_column, divergence_row
        real(dp) :: divergence_weight, geometry_jacobian(2, 2)
        real(dp) :: geometry_point(2), inverse(2, 2), mass_weight
        real(dp) :: physical_column(2), physical_row(2), physical_weight
        real(dp) :: reference_column(2), reference_row(2)
        integer :: entry, local_column, local_count, local_row, local_status
        integer :: max_entries, nx, ny, point_x, point_y, span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric Hdiv assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        divergence_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (degree_x + 1)*degree_y + degree_x*(degree_y + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            qw_x(quadrature_order), qw_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, qw_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, qw_x)
                call build_hdiv_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), vy, dy, &
                        local_status)
                    if (local_status /= 0) return
                    call evaluate_bspline_basis( &
                        knots_y(2:size(knots_y) - 1), degree_y - 1, &
                        nodes_y(point_y), vy_reduced, dy_reduced, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), vx, dx, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_bspline_basis( &
                            knots_x(2:size(knots_x) - 1), degree_x - 1, &
                            nodes_x(point_x), vx_reduced, dx_reduced, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = determinant*qw_x(point_x)*qw_y(point_y)
                        do local_column = 1, local_count
                            call hdiv_local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, vx, dx, vx_reduced, vy, dy, &
                                vy_reduced, reference_column, divergence_column)
                            physical_column = matmul( &
                                geometry_jacobian, reference_column)/determinant
                            divergence_column = divergence_column/determinant
                            do local_row = 1, local_count
                                call hdiv_local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, vx, dx, vx_reduced, vy, dy, &
                                    vy_reduced, reference_row, divergence_row)
                                physical_row = matmul( &
                                    geometry_jacobian, reference_row)/determinant
                                divergence_row = divergence_row/determinant
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*(mass_weight*dot_product( &
                                    physical_row, physical_column) + &
                                    divergence_weight*divergence_row* &
                                    divergence_column)
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            nx*(ny - 1) + (nx - 1)*ny, &
            nx*(ny - 1) + (nx - 1)*ny, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_hdiv_operator_csc

    pure subroutine build_hdiv_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx, ny
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof, x_component_count

        local_dof = 0
        x_component_count = nx*(ny - 1)
        do basis_y = span_y - degree_y, span_y - 1
            do basis_x = span_x - degree_x, span_x
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*nx
            end do
        end do
        do basis_y = span_y - degree_y, span_y
            do basis_x = span_x - degree_x, span_x - 1
                local_dof = local_dof + 1
                local_dofs(local_dof) = &
                    x_component_count + basis_x + (basis_y - 1)*(nx - 1)
            end do
        end do
    end subroutine build_hdiv_local_dofs

    pure subroutine hdiv_local_basis_data( &
            local_dof, span_x, span_y, degree_x, degree_y, vx, dx, &
            vx_reduced, vy, dy, vy_reduced, value, divergence)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: vx(:), dx(:), vx_reduced(:)
        real(dp), intent(in) :: vy(:), dy(:), vy_reduced(:)
        real(dp), intent(out) :: value(2), divergence

        integer :: basis_x, basis_y, offset, x_local_count

        value = 0.0_dp
        x_local_count = (degree_x + 1)*degree_y
        if (local_dof <= x_local_count) then
            offset = local_dof - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x + 1)
            basis_y = span_y - degree_y + offset/(degree_x + 1)
            value(1) = vx(basis_x)*vy_reduced(basis_y)
            divergence = dx(basis_x)*vy_reduced(basis_y)
        else
            offset = local_dof - x_local_count - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x)
            basis_y = span_y - degree_y + offset/degree_x
            value(2) = vx_reduced(basis_x)*vy(basis_y)
            divergence = vx_reduced(basis_x)*dy(basis_y)
        end if
    end subroutine hdiv_local_basis_data

    subroutine assemble_bspline_hcurl_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, curl_coefficient, &
            mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: dx(:), dx_reduced(:), dy(:), dy_reduced(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: qw_x(:), qw_y(:), triplet_values(:)
        real(dp), allocatable :: vx(:), vx_reduced(:), vy(:), vy_reduced(:)
        real(dp) :: curl_column, curl_row, curl_weight, determinant
        real(dp) :: geometry_jacobian(2, 2), geometry_point(2), inverse(2, 2)
        real(dp) :: mass_weight, physical_column(2), physical_row(2)
        real(dp) :: reference_column(2), reference_row(2), physical_weight
        integer :: entry, local_column, local_count, local_row, local_status
        integer :: max_entries, nx, ny, point_x, point_y, span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric Hcurl assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        curl_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = degree_x*(degree_y + 1) + (degree_x + 1)*degree_y
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            qw_x(quadrature_order), qw_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, qw_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, qw_x)
                call build_hcurl_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), vy, dy, &
                        local_status)
                    if (local_status /= 0) return
                    call evaluate_bspline_basis( &
                        knots_y(2:size(knots_y) - 1), degree_y - 1, &
                        nodes_y(point_y), vy_reduced, dy_reduced, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), vx, dx, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_bspline_basis( &
                            knots_x(2:size(knots_x) - 1), degree_x - 1, &
                            nodes_x(point_x), vx_reduced, dx_reduced, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = determinant*qw_x(point_x)*qw_y(point_y)
                        do local_column = 1, local_count
                            call hcurl_local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, vx, dx, vx_reduced, vy, dy, &
                                vy_reduced, reference_column, curl_column)
                            physical_column = &
                                matmul(transpose(inverse), reference_column)
                            curl_column = curl_column/determinant
                            do local_row = 1, local_count
                                call hcurl_local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, vx, dx, vx_reduced, vy, dy, &
                                    vy_reduced, reference_row, curl_row)
                                physical_row = &
                                    matmul(transpose(inverse), reference_row)
                                curl_row = curl_row/determinant
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*(mass_weight*dot_product( &
                                    physical_row, physical_column) + &
                                    curl_weight*curl_row*curl_column)
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            (nx - 1)*ny + nx*(ny - 1), &
            (nx - 1)*ny + nx*(ny - 1), rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_hcurl_operator_csc

    pure subroutine build_hcurl_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx, ny
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof, x_component_count

        local_dof = 0
        x_component_count = (nx - 1)*ny
        do basis_y = span_y - degree_y, span_y
            do basis_x = span_x - degree_x, span_x - 1
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*(nx - 1)
            end do
        end do
        do basis_y = span_y - degree_y, span_y - 1
            do basis_x = span_x - degree_x, span_x
                local_dof = local_dof + 1
                local_dofs(local_dof) = &
                    x_component_count + basis_x + (basis_y - 1)*nx
            end do
        end do
    end subroutine build_hcurl_local_dofs

    pure subroutine hcurl_local_basis_data( &
            local_dof, span_x, span_y, degree_x, degree_y, vx, dx, &
            vx_reduced, vy, dy, vy_reduced, value, curl)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: vx(:), dx(:), vx_reduced(:)
        real(dp), intent(in) :: vy(:), dy(:), vy_reduced(:)
        real(dp), intent(out) :: value(2), curl

        integer :: basis_x, basis_y, offset, x_local_count

        value = 0.0_dp
        x_local_count = degree_x*(degree_y + 1)
        if (local_dof <= x_local_count) then
            offset = local_dof - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x)
            basis_y = span_y - degree_y + offset/degree_x
            value(1) = vx_reduced(basis_x)*vy(basis_y)
            curl = -vx_reduced(basis_x)*dy(basis_y)
        else
            offset = local_dof - x_local_count - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x + 1)
            basis_y = span_y - degree_y + offset/(degree_x + 1)
            value(2) = vx(basis_x)*vy_reduced(basis_y)
            curl = dx(basis_x)*vy_reduced(basis_y)
        end if
    end subroutine hcurl_local_basis_data

    subroutine build_bspline_feec_2d_operators_csc( &
            knots_x, knots_y, degree_x, degree_y, gradient, curl, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y
        type(csc_t), intent(out) :: gradient, curl
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: derivative_x(:, :), derivative_y(:, :)
        real(dp), allocatable :: values(:)
        integer :: column, entry, ix, iy, nx, ny, row, x_component_count
        integer :: local_status

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric FEEC incidence assembly failed")
        call build_bspline_derivative_matrix( &
            knots_x, degree_x, derivative_x, local_status)
        if (local_status /= 0) return
        call build_bspline_derivative_matrix( &
            knots_y, degree_y, derivative_y, local_status)
        if (local_status /= 0) return
        nx = size(derivative_x, 2)
        ny = size(derivative_y, 2)
        x_component_count = (nx - 1)*ny
        allocate(rows(2*(x_component_count + nx*(ny - 1))))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do iy = 1, ny
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = ix, ix + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = column + (iy - 1)*nx
                    values(entry) = derivative_x(ix, column)
                end do
            end do
        end do
        do iy = 1, ny - 1
            do ix = 1, nx
                row = x_component_count + ix + (iy - 1)*nx
                do column = iy, iy + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = ix + (column - 1)*nx
                    values(entry) = derivative_y(iy, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            x_component_count + nx*(ny - 1), nx*ny, rows(:entry), &
            columns(:entry), values(:entry), gradient, status)
        if (status%code /= 0) return

        deallocate(rows, columns, values)
        allocate(rows(4*(nx - 1)*(ny - 1)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do iy = 1, ny - 1
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = iy, iy + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = ix + (column - 1)*(nx - 1)
                    values(entry) = -derivative_y(iy, column)
                end do
                do column = ix, ix + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = x_component_count + &
                        column + (iy - 1)*nx
                    values(entry) = derivative_x(ix, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            (nx - 1)*(ny - 1), x_component_count + nx*(ny - 1), &
            rows, columns, values, curl, status)
    end subroutine build_bspline_feec_2d_operators_csc

    subroutine assemble_bspline_h1_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient, &
            mass_coefficient, stiffness_weight_function, mass_weight_function, &
            stiffness_tensor_function, advecting_coefficients, &
            advection_coefficient, mass_field_coefficients, &
            advection_weight_function, reaction_coefficients, &
            reaction_direction, reaction_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient
        procedure(scalar_weight_2d), optional :: stiffness_weight_function
        procedure(scalar_weight_2d), optional :: mass_weight_function
        procedure(tensor_weight_2d), optional :: stiffness_tensor_function
        real(dp), intent(in), optional :: advecting_coefficients(:)
        real(dp), intent(in), optional :: advection_coefficient
        real(dp), intent(in), optional :: mass_field_coefficients(:)
        procedure(scalar_weight_2d), optional :: advection_weight_function
        real(dp), intent(in), optional :: reaction_coefficients(:)
        real(dp), intent(in), optional :: reaction_direction(2)
        real(dp), intent(in), optional :: reaction_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: quadrature_weights_x(:)
        real(dp), allocatable :: quadrature_weights_y(:), triplet_values(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: determinant, geometry_jacobian(2, 2), geometry_point(2)
        real(dp) :: basis_column, basis_row, gradient_column(2), gradient_row(2)
        real(dp) :: advecting_gradient(2), advection_velocity(2)
        real(dp) :: advection_weight
        real(dp) :: advection_weight_at_point
        real(dp) :: inverse(2, 2)
        real(dp) :: mass_weight, physical_weight, stiffness_weight
        real(dp) :: mass_field_at_point
        real(dp) :: reaction_field_at_point, reaction_weight
        real(dp) :: mass_weight_at_point, stiffness_weight_at_point
        real(dp) :: stiffness_tensor_at_point(2, 2)
        integer :: entry, local_column, local_count, local_row
        integer :: local_status, max_entries, nx, ny, point_x, point_y
        integer :: span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric H1 assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        advection_weight = 0.0_dp
        if (present(advection_coefficient)) then
            advection_weight = advection_coefficient
        end if
        if (present(advecting_coefficients)) then
            if (size(advecting_coefficients) /= nx*ny) return
        else if (advection_weight /= 0.0_dp) then
            return
        end if
        if (present(mass_field_coefficients)) then
            if (size(mass_field_coefficients) /= nx*ny) return
        end if
        reaction_weight = 0.0_dp
        if (present(reaction_coefficient)) then
            reaction_weight = reaction_coefficient
        end if
        if (present(reaction_coefficients)) then
            if (size(reaction_coefficients) /= nx*ny) return
            if (.not. present(reaction_direction)) return
        else if (reaction_weight /= 0.0_dp) then
            return
        end if
        local_count = (degree_x + 1)*(degree_y + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), &
            local_dofs(local_count), local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            quadrature_weights_x(quadrature_order), &
            quadrature_weights_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, quadrature_weights_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, quadrature_weights_x)
                call build_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), value_y, &
                        derivative_y, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), value_x, &
                            derivative_x, local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = determinant* &
                            quadrature_weights_x(point_x)* &
                            quadrature_weights_y(point_y)
                        stiffness_weight_at_point = stiffness_weight
                        mass_weight_at_point = mass_weight
                        advection_weight_at_point = advection_weight
                        stiffness_tensor_at_point = reshape( &
                            [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
                        if (present(stiffness_weight_function)) then
                            call stiffness_weight_function( &
                                geometry_point, stiffness_weight_at_point)
                        end if
                        if (present(mass_weight_function)) then
                            call mass_weight_function( &
                                geometry_point, mass_weight_at_point)
                        end if
                        mass_field_at_point = 1.0_dp
                        if (present(mass_field_coefficients)) then
                            mass_field_at_point = 0.0_dp
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                mass_field_at_point = mass_field_at_point + &
                                    mass_field_coefficients( &
                                    local_dofs(local_row))*basis_row
                            end do
                        end if
                        mass_weight_at_point = &
                            mass_weight_at_point*mass_field_at_point
                        reaction_field_at_point = 0.0_dp
                        if (present(reaction_coefficients)) then
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                reaction_field_at_point = &
                                    reaction_field_at_point + &
                                    reaction_coefficients( &
                                    local_dofs(local_row))*basis_row
                            end do
                        end if
                        if (present(stiffness_tensor_function)) then
                            call stiffness_tensor_function( &
                                geometry_point, stiffness_tensor_at_point)
                        end if
                        if (present(advection_weight_function)) then
                            call advection_weight_function( &
                                geometry_point, advection_weight_at_point)
                        end if
                        advecting_gradient = 0.0_dp
                        if (present(advecting_coefficients)) then
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                gradient_row = matmul( &
                                    transpose(inverse), gradient_row)
                                advecting_gradient = advecting_gradient + &
                                    advecting_coefficients( &
                                    local_dofs(local_row))*gradient_row
                            end do
                        end if
                        advection_velocity = [ &
                            -advecting_gradient(2), advecting_gradient(1)]
                        do local_column = 1, local_count
                            call local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, value_x, derivative_x, value_y, &
                                derivative_y, basis_column, gradient_column)
                            gradient_column = &
                                matmul(transpose(inverse), gradient_column)
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                gradient_row = &
                                    matmul(transpose(inverse), gradient_row)
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*(stiffness_weight_at_point* &
                                    dot_product(gradient_row, matmul( &
                                    stiffness_tensor_at_point, &
                                    gradient_column)) + &
                                    mass_weight_at_point*basis_row*basis_column)
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    0.5_dp*physical_weight* &
                                    advection_weight_at_point*( &
                                    basis_row*dot_product( &
                                    advection_velocity, gradient_column) - &
                                    basis_column*dot_product( &
                                    advection_velocity, gradient_row))
                                if (present(reaction_direction)) then
                                    local_matrix(local_row, local_column) = &
                                        local_matrix(local_row, local_column) + &
                                        physical_weight*reaction_weight* &
                                        reaction_field_at_point*basis_row* &
                                        dot_product( &
                                        reaction_direction, gradient_column)
                                end if
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            nx*ny, nx*ny, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_h1_operator_csc

    subroutine assemble_bspline_h1_operator_csc_jvp( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            control_points_dot, weights_dot, quadrature_order, matrix_dot, &
            status, stiffness_coefficient, mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: control_points_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:, :)
        type(csc_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: local_matrix_dot(:, :)
        real(dp), allocatable :: nodes_x(:), nodes_y(:)
        real(dp), allocatable :: quadrature_weights_x(:)
        real(dp), allocatable :: quadrature_weights_y(:)
        real(dp), allocatable :: triplet_values(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: basis_column, basis_row, contribution_dot
        real(dp) :: geometry_jacobian(2, 2), geometry_jacobian_dot(2, 2)
        real(dp) :: geometry_point(2), geometry_point_dot(2)
        real(dp) :: gradient_column(2), gradient_row(2)
        real(dp) :: mass_weight, quadrature_weight, stiffness_weight
        integer :: entry, local_column, local_count, local_row
        integer :: local_status, max_entries, nx, ny, point_x, point_y
        integer :: span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric H1 assembly JVP failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2) return
        if (any(shape(control_points(1, :, :)) /= [nx, ny])) return
        if (any(shape(control_points_dot) /= shape(control_points))) return
        if (any(shape(weights_dot) /= shape(weights))) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (degree_x + 1)*(degree_y + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate ( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix_dot(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            quadrature_weights_x(quadrature_order), &
            quadrature_weights_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, quadrature_weights_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, quadrature_weights_x)
                call build_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, local_dofs)
                local_matrix_dot = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), value_y, &
                        derivative_y, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), value_x, &
                            derivative_x, local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry_jvp( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, control_points_dot, &
                            weights_dot, nodes_x(point_x), nodes_y(point_y), &
                            geometry_point_dot, geometry_jacobian_dot, &
                            local_status)
                        if (local_status /= 0) return
                        quadrature_weight = quadrature_weights_x(point_x)* &
                            quadrature_weights_y(point_y)
                        do local_column = 1, local_count
                            call local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, value_x, derivative_x, value_y, &
                                derivative_y, basis_column, gradient_column)
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                call generated_bspline_h1_geometry_jvp( &
                                    geometry_jacobian(1, 1), &
                                    geometry_jacobian(2, 1), &
                                    geometry_jacobian(1, 2), &
                                    geometry_jacobian(2, 2), gradient_row(1), &
                                    gradient_row(2), gradient_column(1), &
                                    gradient_column(2), basis_row, basis_column, &
                                    stiffness_weight, mass_weight, &
                                    quadrature_weight, &
                                    geometry_jacobian_dot(1, 1), &
                                    geometry_jacobian_dot(2, 1), &
                                    geometry_jacobian_dot(1, 2), &
                                    geometry_jacobian_dot(2, 2), &
                                    contribution_dot)
                                local_matrix_dot(local_row, local_column) = &
                                    local_matrix_dot(local_row, local_column) + &
                                    contribution_dot
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix_dot(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            nx*ny, nx*ny, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix_dot, status)
    end subroutine assemble_bspline_h1_operator_csc_jvp

    subroutine assemble_bspline_h1_operator_csc_vjp( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix_bar, control_points_bar, weights_bar, &
            status, stiffness_coefficient, mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(in) :: matrix_bar
        real(dp), intent(out) :: control_points_bar(:, :, :)
        real(dp), intent(out) :: weights_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: local_dofs(:)
        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: geometry_control_points_bar(:, :, :)
        real(dp), allocatable :: geometry_weights_bar(:, :)
        real(dp), allocatable :: nodes_x(:), nodes_y(:)
        real(dp), allocatable :: quadrature_weights_x(:)
        real(dp), allocatable :: quadrature_weights_y(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: basis_column, basis_row, contribution_bar
        real(dp) :: geometry_jacobian(2, 2), geometry_jacobian_bar(2, 2)
        real(dp) :: geometry_jacobian_bar_entry(2, 2)
        real(dp) :: geometry_point(2), geometry_point_bar(2)
        real(dp) :: gradient_column(2), gradient_row(2)
        real(dp) :: mass_weight, quadrature_weight, stiffness_weight
        integer :: local_column, local_count, local_row
        integer :: local_status, nx, ny, point_x, point_y, span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric H1 assembly VJP failed")
        control_points_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2) return
        if (any(shape(control_points(1, :, :)) /= [nx, ny])) return
        if (any(shape(control_points_bar) /= shape(control_points))) return
        if (any(shape(weights_bar) /= shape(weights))) return
        if (.not. csc_is_valid(matrix_bar)) return
        if (matrix_bar%nrow /= nx*ny .or. matrix_bar%ncol /= nx*ny) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (degree_x + 1)*(degree_y + 1)
        allocate ( &
            local_dofs(local_count), nodes_x(quadrature_order), &
            nodes_y(quadrature_order), &
            quadrature_weights_x(quadrature_order), &
            quadrature_weights_y(quadrature_order))
        allocate (geometry_control_points_bar, mold=control_points)
        allocate (geometry_weights_bar, mold=weights)
        geometry_point_bar = 0.0_dp
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, quadrature_weights_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, quadrature_weights_x)
                call build_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, local_dofs)
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), value_y, &
                        derivative_y, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), value_x, &
                            derivative_x, local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        quadrature_weight = quadrature_weights_x(point_x)* &
                            quadrature_weights_y(point_y)
                        geometry_jacobian_bar = 0.0_dp
                        do local_column = 1, local_count
                            call local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, value_x, derivative_x, value_y, &
                                derivative_y, basis_column, gradient_column)
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                contribution_bar = csc_entry_value( &
                                    matrix_bar, local_dofs(local_row), &
                                    local_dofs(local_column))
                                call generated_bspline_h1_geometry_vjp( &
                                    geometry_jacobian(1, 1), &
                                    geometry_jacobian(2, 1), &
                                    geometry_jacobian(1, 2), &
                                    geometry_jacobian(2, 2), gradient_row(1), &
                                    gradient_row(2), gradient_column(1), &
                                    gradient_column(2), basis_row, basis_column, &
                                    stiffness_weight, mass_weight, &
                                    quadrature_weight, contribution_bar, &
                                    geometry_jacobian_bar_entry(1, 1), &
                                    geometry_jacobian_bar_entry(2, 1), &
                                    geometry_jacobian_bar_entry(1, 2), &
                                    geometry_jacobian_bar_entry(2, 2))
                                geometry_jacobian_bar = &
                                    geometry_jacobian_bar + &
                                    geometry_jacobian_bar_entry
                            end do
                        end do
                        call evaluate_nurbs_surface_geometry_vjp( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point_bar, &
                            geometry_jacobian_bar, &
                            geometry_control_points_bar, geometry_weights_bar, &
                            local_status)
                        if (local_status /= 0) return
                        control_points_bar = control_points_bar + &
                            geometry_control_points_bar
                        weights_bar = weights_bar + geometry_weights_bar
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_bspline_h1_operator_csc_vjp

    subroutine local_basis_data( &
            local_dof, span_x, span_y, degree_x, degree_y, value_x, &
            derivative_x, value_y, derivative_y, basis_value, gradient)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: value_x(:), derivative_x(:)
        real(dp), intent(in) :: value_y(:), derivative_y(:)
        real(dp), intent(out) :: basis_value
        real(dp), intent(out) :: gradient(2)

        integer :: basis_x, basis_y, offset_x, offset_y

        offset_x = modulo(local_dof - 1, degree_x + 1)
        offset_y = (local_dof - 1)/(degree_x + 1)
        basis_x = span_x - degree_x + offset_x
        basis_y = span_y - degree_y + offset_y
        basis_value = value_x(basis_x)*value_y(basis_y)
        gradient = [ &
            derivative_x(basis_x)*value_y(basis_y), &
            value_x(basis_x)*derivative_y(basis_y)]
    end subroutine local_basis_data

    pure subroutine build_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof

        local_dof = 0
        do basis_y = span_y - degree_y, span_y
            do basis_x = span_x - degree_x, span_x
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*nx
            end do
        end do
    end subroutine build_local_dofs

    pure integer function positive_span_count( &
            knots, degree, basis_count) result(count)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree, basis_count
        integer :: span

        count = 0
        do span = degree + 1, basis_count
            if (knots(span + 1) > knots(span)) count = count + 1
        end do
    end function positive_span_count

    pure real(dp) function csc_entry_value(matrix, row, column) result(value)
        type(csc_t), intent(in) :: matrix
        integer, intent(in) :: row, column

        integer :: entry

        value = 0.0_dp
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                value = matrix%val(entry)
                return
            end if
        end do
    end function csc_entry_value

    pure subroutine inverse_2d(matrix, inverse, determinant, status)
        real(dp), intent(in) :: matrix(2, 2)
        real(dp), intent(out) :: inverse(2, 2), determinant
        integer, intent(out) :: status

        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        status = 1
        inverse = 0.0_dp
        if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(matrix))**2)) return
        inverse = reshape([ &
            matrix(2, 2), -matrix(2, 1), &
            -matrix(1, 2), matrix(1, 1)], [2, 2])/determinant
        status = 0
    end subroutine inverse_2d

end module fortfem_assembly_bspline_2d
