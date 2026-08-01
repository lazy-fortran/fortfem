module fortfem_generalized_eigen_residual
    !! Neutral complex generalized eigen-residual and derivative actions.
    !!
    !! The contract is deliberately independent of any eigensolver or
    !! application convention.  It exposes the residual
    !!
    !!     r(K,M,lambda,u) = K u - lambda M u
    !!
    !! and its analytic tangent and real-complex adjoint actions.  Consumers
    !! such as GLISS, GPEC, or MARS-style adapters can therefore provide their
    !! own assembled operators without importing application-specific formats.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: assemble_generalized_eigen_residual
    public :: assemble_generalized_eigen_residual_jvp
    public :: assemble_generalized_eigen_residual_vjp

    interface finite_complex
        module procedure finite_complex_scalar
        module procedure finite_complex_vector
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine assemble_generalized_eigen_residual( &
            stiffness, mass, eigenvalue, state, residual, status)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :)
        complex(dp), intent(in) :: eigenvalue, state(:)
        complex(dp), intent(out) :: residual(:)
        integer, intent(out) :: status

        residual = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_residual_shapes(stiffness, mass, state, residual) .or. &
            .not. finite_complex(stiffness) .or. .not. finite_complex(mass) .or. &
            .not. finite_complex(eigenvalue) .or. .not. finite_complex(state)) return

        residual = matmul(stiffness, state) - &
            eigenvalue*matmul(mass, state)
        status = 0
    end subroutine assemble_generalized_eigen_residual

    subroutine assemble_generalized_eigen_residual_jvp( &
            stiffness, mass, eigenvalue, state, stiffness_dot, mass_dot, &
            eigenvalue_dot, state_dot, residual_dot, status)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :)
        complex(dp), intent(in) :: eigenvalue, state(:)
        complex(dp), intent(in) :: stiffness_dot(:, :), mass_dot(:, :)
        complex(dp), intent(in) :: eigenvalue_dot, state_dot(:)
        complex(dp), intent(out) :: residual_dot(:)
        integer, intent(out) :: status

        residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_residual_shapes( &
            stiffness, mass, state, residual_dot) .or. &
            .not. valid_tangent_shapes( &
            stiffness, mass, state, stiffness_dot, mass_dot, state_dot) .or. &
            .not. finite_complex(stiffness) .or. .not. finite_complex(mass) .or. &
            .not. finite_complex(eigenvalue) .or. .not. finite_complex(state) .or. &
            .not. finite_complex(stiffness_dot) .or. &
            .not. finite_complex(mass_dot) .or. &
            .not. finite_complex(eigenvalue_dot) .or. &
            .not. finite_complex(state_dot)) return

        residual_dot = matmul(stiffness_dot, state) + &
            matmul(stiffness, state_dot) - eigenvalue_dot*matmul(mass, state) - &
            eigenvalue*matmul(mass_dot, state) - &
            eigenvalue*matmul(mass, state_dot)
        status = 0
    end subroutine assemble_generalized_eigen_residual_jvp

    subroutine assemble_generalized_eigen_residual_vjp( &
            stiffness, mass, eigenvalue, state, residual_bar, stiffness_bar, &
            mass_bar, state_bar, eigenvalue_bar, status)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :)
        complex(dp), intent(in) :: eigenvalue, state(:)
        complex(dp), intent(in) :: residual_bar(:)
        complex(dp), intent(out) :: stiffness_bar(:, :), mass_bar(:, :)
        complex(dp), intent(out) :: state_bar(:), eigenvalue_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: mass_state(:)

        stiffness_bar = cmplx(0.0_dp, 0.0_dp, dp)
        mass_bar = cmplx(0.0_dp, 0.0_dp, dp)
        state_bar = cmplx(0.0_dp, 0.0_dp, dp)
        eigenvalue_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_vjp_shapes( &
            stiffness, mass, state, residual_bar, stiffness_bar, mass_bar, &
            state_bar) .or. &
            .not. finite_complex(stiffness) .or. &
            .not. finite_complex(mass) .or. .not. finite_complex(eigenvalue) .or. &
            .not. finite_complex(state) .or. .not. finite_complex(residual_bar)) return

        allocate(mass_state(size(state)))
        mass_state = matmul(mass, state)
        stiffness_bar = outer_product(residual_bar, conjg(state))
        mass_bar = -conjg(eigenvalue)*outer_product(residual_bar, conjg(state))
        state_bar = matmul(conjg(transpose(stiffness)), residual_bar) - &
            conjg(eigenvalue)*matmul(conjg(transpose(mass)), residual_bar)
        eigenvalue_bar = -sum(residual_bar*conjg(mass_state))
        status = 0
    end subroutine assemble_generalized_eigen_residual_vjp

    logical function valid_residual_shapes(stiffness, mass, state, residual) &
            result(valid)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), state(:), residual(:)

        valid = size(stiffness, 1) > 0 .and. &
            size(stiffness, 1) == size(stiffness, 2) .and. &
            all(shape(mass) == shape(stiffness)) .and. &
            size(state) == size(stiffness, 1) .and. &
            size(residual) == size(state)
    end function valid_residual_shapes

    logical function valid_tangent_shapes( &
            stiffness, mass, state, stiffness_dot, mass_dot, state_dot) &
            result(valid)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), state(:)
        complex(dp), intent(in) :: stiffness_dot(:, :), mass_dot(:, :), state_dot(:)

        valid = all(shape(stiffness_dot) == shape(stiffness)) .and. &
            all(shape(mass_dot) == shape(mass)) .and. size(state_dot) == size(state)
    end function valid_tangent_shapes

    logical function valid_vjp_shapes( &
            stiffness, mass, state, residual_bar, stiffness_bar, mass_bar, state_bar) &
            result(valid)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), state(:)
        complex(dp), intent(in) :: residual_bar(:)
        complex(dp), intent(in) :: stiffness_bar(:, :), mass_bar(:, :), state_bar(:)

        valid = size(stiffness, 1) > 0 .and. &
            size(stiffness, 1) == size(stiffness, 2) .and. &
            all(shape(mass) == shape(stiffness)) .and. &
            size(state) == size(mass, 1) .and. size(residual_bar) == size(state) .and. &
            all(shape(stiffness_bar) == shape(mass)) .and. &
            all(shape(mass_bar) == shape(mass)) .and. size(state_bar) == size(state)
    end function valid_vjp_shapes

    pure function outer_product(left, right) result(product)
        complex(dp), intent(in) :: left(:), right(:)
        complex(dp) :: product(size(left), size(right))
        integer :: i

        do i = 1, size(left)
            product(i, :) = left(i)*right
        end do
    end function outer_product

    logical function finite_complex_scalar(value) result(valid)
        complex(dp), intent(in) :: value

        valid = ieee_is_finite(real(value, dp)) .and. ieee_is_finite(aimag(value))
    end function finite_complex_scalar

    logical function finite_complex_vector(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_vector

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

end module fortfem_generalized_eigen_residual
