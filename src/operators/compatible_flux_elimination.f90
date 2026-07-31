module fortfem_compatible_flux_elimination
    !! Flux-specific Schur reduction with an H(div)-compatible recovery map.
    !!
    !! Given
    !!
    !!   M q + F u = f,     G q + D u = g,
    !!
    !! the primitive returns r=M^(-1)f, X=-M^(-1)F, and
    !! S=D+G X, b=g-G r.  A caller can recover the local flux as
    !! q = r + X u while solving only the condensed state system S u=b.
    !! The flux space, divergence/trace maps, and global skeleton ownership
    !! remain caller-owned, so the same contract applies to RT, BDM, H(curl),
    !! and compatible IGA blocks.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    public :: assemble_compatible_flux_elimination
    public :: assemble_compatible_flux_elimination_jvp
    public :: assemble_compatible_flux_elimination_vjp

contains

    subroutine assemble_compatible_flux_elimination( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
            status)
        real(dp), intent(in) :: flux_mass(:, :), flux_to_state(:, :)
        real(dp), intent(in) :: state_to_flux(:, :), state_matrix(:, :)
        real(dp), intent(in) :: flux_rhs(:), state_rhs(:)
        real(dp), intent(out) :: recovery(:), recovery_matrix(:, :)
        real(dp), intent(out) :: condensed_matrix(:, :), condensed_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: flux_count, state_count, info
        real(dp), allocatable :: flux_solution(:, :), recovery_rhs(:)

        recovery = 0.0_dp
        recovery_matrix = 0.0_dp
        condensed_matrix = 0.0_dp
        condensed_rhs = 0.0_dp
        if (.not. validate_inputs( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
            status)) return
        flux_count = size(flux_mass, 1)
        state_count = size(state_matrix, 1)
        allocate(flux_solution(flux_count, state_count), recovery_rhs(flux_count))
        call dense_solve(flux_mass, flux_to_state, flux_solution, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux mass block is singular")
            return
        end if
        call dense_solve(flux_mass, flux_rhs, recovery_rhs, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux mass block is singular")
            return
        end if
        recovery = recovery_rhs
        recovery_matrix = -flux_solution
        condensed_matrix = state_matrix + matmul( &
            state_to_flux, recovery_matrix)
        condensed_rhs = state_rhs - matmul(state_to_flux, recovery)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_compatible_flux_elimination

    subroutine assemble_compatible_flux_elimination_jvp( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, flux_mass_dot, flux_to_state_dot, state_to_flux_dot, &
            state_matrix_dot, flux_rhs_dot, state_rhs_dot, recovery_dot, &
            recovery_matrix_dot, condensed_matrix_dot, condensed_rhs_dot, status)
        real(dp), intent(in) :: flux_mass(:, :), flux_to_state(:, :)
        real(dp), intent(in) :: state_to_flux(:, :), state_matrix(:, :)
        real(dp), intent(in) :: flux_rhs(:), state_rhs(:)
        real(dp), intent(in) :: flux_mass_dot(:, :), flux_to_state_dot(:, :)
        real(dp), intent(in) :: state_to_flux_dot(:, :), state_matrix_dot(:, :)
        real(dp), intent(in) :: flux_rhs_dot(:), state_rhs_dot(:)
        real(dp), intent(out) :: recovery_dot(:), recovery_matrix_dot(:, :)
        real(dp), intent(out) :: condensed_matrix_dot(:, :), condensed_rhs_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: flux_count, state_count, info
        real(dp), allocatable :: flux_solution(:, :), recovery(:)
        real(dp), allocatable :: flux_solution_dot(:, :), recovery_dot_local(:)

        recovery_dot = 0.0_dp
        recovery_matrix_dot = 0.0_dp
        condensed_matrix_dot = 0.0_dp
        condensed_rhs_dot = 0.0_dp
        if (.not. validate_base_for_jvp( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery_dot, recovery_matrix_dot, condensed_matrix_dot, &
            condensed_rhs_dot, status)) return
        flux_count = size(flux_mass, 1)
        state_count = size(state_matrix, 1)
        if (.not. validate_direction( &
            flux_mass_dot, flux_to_state_dot, state_to_flux_dot, state_matrix_dot, &
            flux_rhs_dot, state_rhs_dot, flux_count, state_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "compatible flux elimination JVP has incompatible increments")
            return
        end if
        allocate(flux_solution(flux_count, state_count), recovery(flux_count), &
            flux_solution_dot(flux_count, state_count), recovery_dot_local(flux_count))
        call dense_solve(flux_mass, flux_to_state, flux_solution, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux mass block is singular")
            return
        end if
        call dense_solve(flux_mass, flux_rhs, recovery, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux mass block is singular")
            return
        end if
        call dense_solve(flux_mass, flux_to_state_dot - matmul( &
            flux_mass_dot, flux_solution), flux_solution_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux elimination JVP system is singular")
            return
        end if
        call dense_solve(flux_mass, flux_rhs_dot - matmul( &
            flux_mass_dot, recovery), recovery_dot_local, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux elimination JVP system is singular")
            return
        end if
        recovery_dot = recovery_dot_local
        recovery_matrix_dot = -flux_solution_dot
        condensed_matrix_dot = state_matrix_dot + matmul( &
            state_to_flux_dot, -flux_solution) + matmul( &
            state_to_flux, recovery_matrix_dot)
        condensed_rhs_dot = state_rhs_dot - matmul( &
            state_to_flux_dot, recovery) - matmul( &
            state_to_flux, recovery_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_compatible_flux_elimination_jvp

    subroutine assemble_compatible_flux_elimination_vjp( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
            recovery_bar, recovery_matrix_bar, condensed_matrix_bar, &
            condensed_rhs_bar, flux_mass_bar, flux_to_state_bar, state_to_flux_bar, &
            state_matrix_bar, flux_rhs_bar, state_rhs_bar, status)
        real(dp), intent(in) :: flux_mass(:, :), flux_to_state(:, :)
        real(dp), intent(in) :: state_to_flux(:, :), state_matrix(:, :)
        real(dp), intent(in) :: flux_rhs(:), state_rhs(:)
        real(dp), intent(in) :: recovery(:), recovery_matrix(:, :)
        real(dp), intent(in) :: condensed_matrix(:, :), condensed_rhs(:)
        real(dp), intent(in) :: recovery_bar(:), recovery_matrix_bar(:, :)
        real(dp), intent(in) :: condensed_matrix_bar(:, :), condensed_rhs_bar(:)
        real(dp), intent(out) :: flux_mass_bar(:, :), flux_to_state_bar(:, :)
        real(dp), intent(out) :: state_to_flux_bar(:, :), state_matrix_bar(:, :)
        real(dp), intent(out) :: flux_rhs_bar(:), state_rhs_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: flux_count, state_count, info, column
        real(dp), allocatable :: flux_solution(:, :), recovery_local(:)
        real(dp), allocatable :: recovery_bar_local(:), solution_bar(:, :)
        real(dp), allocatable :: transpose_bar(:), rhs_bar(:)

        flux_mass_bar = 0.0_dp
        flux_to_state_bar = 0.0_dp
        state_to_flux_bar = 0.0_dp
        state_matrix_bar = 0.0_dp
        flux_rhs_bar = 0.0_dp
        state_rhs_bar = 0.0_dp
        if (.not. validate_inputs( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
            status)) return
        flux_count = size(flux_mass, 1)
        state_count = size(state_matrix, 1)
        if (.not. validate_vjp_outputs( &
            recovery_bar, recovery_matrix_bar, condensed_matrix_bar, &
            condensed_rhs_bar, flux_count, state_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "compatible flux elimination VJP has incompatible cotangents")
            return
        end if
        allocate(flux_solution(flux_count, state_count), recovery_local(flux_count), &
            recovery_bar_local(flux_count), solution_bar(flux_count, state_count), &
            transpose_bar(flux_count), rhs_bar(flux_count))
        call dense_solve(flux_mass, flux_to_state, flux_solution, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux mass block is singular")
            return
        end if
        call dense_solve(flux_mass, flux_rhs, recovery_local, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux mass block is singular")
            return
        end if
        state_matrix_bar = condensed_matrix_bar
        state_rhs_bar = condensed_rhs_bar
        state_to_flux_bar = matmul(condensed_matrix_bar, transpose(-flux_solution)) - &
            outer_product(condensed_rhs_bar, recovery_local)
        recovery_bar_local = recovery_bar - matmul( &
            transpose(state_to_flux), condensed_rhs_bar)
        solution_bar = -recovery_matrix_bar - matmul( &
            transpose(state_to_flux), condensed_matrix_bar)
        do column = 1, state_count
            call dense_solve(transpose(flux_mass), solution_bar(:, column), &
                transpose_bar, info)
            if (info == 0) then
                flux_to_state_bar(:, column) = transpose_bar
                flux_mass_bar = flux_mass_bar - outer_product( &
                    transpose_bar, flux_solution(:, column))
            end if
        end do
        call dense_solve(transpose(flux_mass), recovery_bar_local, rhs_bar, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "compatible flux elimination VJP system is singular")
            return
        end if
        flux_rhs_bar = rhs_bar
        flux_mass_bar = flux_mass_bar - outer_product(rhs_bar, recovery_local)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_compatible_flux_elimination_vjp

    logical function validate_inputs( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
            status) result(valid)
        real(dp), intent(in) :: flux_mass(:, :), flux_to_state(:, :)
        real(dp), intent(in) :: state_to_flux(:, :), state_matrix(:, :)
        real(dp), intent(in) :: flux_rhs(:), state_rhs(:)
        real(dp), intent(in) :: recovery(:), recovery_matrix(:, :)
        real(dp), intent(in) :: condensed_matrix(:, :), condensed_rhs(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: flux_count, state_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "compatible flux elimination has incompatible arrays")
        flux_count = size(flux_mass, 1)
        state_count = size(state_matrix, 1)
        if (flux_count < 1 .or. state_count < 1 .or. &
            size(flux_mass, 2) /= flux_count .or. &
            size(flux_to_state, 1) /= flux_count .or. &
            size(flux_to_state, 2) /= state_count .or. &
            size(state_to_flux, 1) /= state_count .or. &
            size(state_to_flux, 2) /= flux_count .or. &
            size(state_matrix, 2) /= state_count .or. &
            size(flux_rhs) /= flux_count .or. size(state_rhs) /= state_count .or. &
            size(recovery) /= flux_count .or. &
            size(recovery_matrix, 1) /= flux_count .or. &
            size(recovery_matrix, 2) /= state_count .or. &
            size(condensed_matrix, 1) /= state_count .or. &
            size(condensed_matrix, 2) /= state_count .or. &
            size(condensed_rhs) /= state_count) return
        if (.not. finite_matrix(flux_mass) .or. &
            .not. finite_matrix(flux_to_state) .or. &
            .not. finite_matrix(state_to_flux) .or. &
            .not. finite_matrix(state_matrix) .or. &
            .not. finite_vector(flux_rhs) .or. .not. finite_vector(state_rhs) .or. &
            .not. finite_vector(recovery) .or. &
            .not. finite_matrix(recovery_matrix) .or. &
            .not. finite_matrix(condensed_matrix) .or. &
            .not. finite_vector(condensed_rhs)) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_inputs

    logical function validate_base_for_jvp( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
            status) result(valid)
        real(dp), intent(in) :: flux_mass(:, :), flux_to_state(:, :)
        real(dp), intent(in) :: state_to_flux(:, :), state_matrix(:, :)
        real(dp), intent(in) :: flux_rhs(:), state_rhs(:)
        real(dp), intent(in) :: recovery(:), recovery_matrix(:, :)
        real(dp), intent(in) :: condensed_matrix(:, :), condensed_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        valid = validate_inputs( &
            flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
            state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
            status)
    end function validate_base_for_jvp

    logical function validate_direction( &
            flux_mass_dot, flux_to_state_dot, state_to_flux_dot, state_matrix_dot, &
            flux_rhs_dot, state_rhs_dot, flux_count, state_count) result(valid)
        real(dp), intent(in) :: flux_mass_dot(:, :), flux_to_state_dot(:, :)
        real(dp), intent(in) :: state_to_flux_dot(:, :), state_matrix_dot(:, :)
        real(dp), intent(in) :: flux_rhs_dot(:), state_rhs_dot(:)
        integer, intent(in) :: flux_count, state_count

        valid = size(flux_mass_dot, 1) == flux_count .and. &
            size(flux_mass_dot, 2) == flux_count .and. &
            size(flux_to_state_dot, 1) == flux_count .and. &
            size(flux_to_state_dot, 2) == state_count .and. &
            size(state_to_flux_dot, 1) == state_count .and. &
            size(state_to_flux_dot, 2) == flux_count .and. &
            size(state_matrix_dot, 1) == state_count .and. &
            size(state_matrix_dot, 2) == state_count .and. &
            size(flux_rhs_dot) == flux_count .and. &
            size(state_rhs_dot) == state_count .and. &
            finite_matrix(flux_mass_dot) .and. finite_matrix(flux_to_state_dot) .and. &
            finite_matrix(state_to_flux_dot) .and. &
            finite_matrix(state_matrix_dot) .and. &
            finite_vector(flux_rhs_dot) .and. finite_vector(state_rhs_dot)
    end function validate_direction

    logical function validate_vjp_outputs( &
            recovery_bar, recovery_matrix_bar, condensed_matrix_bar, &
            condensed_rhs_bar, flux_count, state_count) result(valid)
        real(dp), intent(in) :: recovery_bar(:), recovery_matrix_bar(:, :)
        real(dp), intent(in) :: condensed_matrix_bar(:, :), condensed_rhs_bar(:)
        integer, intent(in) :: flux_count, state_count

        valid = size(recovery_bar) == flux_count .and. &
            size(recovery_matrix_bar, 1) == flux_count .and. &
            size(recovery_matrix_bar, 2) == state_count .and. &
            size(condensed_matrix_bar, 1) == state_count .and. &
            size(condensed_matrix_bar, 2) == state_count .and. &
            size(condensed_rhs_bar) == state_count .and. &
            finite_vector(recovery_bar) .and. finite_matrix(recovery_matrix_bar) .and. &
            finite_matrix(condensed_matrix_bar) .and. finite_vector(condensed_rhs_bar)
    end function validate_vjp_outputs

    pure logical function finite_matrix(values) result(valid)
        real(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(values))
    end function finite_matrix

    pure logical function finite_vector(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_vector

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: product(size(left), size(right))

        product = spread(left, dim=2, ncopies=size(right))* &
            spread(right, dim=1, ncopies=size(left))
    end function outer_product

end module fortfem_compatible_flux_elimination
