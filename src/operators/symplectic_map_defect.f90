module fortfem_symplectic_map_defect
    !! Structure diagnostic for a fixed linear one-step map.
    !!
    !! For a state map S and a caller-owned skew form Omega, the residual
    !! S^T Omega S - Omega reports failure to preserve the declared
    !! symplectic form.  The block is neutral: it does not choose canonical
    !! coordinates, a time integrator, or a PDE.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_symplectic_map_defect
    public :: assemble_symplectic_map_defect_jvp
    public :: assemble_symplectic_map_defect_vjp

contains

    subroutine assemble_symplectic_map_defect( &
            map, symplectic_form, defect, status)
        real(dp), intent(in) :: map(:, :), symplectic_form(:, :)
        real(dp), intent(out) :: defect(:, :)
        type(fortsparse_status_t), intent(out) :: status

        defect = 0.0_dp
        call validate_inputs(map, symplectic_form, defect, status)
        if (status%code /= FORTSPARSE_OK) return
        defect = matmul(transpose(map), matmul(symplectic_form, map)) - &
            symplectic_form
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_symplectic_map_defect

    subroutine assemble_symplectic_map_defect_jvp( &
            map, symplectic_form, map_dot, symplectic_form_dot, defect_dot, status)
        real(dp), intent(in) :: map(:, :), symplectic_form(:, :)
        real(dp), intent(in) :: map_dot(:, :), symplectic_form_dot(:, :)
        real(dp), intent(out) :: defect_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        defect_dot = 0.0_dp
        call validate_inputs(map, symplectic_form, defect_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_direction( &
            map, symplectic_form, map_dot, symplectic_form_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "symplectic-map JVP has incompatible increments")
            return
        end if
        defect_dot = matmul(transpose(map_dot), matmul(symplectic_form, map)) + &
            matmul(transpose(map), matmul(symplectic_form_dot, map)) + &
            matmul(transpose(map), matmul(symplectic_form, map_dot)) - &
            symplectic_form_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_symplectic_map_defect_jvp

    subroutine assemble_symplectic_map_defect_vjp( &
            map, symplectic_form, defect_bar, map_bar, symplectic_form_bar, status)
        real(dp), intent(in) :: map(:, :), symplectic_form(:, :), defect_bar(:, :)
        real(dp), intent(out) :: map_bar(:, :), symplectic_form_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        map_bar = 0.0_dp
        symplectic_form_bar = 0.0_dp
        call validate_inputs(map, symplectic_form, defect_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(map_bar) /= shape(map)) .or. &
            any(shape(symplectic_form_bar) /= shape(symplectic_form))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "symplectic-map VJP has incompatible cotangents")
            return
        end if
        map_bar = matmul(symplectic_form, matmul(map, transpose(defect_bar))) + &
            matmul(transpose(symplectic_form), matmul(map, defect_bar))
        symplectic_form_bar = matmul(map, matmul(defect_bar, transpose(map))) - &
            defect_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_symplectic_map_defect_vjp

    subroutine validate_inputs(map, symplectic_form, target, status)
        real(dp), intent(in) :: map(:, :), symplectic_form(:, :), target(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: state_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "symplectic-map defect has incompatible arrays")
        state_count = size(map, 1)
        if (state_count < 1 .or. size(map, 2) /= state_count .or. &
            size(symplectic_form, 1) /= state_count .or. &
            size(symplectic_form, 2) /= state_count .or. &
            size(target, 1) /= state_count .or. size(target, 2) /= state_count) return
        if (.not. all(ieee_is_finite(map)) .or. &
            .not. all(ieee_is_finite(symplectic_form)) .or. &
            .not. all(ieee_is_finite(target))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_inputs

    pure logical function valid_direction( &
            map, symplectic_form, map_dot, symplectic_form_dot) result(valid)
        real(dp), intent(in) :: map(:, :), symplectic_form(:, :)
        real(dp), intent(in) :: map_dot(:, :), symplectic_form_dot(:, :)

        valid = all(shape(map_dot) == shape(map)) .and. &
            all(shape(symplectic_form_dot) == shape(symplectic_form)) .and. &
            all(ieee_is_finite(map_dot)) .and. &
            all(ieee_is_finite(symplectic_form_dot))
    end function valid_direction

end module fortfem_symplectic_map_defect
