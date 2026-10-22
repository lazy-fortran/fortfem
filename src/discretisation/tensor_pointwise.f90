module fortfem_tensor_pointwise
    ! Quadrature-point field access for nonlinear matrix-free operators on the
    ! lab tensor meshes: interpolation of a nodal field to all quadrature
    ! points, a mass action with a pointwise weight and a load vector from
    ! pointwise values. All three are sum-factorised and OpenMP-parallel; the
    ! scattering routines reuse the operator's cell colouring.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_operator, only: tensor_operator_t
    use fortfem_tensor_contractions, only: tensor_contract
    implicit none
    private
    public :: tensor_quadrature_values, tensor_weighted_mass, tensor_weighted_load

contains

    subroutine tensor_quadrature_values(op, x, values)
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:, :) ! (quadrature points, cells)
        integer :: e, nl, mq, ws
        nl = product(op%ni); mq = product(op%no)
        ws = maxval([op%ni, op%no])**3
        !$omp parallel default(shared) private(e)
        block
            real(dp), allocatable :: u(:), t1(:), t2(:), b(:, :)
            real(dp) :: identity(1, 1)
            allocate (u(nl), t1(ws), t2(ws))
            b = op%line%b(:, :, 1)
            identity = 1
            !$omp do schedule(static)
            do e = 1, op%mesh%nelem
                u = x(op%mesh%dof(:, e))
                if (op%mesh%dim == 2) then
                    call tensor_contract(op%ni, op%no, b, b, identity, u, &
                                         values(:, e), t1, t2)
                else
                    call tensor_contract(op%ni, op%no, b, b, b, u, values(:, e), &
                                         t1, t2)
                end if
            end do
            !$omp end do
        end block
        !$omp end parallel
    end subroutine

    subroutine tensor_weighted_mass(op, weight, x, y, eliminate)
        ! y = sum over cells of B^T diag(weight * volume) B x. With eliminate,
        ! boundary inputs are masked and boundary rows of y are zero.
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: weight(:, :), x(:)
        real(dp), intent(out) :: y(:)
        logical, intent(in), optional :: eliminate
        logical :: mask
        integer :: e, ce, colour, a, nl, mq, ws
        mask = .true.
        if (present(eliminate)) mask = eliminate
        nl = product(op%ni); mq = product(op%no)
        ws = maxval([op%ni, op%no])**3
        y = 0
        !$omp parallel default(shared) private(e,ce,colour,a)
        block
            real(dp), allocatable :: u(:), q(:), v(:), t1(:), t2(:), b(:, :), &
                                     bt(:, :)
            real(dp) :: identity(1, 1)
            allocate (u(nl), q(mq), v(nl), t1(ws), t2(ws))
            b = op%line%b(:, :, 1)
            bt = transpose(b)
            identity = 1
            do colour = 1, op%ncolours
                !$omp do schedule(static)
                do ce = op%colour_ptr(colour), op%colour_ptr(colour + 1) - 1
                    e = op%colour_cells(ce)
                    do a = 1, nl
                        u(a) = x(op%mesh%dof(a, e))
                        if (mask) then
                            if (op%mesh%boundary(op%mesh%dof(a, e))) u(a) = 0
                        end if
                    end do
                    if (op%mesh%dim == 2) then
                        call tensor_contract(op%ni, op%no, b, b, identity, u, q, t1, t2)
                        q = q*weight(:, e)*op%volume(:, e)
                        call tensor_contract(op%no, op%ni, bt, bt, identity, q, v, &
                                             t1, t2)
                    else
                        call tensor_contract(op%ni, op%no, b, b, b, u, q, t1, t2)
                        q = q*weight(:, e)*op%volume(:, e)
                        call tensor_contract(op%no, op%ni, bt, bt, bt, q, v, t1, t2)
                    end if
                    do a = 1, nl
                        y(op%mesh%dof(a, e)) = y(op%mesh%dof(a, e)) + v(a)
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        if (mask) then
            where (op%mesh%boundary) y = 0
        end if
    end subroutine

    subroutine tensor_weighted_load(op, values, rhs, eliminate)
        ! rhs = sum over cells of B^T (values * volume); values at quadrature
        ! points. With eliminate (default), boundary rows are zero.
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: values(:, :)
        real(dp), intent(out) :: rhs(:)
        logical, intent(in), optional :: eliminate
        logical :: mask
        integer :: e, ce, colour, a, nl, mq, ws
        mask = .true.
        if (present(eliminate)) mask = eliminate
        nl = product(op%ni); mq = product(op%no)
        ws = maxval([op%ni, op%no])**3
        rhs = 0
        !$omp parallel default(shared) private(e,ce,colour,a)
        block
            real(dp), allocatable :: q(:), v(:), t1(:), t2(:), bt(:, :)
            real(dp) :: identity(1, 1)
            allocate (q(mq), v(nl), t1(ws), t2(ws))
            bt = transpose(op%line%b(:, :, 1))
            identity = 1
            do colour = 1, op%ncolours
                !$omp do schedule(static)
                do ce = op%colour_ptr(colour), op%colour_ptr(colour + 1) - 1
                    e = op%colour_cells(ce)
                    q = values(:, e)*op%volume(:, e)
                    if (op%mesh%dim == 2) then
                        call tensor_contract(op%no, op%ni, bt, bt, identity, q, v, &
                                             t1, t2)
                    else
                        call tensor_contract(op%no, op%ni, bt, bt, bt, q, v, t1, t2)
                    end if
                    do a = 1, nl
                        rhs(op%mesh%dof(a, e)) = rhs(op%mesh%dof(a, e)) + v(a)
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        if (mask) then
            where (op%mesh%boundary) rhs = 0
        end if
    end subroutine
end module fortfem_tensor_pointwise
