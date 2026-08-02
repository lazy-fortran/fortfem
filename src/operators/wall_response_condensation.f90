module fortfem_wall_response_condensation
    !! Neutral Schur reduction for field/wall response blocks.
    !!
    !! Given blocks M_ee, M_ey, M_ye, and M_yy, this module returns
    !! M_ee - M_ey M_yy^{-1} M_ye.  It owns no wall basis or normalization;
    !! the operation is reusable for FEM/BEM, DtN, PML, and external response
    !! adapters.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    implicit none
    private

    public :: condense_wall_response_blocks
    public :: condense_wall_response_blocks_jvp
    public :: condense_wall_response_blocks_vjp

contains

    subroutine condense_wall_response_blocks( &
            exterior_exterior, exterior_wall, wall_exterior, wall_wall, &
            effective, status)
        complex(dp), intent(in) :: exterior_exterior(:, :), exterior_wall(:, :)
        complex(dp), intent(in) :: wall_exterior(:, :), wall_wall(:, :)
        complex(dp), intent(out) :: effective(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: wall_to_exterior(:, :)
        integer :: info

        effective = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_blocks( &
            exterior_exterior, exterior_wall, wall_exterior, wall_wall, effective)) return
        allocate(wall_to_exterior(size(wall_wall, 1), size(exterior_exterior, 2)))
        call dense_solve(wall_wall, wall_exterior, wall_to_exterior, info)
        if (info /= 0) then
            status = 2
            return
        end if
        effective = exterior_exterior - matmul(exterior_wall, wall_to_exterior)
        status = 0
    end subroutine condense_wall_response_blocks

    subroutine condense_wall_response_blocks_jvp( &
            exterior_exterior, exterior_wall, wall_exterior, wall_wall, &
            exterior_exterior_dot, exterior_wall_dot, wall_exterior_dot, &
            wall_wall_dot, effective_dot, status)
        complex(dp), intent(in) :: exterior_exterior(:, :), exterior_wall(:, :)
        complex(dp), intent(in) :: wall_exterior(:, :), wall_wall(:, :)
        complex(dp), intent(in) :: exterior_exterior_dot(:, :), exterior_wall_dot(:, :)
        complex(dp), intent(in) :: wall_exterior_dot(:, :), wall_wall_dot(:, :)
        complex(dp), intent(out) :: effective_dot(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: wall_to_exterior(:, :), wall_to_exterior_dot(:, :)
        complex(dp), allocatable :: rhs_dot(:, :)
        integer :: info

        effective_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_blocks( &
            exterior_exterior, exterior_wall, wall_exterior, wall_wall, effective_dot)) return
        if (any(shape(exterior_exterior_dot) /= shape(exterior_exterior)) .or. &
            any(shape(exterior_wall_dot) /= shape(exterior_wall)) .or. &
            any(shape(wall_exterior_dot) /= shape(wall_exterior)) .or. &
            any(shape(wall_wall_dot) /= shape(wall_wall))) return
        allocate( &
            wall_to_exterior(size(wall_wall, 1), size(exterior_exterior, 2)), &
            wall_to_exterior_dot(size(wall_wall, 1), size(exterior_exterior, 2)), &
            rhs_dot(size(wall_wall, 1), size(exterior_exterior, 2)))
        call dense_solve(wall_wall, wall_exterior, wall_to_exterior, info)
        if (info /= 0) then
            status = 2
            return
        end if
        rhs_dot = wall_exterior_dot - matmul(wall_wall_dot, wall_to_exterior)
        call dense_solve(wall_wall, rhs_dot, wall_to_exterior_dot, info)
        if (info /= 0) then
            status = 2
            return
        end if
        effective_dot = exterior_exterior_dot - &
            matmul(exterior_wall_dot, wall_to_exterior) - &
            matmul(exterior_wall, wall_to_exterior_dot)
        status = 0
    end subroutine condense_wall_response_blocks_jvp

    subroutine condense_wall_response_blocks_vjp( &
            exterior_exterior, exterior_wall, wall_exterior, wall_wall, effective_bar, &
            exterior_exterior_bar, exterior_wall_bar, wall_exterior_bar, &
            wall_wall_bar, status)
        complex(dp), intent(in) :: exterior_exterior(:, :), exterior_wall(:, :)
        complex(dp), intent(in) :: wall_exterior(:, :), wall_wall(:, :)
        complex(dp), intent(in) :: effective_bar(:, :)
        complex(dp), allocatable, intent(out) :: exterior_exterior_bar(:, :)
        complex(dp), allocatable, intent(out) :: exterior_wall_bar(:, :)
        complex(dp), allocatable, intent(out) :: wall_exterior_bar(:, :)
        complex(dp), allocatable, intent(out) :: wall_wall_bar(:, :)
        complex(dp), allocatable :: wall_to_exterior(:, :), wall_to_exterior_bar(:, :)
        complex(dp), allocatable :: adjoint(:, :)
        integer, intent(out) :: status
        integer :: info

        if (allocated(exterior_exterior_bar)) deallocate(exterior_exterior_bar)
        if (allocated(exterior_wall_bar)) deallocate(exterior_wall_bar)
        if (allocated(wall_exterior_bar)) deallocate(wall_exterior_bar)
        if (allocated(wall_wall_bar)) deallocate(wall_wall_bar)
        allocate( &
            exterior_exterior_bar(size(exterior_exterior, 1), size(exterior_exterior, 2)), &
            exterior_wall_bar(size(exterior_wall, 1), size(exterior_wall, 2)), &
            wall_exterior_bar(size(wall_exterior, 1), size(wall_exterior, 2)), &
            wall_wall_bar(size(wall_wall, 1), size(wall_wall, 2)))
        exterior_exterior_bar = cmplx(0.0_dp, 0.0_dp, dp)
        exterior_wall_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wall_exterior_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wall_wall_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_blocks( &
            exterior_exterior, exterior_wall, wall_exterior, wall_wall, effective_bar)) return
        allocate(wall_to_exterior(size(wall_wall, 1), size(exterior_exterior, 2)))
        call dense_solve(wall_wall, wall_exterior, wall_to_exterior, info)
        if (info /= 0) then
            status = 2
            return
        end if
        exterior_exterior_bar = effective_bar
        exterior_wall_bar = -matmul( &
            effective_bar, conjg(transpose(wall_to_exterior)))
        allocate(wall_to_exterior_bar(size(wall_wall, 1), size(exterior_exterior, 2)))
        wall_to_exterior_bar = -matmul( &
            conjg(transpose(exterior_wall)), effective_bar)
        allocate(adjoint(size(wall_wall, 1), size(exterior_exterior, 2)))
        call dense_solve(conjg(transpose(wall_wall)), wall_to_exterior_bar, adjoint, info)
        if (info /= 0) then
            status = 2
            return
        end if
        wall_exterior_bar = adjoint
        wall_wall_bar = -matmul( &
            adjoint, conjg(transpose(wall_to_exterior)))
        status = 0
    end subroutine condense_wall_response_blocks_vjp

    logical function valid_blocks( &
            exterior_exterior, exterior_wall, wall_exterior, wall_wall, output) &
            result(valid)
        complex(dp), intent(in) :: exterior_exterior(:, :), exterior_wall(:, :)
        complex(dp), intent(in) :: wall_exterior(:, :), wall_wall(:, :), output(:, :)
        integer :: exterior_count, wall_count

        exterior_count = size(exterior_exterior, 1)
        wall_count = size(wall_wall, 1)
        valid = exterior_count > 0 .and. wall_count > 0 .and. &
            size(exterior_exterior, 2) == exterior_count .and. &
            size(wall_wall, 2) == wall_count .and. &
            size(exterior_wall, 1) == exterior_count .and. &
            size(exterior_wall, 2) == wall_count .and. &
            size(wall_exterior, 1) == wall_count .and. &
            size(wall_exterior, 2) == exterior_count .and. &
            size(output, 1) == exterior_count .and. &
            size(output, 2) == exterior_count
    end function valid_blocks

end module fortfem_wall_response_condensation
