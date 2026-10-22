module fortfem_bem_circulant
    !! A(r,s)=B(modulo(s-r,nphi)), with caller-visible global DOF maps.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_fft, only: fft_c2c, fft_c2c_plan_init, fortnum_fft_plan_t, &
                           fortnum_fft_workspace_t, fft_workspace_init
    use fortnum_linalg, only: dense_solve
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t
    implicit none
    private
    public :: circulant_matrix_t, bem_circulant_t, assignment(=)

    type :: circulant_matrix_t
        real(dp), allocatable :: blocks(:, :, :)
        integer, allocatable :: row_dofs(:, :), column_dofs(:, :)
    contains
        procedure :: valid => valid_matrix
        procedure :: expand => expand_matrix
        procedure :: apply => apply_matrix
        procedure :: harmonics => matrix_harmonics
        procedure :: solve => solve_matrix
        procedure :: symmetrize => symmetrize_matrix
    end type

    type :: bem_circulant_t
        type(circulant_matrix_t) :: v, k, w, mass
        real(dp), allocatable :: area(:)
        integer :: genus = -1
    contains
        procedure :: expand => expand_bem
    end type

    interface assignment(=)
        module procedure assign_matrix, assign_bem
    end interface
contains
    subroutine assign_matrix(lhs, rhs)
        type(circulant_matrix_t), intent(out) :: lhs
        type(circulant_matrix_t), intent(in) :: rhs
        if (allocated(rhs%blocks)) lhs%blocks = rhs%blocks
        if (allocated(rhs%row_dofs)) lhs%row_dofs = rhs%row_dofs
        if (allocated(rhs%column_dofs)) lhs%column_dofs = rhs%column_dofs
    end subroutine

    subroutine assign_bem(lhs, rhs)
        type(bem_circulant_t), intent(out) :: lhs
        type(bem_circulant_t), intent(in) :: rhs
        lhs%v = rhs%v
        lhs%k = rhs%k
        lhs%w = rhs%w
        lhs%mass = rhs%mass
        if (allocated(rhs%area)) lhs%area = rhs%area
        lhs%genus = rhs%genus
    end subroutine

    logical function valid_map(dofs) result(valid)
        integer, intent(in) :: dofs(:, :)
        integer :: i, j
        integer :: counts(size(dofs))
        valid = .false.
        if (any(dofs < 1)) return
        if (any(dofs > size(dofs))) return
        counts = 0
        do j = 1, size(dofs, 2)
            do i = 1, size(dofs, 1)
                counts(dofs(i, j)) = counts(dofs(i, j)) + 1
            end do
        end do
        valid = all(counts == 1)
    end function

    logical function valid_matrix(self) result(valid)
        class(circulant_matrix_t), intent(in) :: self
        integer :: nr, nc, ns
        valid = .false.
        if (.not. allocated(self%blocks)) return
        if (.not. allocated(self%row_dofs)) return
        if (.not. allocated(self%column_dofs)) return
        nr = size(self%blocks, 1); nc = size(self%blocks, 2)
        ns = size(self%blocks, 3)
        if (min(nr, nc, ns) < 1) return
        if (size(self%row_dofs, 1) /= nr) return
        if (size(self%row_dofs, 2) /= ns) return
        if (size(self%column_dofs, 1) /= nc) return
        if (size(self%column_dofs, 2) /= ns) return
        if (.not. valid_map(self%row_dofs)) return
        if (.not. valid_map(self%column_dofs)) return
        if (.not. all(ieee_is_finite(self%blocks))) return
        valid = .true.
    end function

    subroutine expand_matrix(self, matrix, status)
        class(circulant_matrix_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        integer :: r, s, ns
        status = 1
        if (.not. self%valid()) return
        ns = size(self%blocks, 3)
        allocate (matrix(size(self%row_dofs), size(self%column_dofs)))
        !$omp parallel do collapse(2) default(shared) private(r,s)
        do s = 1, ns
            do r = 1, ns
                matrix(self%row_dofs(:, r), self%column_dofs(:, s)) = &
                    self%blocks(:, :, modulo(s - r, ns) + 1)
            end do
        end do
        !$omp end parallel do
        status = 0
    end subroutine

    subroutine apply_matrix(self, x, y, status)
        class(circulant_matrix_t), intent(in) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        real(dp), allocatable :: sector_x(:, :), sector_y(:, :)
        integer :: r, s, ns, i, j
        status = 1
        if (.not. self%valid()) return
        if (size(x) /= size(self%column_dofs)) return
        if (size(y) /= size(self%row_dofs)) return
        if (.not. all(ieee_is_finite(x))) return
        ns = size(self%blocks, 3)
        allocate (sector_x(size(self%blocks, 2), ns), &
                  sector_y(size(self%blocks, 1), ns))
        do s = 1, ns
            do i = 1, size(sector_x, 1)
                sector_x(i, s) = x(self%column_dofs(i, s))
            end do
        end do
        sector_y = 0._dp
        !$omp parallel do default(shared) private(r,s,i,j)
        do r = 1, ns
            do s = 1, ns
                do j = 1, size(sector_x, 1)
                    do i = 1, size(sector_y, 1)
                        sector_y(i, r) = sector_y(i, r) + &
                                 self%blocks(i, j, modulo(s - r, ns) + 1)*sector_x(j, s)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        do r = 1, ns
            do i = 1, size(sector_y, 1)
                y(self%row_dofs(i, r)) = sector_y(i, r)
            end do
        end do
        status = 0
    end subroutine

    subroutine matrix_harmonics(self, modes, status)
        !! Forward-transforming sector vectors (negative exponent) gives
        !! A_hat(k)=sum_s B(s)*exp(+2*pi*i*k*s/nphi), without normalization.
        class(circulant_matrix_t), intent(in) :: self
        complex(dp), allocatable, intent(out) :: modes(:, :, :)
        integer, intent(out) :: status
        type(fortnum_fft_plan_t) :: plan
        type(fortnum_fft_workspace_t) :: workspace
        integer :: i, j
        status = 1
        if (.not. self%valid()) return
        call fft_c2c_plan_init(plan, size(self%blocks, 3))
        modes = cmplx(self%blocks, 0._dp, dp)
        !$omp parallel default(shared) private(i,j,workspace)
        call fft_workspace_init(workspace, plan)
        !$omp do collapse(2)
        do j = 1, size(modes, 2)
            do i = 1, size(modes, 1)
                call fft_c2c(modes(i, j, :), 1, plan, workspace)
            end do
        end do
        !$omp end do
        !$omp end parallel
        status = 0
    end subroutine

    subroutine solve_matrix(self, rhs, x, status)
        !! Independent complex dense solves, one per Fourier harmonic.
        class(circulant_matrix_t), intent(in) :: self
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: modes(:, :, :), b(:, :), solution(:, :)
        type(fortnum_fft_plan_t) :: plan
        integer, allocatable :: info(:)
        integer :: i, s, n, ns
        status = 1
        if (.not. self%valid()) return
        n = size(self%blocks, 1); ns = size(self%blocks, 3)
        if (size(self%blocks, 2) /= n) return
        if (size(rhs) /= n*ns .or. size(x) /= n*ns) return
        if (.not. all(ieee_is_finite(rhs))) return
        call self%harmonics(modes, status)
        if (status /= 0) return
        call fft_c2c_plan_init(plan, ns)
        allocate (b(n, ns), solution(n, ns), info(ns))
        do s = 1, ns
            b(:, s) = cmplx(rhs(self%row_dofs(:, s)), 0._dp, dp)
        end do
        !$omp parallel do default(shared) private(i)
        do i = 1, n
            call fft_c2c(b(i, :), -1, plan)
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(s)
        do s = 1, ns
            call dense_solve(modes(:, :, s), b(:, s), solution(:, s), info(s))
        end do
        !$omp end parallel do
        status = maxval(abs(info))
        if (status /= 0) return
        !$omp parallel do default(shared) private(i)
        do i = 1, n
            call fft_c2c(solution(i, :), 1, plan)
        end do
        !$omp end parallel do
        do s = 1, ns
            x(self%column_dofs(:, s)) = real(solution(:, s), dp)/real(ns, dp)
        end do
    end subroutine

    subroutine symmetrize_matrix(self, raw_asymmetry)
        class(circulant_matrix_t), intent(inout) :: self
        real(dp), intent(out) :: raw_asymmetry
        real(dp), allocatable :: original(:, :, :)
        integer :: s, ns
        ns = size(self%blocks, 3)
        original = self%blocks
        raw_asymmetry = 0._dp
        do s = 1, ns
            raw_asymmetry = max(raw_asymmetry, maxval(abs(original(:, :, s) - &
                                     transpose(original(:, :, modulo(1 - s, ns) + 1)))))
            self%blocks(:, :, s) = .5_dp*(original(:, :, s) + &
                                       transpose(original(:, :, modulo(1 - s, ns) + 1)))
        end do
        raw_asymmetry = raw_asymmetry/ &
                        max(tiny(1._dp), maxval(abs(original)))
    end subroutine

    subroutine expand_bem(self, op, status)
        class(bem_circulant_t), intent(in) :: self
        type(bem_operators_t), intent(out) :: op
        integer, intent(out) :: status
        integer :: i, j
        status = 1
        if (.not. allocated(self%area)) return
        call self%v%expand(op%v, status)
        if (status /= 0) return
        call self%k%expand(op%k, status)
        if (status /= 0) return
        call self%w%expand(op%w, status)
        if (status /= 0) return
        call self%mass%expand(op%mass, status)
        if (status /= 0) return
        allocate (op%adjoint(size(op%k, 2), size(op%k, 1)))
        do j = 1, size(op%k, 1)
            do i = 1, size(op%k, 2)
                op%adjoint(i, j) = op%k(j, i)
            end do
        end do
        op%area = self%area
        op%genus = self%genus
    end subroutine
end module
