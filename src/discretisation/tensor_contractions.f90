module fortfem_tensor_contractions
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: tensor_contract
contains
    pure subroutine tensor_contract(ni, no, a, b, c, u, v, t1, t2)
        integer, intent(in) :: ni(3), no(3)
        real(dp), intent(in) :: a(:, :), b(:, :), c(:, :), u(:)
        real(dp), intent(out) :: v(:)
        real(dp), intent(inout) :: t1(:), t2(:)
        integer :: i, j, k, r, base, dest
        ! Three one-dimensional contractions, with the leftmost index contiguous.
        t1(1:no(1)*ni(2)*ni(3)) = 0.0_dp
        do k = 1, ni(3)
            do j = 1, ni(2)
                base = ni(1)*(j - 1 + ni(2)*(k - 1))
                dest = no(1)*(j - 1 + ni(2)*(k - 1))
                do r = 1, ni(1)
                    !$omp simd
                    do i = 1, no(1)
                        t1(dest + i) = t1(dest + i) + a(i, r)*u(base + r)
                    end do
                end do
            end do
        end do
        t2(1:no(1)*no(2)*ni(3)) = 0.0_dp
        do k = 1, ni(3)
            do j = 1, no(2)
                do r = 1, ni(2)
                    !$omp simd
                    do i = 1, no(1)
                        dest = i + no(1)*(j - 1 + no(2)*(k - 1))
                        base = i + no(1)*(r - 1 + ni(2)*(k - 1))
                        t2(dest) = t2(dest) + b(j, r)*t1(base)
                    end do
                end do
            end do
        end do
        v(1:product(no)) = 0.0_dp
        do k = 1, no(3)
            do r = 1, ni(3)
                !$omp simd
                do i = 1, no(1)*no(2)
                    dest = i + no(1)*no(2)*(k - 1)
                    base = i + no(1)*no(2)*(r - 1)
                    v(dest) = v(dest) + c(k, r)*t2(base)
                end do
            end do
        end do
    end subroutine
end module
