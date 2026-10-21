module fortfem_codegen_triangle_basis
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*)
    implicit none
    private
    public :: triangle_scalar_bases
contains
    subroutine triangle_scalar_bases(x, y, p1, p2)
        type(expr_t), intent(in) :: x, y
        type(expr_t), intent(out) :: p1(3), p2(6)
        integer :: i
        p1(1) = 1 - x - y; p1(2) = x; p1(3) = y
        do i = 1, 3
            p2(i) = p1(i)*(2*p1(i) - 1)
        end do
        p2(4) = 4*p1(1)*p1(2)
        p2(5) = 4*p1(2)*p1(3)
        p2(6) = 4*p1(3)*p1(1)
    end subroutine triangle_scalar_bases
end module fortfem_codegen_triangle_basis
