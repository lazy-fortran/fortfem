module fortfem_generated_lab_guenter
    use fortfem_generated_lab_guenter_kernel, only: lab_guenter_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: lab_guenter
contains
    pure subroutine lab_guenter(x, y, values)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: values(8)
        call lab_guenter_kernel(x, y, acos(-1.0_dp), values)
    end subroutine
end module
