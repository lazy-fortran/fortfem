module fortfem_kinds
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter, public :: dp = real64

    real(dp), parameter, public :: pi = 3.14159265358979323846_dp

end module fortfem_kinds
