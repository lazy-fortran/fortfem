module fortfem_tensor_line_tables
    ! Lab contract (fortfem_lab controller, 2026-10-22). Copy verbatim to
    ! src/discretisation/tensor_line_tables.f90. Per-direction 1D tables that every
    ! tensor-product (sum-factorised) operator consumes. One instance per direction;
    ! GLL Lagrange, B-spline/IGA and Fourier builders all fill this type.
    ! Extensions go in separate modules; change this type only via the controller.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: LINE_GLL = 1
    integer, parameter, public :: LINE_BSPLINE = 2
    integer, parameter, public :: LINE_FOURIER = 3

    public :: tensor_line_tables_t

    type :: tensor_line_tables_t
        integer :: family = 0 ! LINE_GLL, LINE_BSPLINE or LINE_FOURIER
        integer :: degree = 0 ! polynomial degree; Fourier: max mode number
        integer :: nelem = 0 ! 1D cells (Fourier: 1)
        integer :: nloc = 0 ! local functions per cell
        integer :: nq = 0 ! quadrature points per cell
        integer :: ndof = 0 ! global 1D dofs after periodic identification
        logical :: periodic = .false.
        real(dp) :: x0 = 0.0_dp ! parameter interval start
        real(dp) :: x1 = 1.0_dp ! parameter interval end
        real(dp), allocatable :: xq(:, :) ! (nq, nelem) quadrature points in [x0, x1]
        real(dp), allocatable :: wq(:, :) ! (nq, nelem) weights incl. cell length
        real(dp), allocatable :: b(:, :, :) ! (nq, nloc, nelem) basis values
        real(dp), allocatable :: d(:, :, :) ! (nq, nloc, nelem) d/dx in parameter
        integer, allocatable :: dof(:, :) ! (nloc, nelem) global 1D index, 1-based
    end type tensor_line_tables_t
end module fortfem_tensor_line_tables
