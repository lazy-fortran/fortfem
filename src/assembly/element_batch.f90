module fortfem_element_batch
    !! SoA element blocks: contiguous elements, then local row, then column.
    !! Storage is selected by input kind; output kind selects arithmetic.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: p1_element_batch
    interface p1_element_batch
        module procedure p1_batch_fp64, p1_batch_fp32, p1_batch_mixed
    end interface
contains
    subroutine p1_batch_fp64(jac, diffusion, reaction, blocks, status)
        use fortfem_generated_p1_cell_fp64, only: cell => generated_p1_cell_fp64
        use fortfem_generated_triangle_measure_fp64, only: &
            measure => generated_triangle_measure_fp64
        integer, parameter :: ck = real64
        real(real64), contiguous, intent(in) :: jac(:, :), diffusion(:), reaction(:)
        real(real64), contiguous, intent(out) :: blocks(:, :, :)
        integer, contiguous, intent(out) :: status(:)
        include 'fortfem_p1_element_batch.inc'
    end subroutine p1_batch_fp64

    subroutine p1_batch_fp32(jac, diffusion, reaction, blocks, status)
        use fortfem_generated_p1_cell_fp32, only: cell => generated_p1_cell_fp32
        use fortfem_generated_triangle_measure_fp32, only: &
            measure => generated_triangle_measure_fp32
        integer, parameter :: ck = real32
        real(real32), contiguous, intent(in) :: jac(:, :), diffusion(:), reaction(:)
        real(real32), contiguous, intent(out) :: blocks(:, :, :)
        integer, contiguous, intent(out) :: status(:)
        include 'fortfem_p1_element_batch.inc'
    end subroutine p1_batch_fp32

    subroutine p1_batch_mixed(jac, diffusion, reaction, blocks, status)
        use fortfem_generated_p1_cell_fp64, only: cell => generated_p1_cell_fp64
        use fortfem_generated_triangle_measure_fp64, only: &
            measure => generated_triangle_measure_fp64
        integer, parameter :: ck = real64
        real(real32), contiguous, intent(in) :: jac(:, :), diffusion(:), reaction(:)
        real(real64), contiguous, intent(out) :: blocks(:, :, :)
        integer, contiguous, intent(out) :: status(:)
        include 'fortfem_p1_element_batch.inc'
    end subroutine p1_batch_mixed
end module fortfem_element_batch
