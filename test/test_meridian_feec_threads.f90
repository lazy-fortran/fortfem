program test_meridian_feec_threads
    !! Thread invariance of the element-parallel assembly: operators and
    !! masses agree between thread counts to roundoff (they are in fact
    !! bitwise equal, the triplet slices being fixed per element); the
    !! reductions of loads and error norms agree to a relative 1e-12.
    !! Vacuous in a build without OpenMP.
    !$ use omp_lib, only: omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_load, meridian_error_norm, meridian_gs_stiffness
    use fortsparse, only: csc_t
    use meridian_test_fields, only: mode_n, vec_a, fun_astar
    use meridian_test_meshes, only: make_annulus_mesh
    implicit none

    integer, parameter :: counts(4) = [1, 2, 4, 8]
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: ref, fs
    type(csc_t) :: gs_ref, gs
    complex(dp), allocatable :: b_ref(:), b(:), c(:)
    real(dp) :: e_ref, e, dmax, scale
    logical :: all_passed
    integer :: i, status, sp

    all_passed = .true.
    mode_n = 3
    call make_annulus_mesh(16, 64, 0.4_dp, 1.0_dp, 3.0_dp, 2, mm)
    !$ call omp_set_num_threads(1)
    call meridian_feec_build(mm, 2, ref, status)
    call meridian_gs_stiffness(ref, gs_ref, status)
    call meridian_load(ref, 1, vec_a, b_ref, status)
    allocate (c(size(b_ref)))
    c = b_ref
    call meridian_error_norm(ref, 1, c, fun_astar, e_ref, status)
    do i = 2, size(counts)
        !$ call omp_set_num_threads(counts(i))
        call meridian_feec_build(mm, 2, fs, status)
        dmax = 0.0_dp
        scale = 0.0_dp
        do sp = 0, 3
            call diff(ref%mass(sp), fs%mass(sp), dmax, scale)
        end do
        call diff(ref%grad%a, fs%grad%a, dmax, scale)
        call diff(ref%grad%b, fs%grad%b, dmax, scale)
        call diff(ref%curl%a, fs%curl%a, dmax, scale)
        call diff(ref%curl%b, fs%curl%b, dmax, scale)
        call diff(ref%div%a, fs%div%a, dmax, scale)
        call diff(ref%div%b, fs%div%b, dmax, scale)
        call meridian_gs_stiffness(fs, gs, status)
        call diff(gs_ref, gs, dmax, scale)
        call record(dmax <= 1.0e-13_dp*scale, "operators and masses invariant")
        call meridian_load(fs, 1, vec_a, b, status)
        call meridian_error_norm(fs, 1, c, fun_astar, e, status)
        call record(maxval(abs(b - b_ref)) < 1.0e-12_dp*maxval(abs(b_ref)), &
            "load vector invariant")
        call record(abs(e - e_ref) < 1.0e-12_dp*e_ref, "error norm invariant")
        print '(a,i0,a,es10.2)', "  threads ", counts(i), &
            " max matrix difference ", dmax
    end do
    call check_summary("Meridian per-n FEEC: thread invariance")
    if (.not. all_passed) error stop 1

contains

    subroutine record(cond, d)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. cond
        call check_condition(cond, d)
    end subroutine record

    subroutine diff(a, b, dmax, scale)
        type(csc_t), intent(in) :: a, b
        real(dp), intent(inout) :: dmax, scale
        logical :: same

        same = a%nnz == b%nnz
        if (same) same = all(a%col_ptr == b%col_ptr)
        if (same) same = all(a%row_idx == b%row_idx)
        if (.not. same) then
            dmax = huge(1.0_dp)
            return
        end if
        if (a%nnz > 0) then
            dmax = max(dmax, maxval(abs(a%val - b%val)))
            scale = max(scale, maxval(abs(a%val)))
        end if
    end subroutine diff

end program test_meridian_feec_threads
