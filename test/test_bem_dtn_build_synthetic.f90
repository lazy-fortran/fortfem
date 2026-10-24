program test_bem_dtn_build_synthetic
    !! build_exterior_dtn on synthetic rectangular operators (nf /= nv):
    !! flux map solves V q = -(M/2-K) u and the Steklov matrix equals the
    !! explicit-loop Schur complement W - J^T Q, symmetrised.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: build_exterior_dtn, &
        bem_operators_t, bem_exterior_dtn_t
    implicit none
    integer, parameter :: nf = 37, nv = 23
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    real(dp), allocatable :: flux_map(:, :), steklov(:, :), reference(:, :)
    real(dp), allocatable :: jump(:, :), unit(:), flux(:), load(:)
    real(dp) :: scale
    integer :: i, j, p, status
    logical :: passed

    passed = .true.
    allocate (op%v(nf, nf), op%k(nf, nv), op%mass(nf, nv), op%w(nv, nv))
    do j = 1, nf
        do i = 1, nf
            op%v(i, j) = 1.0_dp/(1.0_dp + real((i - j)**2, dp))
        end do
        op%v(j, j) = op%v(j, j) + 2.0_dp
    end do
    do j = 1, nv
        do i = 1, nf
            op%k(i, j) = sin(0.3_dp*real(i*j, dp))/real(i + j, dp)
            op%mass(i, j) = merge(1.0_dp, 0.0_dp, i == j)
        end do
    end do
    do j = 1, nv
        do i = 1, nv
            op%w(i, j) = 1.0_dp/(1.0_dp + real(abs(i - j), dp))
        end do
    end do

    call build_exterior_dtn(op, exterior, status)
    call record(status == 0, 'rectangular synthetic DtN builds')
    if (status /= 0) error stop 1
    call exterior%matrix(steklov, status)

    allocate (flux_map(nf, nv), unit(nv), flux(nf), load(nv))
    do j = 1, nv
        unit = 0.0_dp
        unit(j) = 1.0_dp
        call exterior%apply(unit, flux, load, status)
        flux_map(:, j) = flux
    end do
    jump = 0.5_dp*op%mass - op%k
    call record(maxval(abs(matmul(op%v, flux_map) + jump)) < 1.e-12_dp, &
        'flux map solves V q = -(M/2-K) u')

    allocate (reference(nv, nv))
    do j = 1, nv
        do i = 1, nv
            reference(i, j) = op%w(i, j)
            do p = 1, nf
                reference(i, j) = reference(i, j) - jump(p, i)*flux_map(p, j)
            end do
        end do
    end do
    reference = 0.5_dp*(reference + transpose(reference))
    scale = maxval(abs(reference))
    call record(maxval(abs(steklov - reference)) < 1.e-13_dp*scale, &
        'Steklov matrix equals the explicit Schur complement')
    call record(maxval(abs(steklov - transpose(steklov))) == 0.0_dp, &
        'Steklov matrix is exactly symmetric')
    call check_summary('Synthetic rectangular DtN build')
    if (.not. passed) error stop 1
contains
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
