module fortfem_meridian_fourier_equilibration
    !! A residual majorant for the retained real Fourier-Galerkin space.
    !! j = K grad u_h - f is projected with the assembly's phi rule. Retain
    !! j_phi and reconstruct R j_pol in RT_k separately in each real mode:
    !! -div_2 sigma_b = Pi_k (R f0 + d_phi j_phi)_b, up to the measured
    !! divergence defect. Normal continuity is inherited from RT_k.
    !!
    !! If K >= chi_lower I everywhere, the residual is bounded by the sum
    !! of the mode-wise flux, oscillation and divergence norms, using
    !! A = R chi_lower and Parseval weights (1, 1/2, 1/2, ...). Thus the
    !! majorant controls the FULL coupled energy, not an n=0 restriction.
    !! The toroidal residual vanishes because j_phi is retained. No inverse
    !! K or uncoupled-coefficient assumption is used. Efficiency may degrade
    !! with anisotropy. Meridional and phi quadrature are NOT enclosed;
    !! the result is practical for the sampled Fourier-Galerkin model.
    !! The finite Fourier band and discrete Dirichlet trace stay fixed.
    !! Storage is O(cells * nq * modes); FFT sampling is O(cells * nq *
    !! nphi log nphi), followed by O(cells * modes) patch reconstruction at
    !! fixed spatial degree. Scratch is allocated once per OpenMP thread.
    use fortfem_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_qgrad
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t, &
        meridian_fourier_source_t
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK, EQUILIBRATION_BAD_INPUT
    use fortfem_generated_fourier_equil_phi_gradient, only: &
        generated_fourier_equil_phi_gradient
    use fortfem_generated_fourier_equil_flux_component, only: &
        generated_fourier_equil_flux_component
    use fortfem_generated_fourier_equil_metric, only: generated_fourier_equil_metric
    use fortfem_generated_fourier_equil_source, only: generated_fourier_equil_source
    use fortnum_fft, only: fortnum_fft_plan_t, fortnum_fft_workspace_t, &
        fft_c2c_plan_init, fft_workspace_init, fft_c2c
    implicit none
    private
    public :: meridian_fourier_equilibrated_estimate, fourier_equilibrated_t

    type :: fourier_equilibrated_t
        real(dp) :: bound = huge(1.0_dp)
        real(dp) :: flux_defect = 0, oscillation = 0, divergence_defect = 0
        ! Unweighted norms per real mode (constant, cos1, sin1, ...).
        real(dp), allocatable :: mode_flux(:), mode_oscillation(:), mode_divergence(:)
        ! Physical RT coefficients of R j_pol; toroidal j_phi at the rule.
        real(dp), allocatable :: flux_coefficients(:, :, :), toroidal_flux(:, :, :)
    end type
contains
    subroutine meridian_fourier_equilibrated_estimate(fs, modes, nphi, coef, &
            source, u, chi_lower, vertices, triangles, xi, eta, weights, &
            coefficient_lower, friedrichs_constant, estimate, status, geometry_nodes)
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        class(meridian_fourier_tensor_t), intent(in) :: coef
        class(meridian_fourier_source_t), intent(in) :: source
        real(dp), intent(in) :: u(:), chi_lower, vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: xi(:), eta(:), weights(:), coefficient_lower(:)
        real(dp), intent(in) :: friedrichs_constant
        type(fourier_equilibrated_t), intent(out) :: estimate
        integer, intent(out) :: status
        real(dp), intent(in), optional :: geometry_nodes(:, :, :)
        type(equilibrated_flux_estimate_t) :: planar
        type(fortnum_fft_plan_t) :: plan
        type(fortnum_fft_workspace_t) :: work
        complex(dp), allocatable :: z(:)
        real(dp), allocatable :: norm(:)
        real(dp), allocatable :: loads(:, :, :), flux(:, :, :, :), scalar(:, :)
        real(dp), allocatable :: xq(:, :), kt(:, :, :), f0(:), fl(:, :)
        real(dp), allocatable :: val(:), grad(:, :), um(:, :), gm(:, :)
        real(dp), allocatable :: samples(:, :), jm(:, :), fm(:)
        real(dp) :: x(2), wq, rr, g(3), sc, ss
        integer :: nb, nt, nl, nq, np, t, q, j, b, n, off, st, failed, c
        real(dp), parameter :: pi = acos(-1.0_dp)

        status = EQUILIBRATION_BAD_INPUT
        if (modes < 0 .or. nphi < 2*modes + 1) return
        if (.not. ieee_is_finite(chi_lower)) return
        if (chi_lower <= 0) return
        nb = 2*modes + 1; nt = fs%mm%topology%n_triangles
        nl = fs%nl0; nq = fs%nq; np = nq*nphi
        if (size(u) /= fs%n0*nb) return
        if (.not. all(ieee_is_finite(u))) return
        if (size(weights) /= nq) return
        if (size(coefficient_lower) /= nt) return
        if (any(coefficient_lower <= 0)) return
        if (size(triangles, 2) /= nt) return
        allocate (norm(nb))
        norm = 0.5_dp; norm(1) = 1
        call fft_c2c_plan_init(plan, nphi)
        allocate (loads(nq, nt, nb), flux(2, nq, nt, nb), scalar(nq, nt))
        allocate (estimate%toroidal_flux(nq, nt, nb))
        failed = 0
        !$omp parallel default(shared) private(t, q, j, b, n, off, st, x, wq, &
        !$omp rr, g, sc, ss, xq, kt, f0, fl, val, grad, um, gm, samples, jm, fm, &
        !$omp work, z, c) &
        !$omp reduction(max:failed)
        allocate (xq(3, np), kt(3, 3, np), f0(np), fl(3, np))
        allocate (val(nl), grad(2, nl), um(nl, nb), gm(3, nb))
        allocate (samples(3, nphi), jm(3, nb), fm(nb), z(nphi))
        call fft_workspace_init(work, plan)
        !$omp do schedule(static)
        do t = 1, nt
            st = 0
            do q = 1, nq
                call meridian_qgrad(fs, t, q, x, wq, val, grad, st)
                if (st /= 0) exit
                do j = 1, nphi
                    xq(:, (q - 1)*nphi + j) = [x, 2*pi*(j - 1)/nphi]
                end do
            end do
            if (st == 0) call coef%tensor(t, xq, kt, st)
            if (st == 0) call source%source(t, xq, f0, fl, st)
            if (st /= 0) then
                failed = 1
                cycle
            end if
            if (.not. all(ieee_is_finite(kt)) .or. &
                .not. all(ieee_is_finite(f0)) .or. &
                .not. all(ieee_is_finite(fl))) then
                failed = 1
                cycle
            end if
            do b = 1, nb
                um(:, b) = u((b - 1)*fs%n0 + fs%sp(0)%dof(:, t))
            end do
            do q = 1, nq
                call meridian_qgrad(fs, t, q, x, wq, val, grad, st)
                rr = x(1)
                if (fs%cyl_r0 > 0) rr = fs%cyl_r0
                gm(1, :) = matmul(grad(1, :), um)
                gm(3, :) = matmul(grad(2, :), um)
                gm(2, 1) = 0
                do n = 1, modes
                    call generated_fourier_equil_phi_gradient(rr, real(n, dp), &
                        dot_product(val, um(:, 2*n)), &
                        dot_product(val, um(:, 2*n + 1)), &
                        gm(2, 2*n), gm(2, 2*n + 1))
                end do
                do c = 1, 3
                    call synthesize(gm(c, :), z, plan, work)
                    samples(c, :) = real(z, dp)
                end do
                do j = 1, nphi
                    off = (q - 1)*nphi + j
                    g = samples(:, j)
                    do c = 1, 3
                        call generated_fourier_equil_flux_component( &
                            kt(c, 1, off), kt(c, 2, off), kt(c, 3, off), &
                            g(1), g(2), g(3), fl(c, off), samples(c, j))
                    end do
                end do
                do c = 1, 3
                    z = cmplx(samples(c, :), 0.0_dp, dp)
                    call analyze(z, jm(c, :), plan, work)
                end do
                z = cmplx(f0((q - 1)*nphi + 1:q*nphi), 0.0_dp, dp)
                call analyze(z, fm, plan, work)
                do b = 1, nb
                    call generated_fourier_equil_metric(rr, chi_lower, &
                        jm(1, b), jm(3, b), fm(b), scalar(q, t), &
                        flux(1, q, t, b), flux(2, q, t, b), loads(q, t, b))
                end do
                estimate%toroidal_flux(q, t, :) = jm(2, :)
                do n = 1, modes
                    call generated_fourier_equil_source(rr, real(n, dp), &
                        fm(2*n), fm(2*n + 1), jm(2, 2*n), jm(2, 2*n + 1), sc, ss)
                    loads(q, t, 2*n) = sc
                    loads(q, t, 2*n + 1) = ss
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (failed /= 0) return
        allocate (estimate%mode_flux(nb), estimate%mode_oscillation(nb), &
            estimate%mode_divergence(nb))
        do b = 1, nb
            call equilibrated_flux_estimate(vertices, triangles, fs%k, xi, eta, &
                weights, scalar, loads(:, :, b), flux(:, :, :, b), &
                coefficient_lower, friedrichs_constant, planar, status, &
                geometry_nodes=geometry_nodes)
            if (status /= EQUILIBRATION_OK) return
            if (b == 1) allocate (estimate%flux_coefficients( &
                size(planar%flux_coefficients, 1), nt, nb))
            estimate%flux_coefficients(:, :, b) = planar%flux_coefficients
            estimate%mode_flux(b) = planar%flux_defect
            estimate%mode_oscillation(b) = planar%oscillation
            estimate%mode_divergence(b) = planar%divergence_defect
        end do
        estimate%flux_defect = sqrt(sum(norm*estimate%mode_flux**2))
        estimate%oscillation = sqrt(sum(norm*estimate%mode_oscillation**2))
        estimate%divergence_defect = sqrt(sum(norm*estimate%mode_divergence**2))
        estimate%bound = estimate%flux_defect + estimate%oscillation &
            + estimate%divergence_defect
        status = EQUILIBRATION_OK
    end subroutine
    subroutine synthesize(amplitude, z, plan, work)
        real(dp), intent(in) :: amplitude(:)
        complex(dp), intent(out) :: z(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer :: n

        z = cmplx(0.0_dp, 0.0_dp, dp)
        z(1) = cmplx(amplitude(1), 0.0_dp, dp)
        do n = 1, (size(amplitude) - 1)/2
            z(n + 1) = cmplx(amplitude(2*n), -amplitude(2*n + 1), dp)
        end do
        call fft_c2c(z, 1, plan, work)
    end subroutine

    subroutine analyze(z, amplitude, plan, work)
        complex(dp), intent(inout) :: z(:)
        real(dp), intent(out) :: amplitude(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        type(fortnum_fft_workspace_t), intent(inout) :: work
        integer :: n

        call fft_c2c(z, -1, plan, work)
        amplitude(1) = real(z(1), dp)/size(z)
        do n = 1, (size(amplitude) - 1)/2
            amplitude(2*n) = 2*real(z(n + 1), dp)/size(z)
            amplitude(2*n + 1) = -2*aimag(z(n + 1))/size(z)
        end do
    end subroutine
end module fortfem_meridian_fourier_equilibration
