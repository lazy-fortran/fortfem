program test_exterior_harmonic_tangent_slow
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t, &
                                                 bem_exterior_dtn_t, build_exterior_dtn
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
                                           bem_pair_report_t, assemble_parametric_bem_p2
    use fortfem_exterior_harmonic_fields, only: assemble_exterior_harmonics, &
                                      tangent_exterior_harmonics, evaluate_current_loops
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use parametric_bem_test_geometry, only: test_torus_t
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(test_torus_t) :: torus
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    integer, allocatable :: td(:, :), fd(:, :)
    real(dp),allocatable :: vertices(:,:),nload(:,:),gram(:,:),s(:,:),u(:,:),gt(:,:)
  real(dp), allocatable :: flux(:), load(:), xi(:), eta(:), qw(:), x(:, :), normal(:, :)
    real(dp), allocatable :: potential(:, :, :), field(:, :, :), weight(:), basis(:, :)
    real(dp) :: wire(3, 128), angle, asym, a(3), b(3), n(3), gradient(2), hess(2, 2)
    real(dp) :: total, reference, value, assembly, solve, start, gram_time
    integer :: i, p, t, q, nq, status
    logical :: passed
    passed = .true.
    do i = 1, 128
        angle = 2*acos(-1.0_dp)*real(i - 1, dp)/128
        wire(:, i) = [2*cos(angle), 2*sin(angle), 0.0_dp]
    end do
    call generate_torus_surface_mesh(2.0_dp, 0.6_dp, 8, 16, &
                                     vertices, torus%faces, torus%parameters)
    options%maximum_degree = 80
    start = omp_get_wtime()
    call assemble_parametric_bem_p2(torus,torus%faces,options,op,td,fd, &
        report,status)
    assembly = omp_get_wtime() - start
    write(*,'(a,4i12,es18.9)') 'Tangent pair status failed accepted degree ratio ', &
        status,report%failed_pairs,report%accepted_pairs, &
        report%maximum_degree_used,report%error_ratio
    if (status /= 0) error stop 'tangent P2 BEM assembly'
    start = omp_get_wtime()
    call build_exterior_dtn(op, exterior, status)
    if (status /= 0) error stop 'tangent DtN'
    call exterior%matrix(s, status)
    solve = omp_get_wtime() - start
    start = omp_get_wtime()
    call assemble_exterior_harmonics(torus, torus%faces, td, 16, wire, [1, 129], &
                                     nload, gram, asym, status)
    gram_time = omp_get_wtime() - start
    if (status /= 0) error stop 'tangent raw harmonic'
    call tangent_exterior_harmonics(s, nload, gram, u, gt, status)
  call record(status==0,'tangent harmonic construction solves exterior Neumann problem')
    call record(gt(1, 1) > 0 .and. gt(1, 1) < gram(1, 1), &
            'tangent representative minimizes positive exterior energy at fixed period')
    call record(maxval(abs(matmul(s, u) - nload)) < 1.e-11_dp, &
                'scalar cross energy of tangent representative vanishes')
    allocate (flux(size(op%v, 1)), load(size(s, 1)))
    call exterior%apply(u(:, 1), flux, load, status)
    if (status /= 0) error stop 'tangent flux'
    call triangle_duffy_quadrature(20, xi, eta, qw, status)
    nq = size(qw)
    allocate (x(3, nq*size(torus%faces, 2)), normal(3, nq*size(torus%faces, 2)), &
              weight(nq*size(torus%faces, 2)), basis(3, nq))
    do t = 1, nq
        do i = 1, 3
            call generated_p1_jet(i, xi(t), eta(t), basis(i, t), gradient, hess)
        end do
    end do
    do p = 1, size(torus%faces, 2)
        do t = 1, nq
            q = (p - 1)*nq + t
            call torus%evaluate(p, [xi(t), eta(t)], x(:, q), a, b, status)
            n = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
            normal(:, q) = n/norm2(n); weight(q) = qw(t)*norm2(n)
        end do
    end do
    allocate (potential(3, size(x, 2), 1), field(3, size(x, 2), 1))
    call evaluate_current_loops(wire, [1, 129], x, potential, field, status)
    if (status /= 0) error stop 'tangent independent field samples'
    total = 0; reference = 0
    do p = 1, size(torus%faces, 2)
        do t = 1, nq
            q = (p - 1)*nq + t
            value = dot_product(flux(fd(:, p)), basis(:, t)) + &
                    dot_product(field(:, q, 1), normal(:, q))
            total = total + weight(q)*value**2
            reference = reference + weight(q)*sum(field(:, q, 1)**2)
        end do
    end do
write(*,'(a,es18.9)') 'Exterior tangent normal-field relative L2 ',sqrt(total/reference)
    call record(sqrt(total/reference) < 0.03_dp, &
              'physical tangent harmonic field has less than three percent normal leak')
 write(*,'(a,i4,5es18.9)') 'Exterior scaling threads assembly solve gram G tangentG ', &
        omp_get_max_threads(), assembly, solve, gram_time, gram(1, 1), gt(1, 1)
    call check_summary('Exterior tangent harmonic field and BEM scaling')
    if (.not. passed) error stop 1
contains
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
