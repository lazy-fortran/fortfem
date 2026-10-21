program test_exterior_harmonic_energy
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_elliptic, only: ellip_k, ellip_e
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_exterior_harmonic_fields, only: assemble_exterior_harmonics
    use parametric_bem_test_geometry, only: test_torus_t
    implicit none
    type(test_torus_t) :: torus
    real(dp), allocatable :: vertices(:, :), normal_load(:, :), gram(:, :)
    real(dp) :: wire(3, 128), theta, asym, energy1, energy2, energy3
    integer :: i, status
    logical :: passed
    passed = .true.
    do i = 1, 128
        theta = 2*acos(-1.0_dp)*real(i - 1, dp)/128
        wire(:, i) = [2*cos(theta), 2*sin(theta), 0.0_dp]
    end do
    call generate_torus_surface_mesh(2.0_dp, 0.6_dp, 8, 16, &
                                     vertices, torus%faces, torus%parameters)
    call assemble_exterior_harmonics(torus, torus%faces, torus%faces, 12, &
                                     wire, [1, 129], normal_load, gram, asym, status)
    call record(status == 0, 'general curved boundary harmonic energy construction')
    if (status /= 0) error stop 'exterior Gram construction'
    energy1 = volume_oracle(32, 12.0_dp)
    energy2 = volume_oracle(64, 12.0_dp)
    energy3 = volume_oracle(64, 24.0_dp)
    write (*, '(a,4es18.9)') 'Exterior loop Gram / volume 32 / 64 / tail ', &
        gram(1, 1), energy1, energy2, energy3
   call record(abs(energy1/energy2-1)<0.003_dp,'independent volume rule resolves field')
  call record(abs(energy3/energy2 - 1) < 0.003_dp, 'independent exterior tail is small')
    call record(abs(gram(1, 1)/energy3 - 1) < 0.01_dp, &
           'boundary Green energy agrees with independent elliptic-field volume oracle')
    call record(maxval(abs(normal_load)) > 1.e-5_dp, &
                'non-tangent representative retains nonzero scalar cross energy')
  call record(abs(sum(normal_load))<1.e-8_dp,'closed wire has zero total magnetic flux')
    call check_summary('Finite-energy exterior poloidal harmonic sector')
    if (.not. passed) error stop 1
contains
    real(dp) function volume_oracle(n, extent) result(energy)
        integer, intent(in) :: n
        real(dp), intent(in) :: extent
        real(dp) :: nodes(n), w(n), zn(n), zw(n), r, z, zmin, jac, field(2), t
        real(dp) :: lo, hi
        integer :: section, i, j
        energy = 0
        do section = 1, 3
            lo = 0; hi = 1.4_dp
            if (section == 2) then
                lo = 0; hi = acos(-1.0_dp)
            else if (section == 3) then
                lo = 2.6_dp; hi = extent
            end if
            call gauss_legendre_ab(n, lo, hi, nodes, w)
            do i = 1, n
                r = nodes(i); zmin = 0; jac = 1
                if (section == 2) then
                    t = nodes(i); r = 2 + 0.6_dp*cos(t)
                    zmin = 0.6_dp*sin(t); jac = 0.6_dp*sin(t)
                end if
                call gauss_legendre_ab(n, zmin, extent, zn, zw)
                do j = 1, n
                    z = zn(j)
                    call circular_field(r, z, field)
                    energy = energy + 4*acos(-1.0_dp)*w(i)*zw(j)*jac*r*sum(field**2)
                end do
            end do
        end do
    end function
    subroutine circular_field(r, z, b)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: b(2)
        real(dp) :: d, q, m, k, e
        ! Independent circular-loop closed form, loop radius 2 and mu0 I=1.
        d = (2 + r)**2 + z*z; q = (2 - r)**2 + z*z; m = 8*r/d
        k = ellip_k(m); e = ellip_e(m)
        b(1) = z/(2*acos(-1.0_dp)*r*sqrt(d))*(-k + (4 + r*r + z*z)*e/q)
        b(2) = 1/(2*acos(-1.0_dp)*sqrt(d))*(k + (4 - r*r - z*z)*e/q)
    end subroutine
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
