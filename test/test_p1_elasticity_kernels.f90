program test_p1_elasticity_kernels
    !! Generated vector P1 elasticity kernels: the 6x6 cell matrix against the
    !! handwritten assembly of elasticity_p1_2d_core (both orientations), symmetry
    !! and rigid-body null space; the 4x4 impedance edge matrix against the
    !! absorbing block of elasticity_planar_acoustic_dtn_2d; JVP against central
    !! differences and VJP against JVP by the (bilinear) dot-product test.
    use check, only: check_condition, check_summary
    use fortfem_elasticity_p1_2d_core, only: assemble_elasticity_dynamic_p1_2d
    use fortfem_generated_p1_elastic_impedance_edge, only: &
        generated_p1_elastic_impedance_edge
    use fortfem_generated_p1_elastic_impedance_edge_jvp, only: &
        generated_p1_elastic_impedance_edge_jvp
    use fortfem_generated_p1_elastic_impedance_edge_vjp, only: &
        generated_p1_elastic_impedance_edge_vjp
    use fortfem_generated_p1_elasticity_cell, only: generated_p1_elasticity_cell
    use fortfem_generated_p1_elasticity_cell_jvp, only: &
        generated_p1_elasticity_cell_jvp
    use fortfem_generated_p1_elasticity_cell_vjp, only: &
        generated_p1_elasticity_cell_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer :: trial

    call random_seed_fixed()
    do trial = 1, 20
        call check_cell_against_core(trial)
    end do
    call check_cell_structure()
    call check_edge_against_absorbing_block()
    call check_cell_derivatives()
    call check_edge_derivatives()
    call check_summary("p1_elasticity_kernels")

contains

    subroutine check_cell_against_core(trial)
        integer, intent(in) :: trial
        real(dp) :: vertices(2, 3), r(4), omega, nu, rho, orientation
        real(dp) :: a, b, c, d
        complex(dp) :: young, lambda, mu, reference(6, 6), cell(6, 6)
        integer :: status
        character(48) :: label

        call random_number(vertices)
        call affine(vertices, a, b, c, d)
        ! Odd trials clockwise, even trials counterclockwise.
        if ((mod(trial, 2) == 1) .eqv. (a*d - b*c > 0.0_dp)) &
            vertices(:, [2, 3]) = vertices(:, [3, 2])
        call random_number(r)
        omega = 0.5_dp + 3.0_dp*r(1)
        nu = 0.45_dp*r(2)
        rho = 0.5_dp + r(3)
        young = cmplx(1.0_dp + r(4), 0.1_dp*r(4), dp)
        call assemble_elasticity_dynamic_p1_2d(vertices, reshape([1, 2, 3], &
            [3, 1]), omega, young, nu, rho, reference, status)
        mu = young/(2.0_dp*(1.0_dp + nu))
        lambda = young*nu/((1.0_dp + nu)*(1.0_dp - 2.0_dp*nu))
        call affine(vertices, a, b, c, d)
        orientation = sign(1.0_dp, a*d - b*c)
        call generated_p1_elasticity_cell(a, b, c, d, orientation, lambda, mu, &
            cmplx(rho*omega**2, 0.0_dp, dp), cell)
        write (label, "(a,i0)") "cell equals core local matrix, trial ", trial
        call check_condition(status == 0 .and. maxval(abs(cell - reference)) <= &
            1.0e-12_dp*maxval(abs(reference)), trim(label))
    end subroutine check_cell_against_core

    subroutine check_cell_structure()
        real(dp) :: vertices(2, 3), rigid(6, 3), a, b, c, d
        complex(dp) :: lambda, mu, inertia, cell(6, 6), static(6, 6)
        integer :: k

        call random_number(vertices)
        call affine(vertices, a, b, c, d)
        lambda = (1.3_dp, 0.2_dp)
        mu = (0.7_dp, -0.1_dp)
        inertia = (2.1_dp, 0.4_dp)
        call generated_p1_elasticity_cell(a, b, c, d, sign(1.0_dp, a*d - b*c), &
            lambda, mu, inertia, cell)
        call check_condition(maxval(abs(cell - transpose(cell))) <= &
            1.0e-14_dp*maxval(abs(cell)), "cell is complex symmetric")
        call generated_p1_elasticity_cell(a, b, c, d, sign(1.0_dp, a*d - b*c), &
            lambda, mu, (0.0_dp, 0.0_dp), static)
        do k = 1, 3
            rigid(2*k - 1:2*k, 1) = [1.0_dp, 0.0_dp]
            rigid(2*k - 1:2*k, 2) = [0.0_dp, 1.0_dp]
            rigid(2*k - 1:2*k, 3) = [-vertices(2, k), vertices(1, k)]
        end do
        call check_condition(maxval(abs(matmul(static, rigid))) <= &
            1.0e-12_dp*maxval(abs(static)), "rigid motions span the static null space")
    end subroutine check_cell_structure

    subroutine check_edge_against_absorbing_block()
        real(dp) :: r(5), length, normal(2), omega, angle
        complex(dp) :: edge(4, 4), zp, zs
        integer :: k
        logical :: ok

        ok = .true.
        do k = 1, 10
            call random_number(r)
            length = 0.1_dp + r(1)
            angle = 6.0_dp*r(2)
            normal = [cos(angle), sin(angle)]
            zp = cmplx(1.0_dp + r(3), 0.0_dp, dp)
            zs = cmplx(0.3_dp + r(4), 0.0_dp, dp)
            omega = 0.5_dp + 4.0_dp*r(5)
            call generated_p1_elastic_impedance_edge(length, normal(1), normal(2), &
                zp, zs, omega, edge)
            ok = ok .and. maxval(abs(edge - absorbing_block(length, normal, zp, zs, &
                omega))) <= 1.0e-14_dp*maxval(abs(edge))
        end do
        call check_condition(ok, "edge equals absorbing block for real impedances")
        zp = (1.5_dp, 0.2_dp)
        zs = (0.6_dp, -0.3_dp)
        call generated_p1_elastic_impedance_edge(length, normal(1), normal(2), zp, &
            zs, omega, edge)
        call check_condition(maxval(abs(edge - absorbing_block(length, normal, zp, &
            zs, omega))) <= 1.0e-14_dp*maxval(abs(edge)) .and. &
            maxval(abs(edge - transpose(edge))) <= 1.0e-15_dp*maxval(abs(edge)), &
            "edge keeps complex impedances and is symmetric")
    end subroutine check_edge_against_absorbing_block

    function absorbing_block(length, normal, p_impedance, s_impedance, omega) &
            result(local)
        !! Transcription of assemble_absorbing_boundary (elasticity_planar_
        !! acoustic_dtn_2d.f90) without its real() on the impedances.
        real(dp), intent(in) :: length, normal(2), omega
        complex(dp), intent(in) :: p_impedance, s_impedance
        complex(dp) :: local(4, 4), coefficient(2, 2)
        real(dp) :: edge_mass(2, 2)
        integer :: component, first, second

        edge_mass = length/6.0_dp*reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        coefficient = (0.0_dp, 0.0_dp)
        do component = 1, 2
            coefficient(component, component) = s_impedance
        end do
        do first = 1, 2
            do second = 1, 2
                coefficient(first, second) = coefficient(first, second) + &
                    (p_impedance - s_impedance)*normal(first)*normal(second)
            end do
        end do
        do first = 1, 2
            do second = 1, 2
                local(2*first - 1:2*first, 2*second - 1:2*second) = &
                    cmplx(0.0_dp, omega*edge_mass(first, second), dp)*coefficient
            end do
        end do
    end function absorbing_block

    subroutine check_cell_derivatives()
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp) :: geometry(4), geometry_dot(4), orientation, r(6)
        complex(dp) :: material(3), material_dot(3), bar(6, 6)
        complex(dp) :: tangent(6, 6), plus(6, 6), minus(6, 6), cotangent(7)
        complex(dp) :: left, right

        call random_number(geometry)
        geometry = geometry + [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        orientation = sign(1.0_dp, geometry(1)*geometry(4) - geometry(2)*geometry(3))
        material = [(1.2_dp, 0.1_dp), (0.8_dp, -0.05_dp), (3.0_dp, 0.3_dp)]
        call random_number(geometry_dot)
        call random_number(r)
        material_dot = cmplx(r(1:3) - 0.5_dp, r(4:6) - 0.5_dp, dp)
        call generated_p1_elasticity_cell_jvp(geometry(1), geometry(2), &
            geometry(3), geometry(4), orientation, material(1), material(2), &
            material(3), geometry_dot(1), geometry_dot(2), geometry_dot(3), &
            geometry_dot(4), material_dot(1), material_dot(2), material_dot(3), &
            tangent)
        call cell_at(geometry + h*geometry_dot, orientation, &
            material + h*material_dot, plus)
        call cell_at(geometry - h*geometry_dot, orientation, &
            material - h*material_dot, minus)
        call check_condition(maxval(abs((plus - minus)/(2.0_dp*h) - tangent)) <= &
            1.0e-7_dp*maxval(abs(tangent)), "cell JVP matches central differences")

        bar = cmplx(random_matrix(6), random_matrix(6), dp)
        call generated_p1_elasticity_cell_vjp(geometry(1), geometry(2), &
            geometry(3), geometry(4), orientation, material(1), material(2), &
            material(3), bar, cotangent(1), cotangent(2), cotangent(3), &
            cotangent(4), cotangent(5), cotangent(6), cotangent(7))
        left = sum(bar*tangent)
        right = sum(cotangent(1:4)*geometry_dot) + sum(cotangent(5:7)*material_dot)
        call check_condition(abs(left - right) <= 1.0e-12_dp*abs(left), &
            "cell VJP is the transpose of the JVP")
    end subroutine check_cell_derivatives

    subroutine check_edge_derivatives()
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp) :: geometry(4), geometry_dot(4), r(4)
        complex(dp) :: impedance(2), impedance_dot(2), bar(4, 4)
        complex(dp) :: tangent(4, 4), plus(4, 4), minus(4, 4), cotangent(6)
        complex(dp) :: left, right

        ! geometry = (length, normal_1, normal_2, angular frequency)
        geometry = [0.7_dp, 0.6_dp, 0.8_dp, 2.5_dp]
        impedance = [(1.4_dp, 0.1_dp), (0.5_dp, -0.2_dp)]
        call random_number(geometry_dot)
        call random_number(r)
        impedance_dot = cmplx(r(1:2) - 0.5_dp, r(3:4) - 0.5_dp, dp)
        call generated_p1_elastic_impedance_edge_jvp(geometry(1), geometry(2), &
            geometry(3), impedance(1), impedance(2), geometry(4), &
            geometry_dot(1), geometry_dot(2), geometry_dot(3), impedance_dot(1), &
            impedance_dot(2), geometry_dot(4), tangent)
        call edge_at(geometry + h*geometry_dot, impedance + h*impedance_dot, plus)
        call edge_at(geometry - h*geometry_dot, impedance - h*impedance_dot, minus)
        call check_condition(maxval(abs((plus - minus)/(2.0_dp*h) - tangent)) <= &
            1.0e-7_dp*maxval(abs(tangent)), "edge JVP matches central differences")

        bar = cmplx(random_matrix(4), random_matrix(4), dp)
        call generated_p1_elastic_impedance_edge_vjp(geometry(1), geometry(2), &
            geometry(3), impedance(1), impedance(2), geometry(4), bar, &
            cotangent(1), cotangent(2), cotangent(3), cotangent(4), cotangent(5), &
            cotangent(6))
        left = sum(bar*tangent)
        right = sum(cotangent([1, 2, 3, 6])*geometry_dot) + &
            sum(cotangent(4:5)*impedance_dot)
        call check_condition(abs(left - right) <= 1.0e-12_dp*abs(left), &
            "edge VJP is the transpose of the JVP")
    end subroutine check_edge_derivatives

    subroutine cell_at(g, orientation, m, cell)
        real(dp), intent(in) :: g(4), orientation
        complex(dp), intent(in) :: m(3)
        complex(dp), intent(out) :: cell(6, 6)

        call generated_p1_elasticity_cell(g(1), g(2), g(3), g(4), orientation, &
            m(1), m(2), m(3), cell)
    end subroutine cell_at

    subroutine edge_at(g, z, edge)
        real(dp), intent(in) :: g(4)
        complex(dp), intent(in) :: z(2)
        complex(dp), intent(out) :: edge(4, 4)

        call generated_p1_elastic_impedance_edge(g(1), g(2), g(3), z(1), z(2), &
            g(4), edge)
    end subroutine edge_at

    pure subroutine affine(vertices, a, b, c, d)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp), intent(out) :: a, b, c, d

        a = vertices(1, 2) - vertices(1, 1)
        b = vertices(2, 2) - vertices(2, 1)
        c = vertices(1, 3) - vertices(1, 1)
        d = vertices(2, 3) - vertices(2, 1)
    end subroutine affine

    function random_matrix(n) result(matrix)
        integer, intent(in) :: n
        real(dp) :: matrix(n, n)

        call random_number(matrix)
        matrix = 2.0_dp*matrix - 1.0_dp
    end function random_matrix

    subroutine random_seed_fixed()
        integer, allocatable :: state(:)
        integer :: k, n

        call random_seed(size=n)
        allocate(state(n))
        state = [(1009*k + 17, k=1, n)]
        call random_seed(put=state)
    end subroutine random_seed_fixed

end program test_p1_elasticity_kernels
