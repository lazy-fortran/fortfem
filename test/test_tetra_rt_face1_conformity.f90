program test_tetra_rt_face1_conformity
    !! Regression for the RT global dof map on local face 1. Two positively
    !! oriented tetrahedra share their local face 1 (opposite local vertex 2)
    !! in every combination of vertex orderings. A smooth P_k^3 field lies in
    !! the conforming RT_k space, so its global L2 projection must reproduce it
    !! and have a continuous normal trace across the shared face. A face-1
    !! permutation taken from the outward order [1, 4, 3] instead of the RT
    !! moment parametrisation order [1, 3, 4] breaks both for k >= 1.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_rt_arbitrary_order, only: tetra_rt_t, initialize_tetra_rt, &
        evaluate_tetra_rt
    use fortfem_tetra_rt_global_dof_map, only: build_tetra_rt_dof_map, &
        build_tetra_rt_basis_transform
    use fortfem_tetra_piola_maps, only: map_tetra_rt_contravariant
    implicit none

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
        end subroutine dgesv
    end interface

    real(dp), parameter :: tolerance = 1.0e-7_dp
    !> Orderings of the shared face (1, 3, 4) in local slots (1, 3, 4).
    integer, parameter :: even_order(3, 3) = reshape( &
        [1, 3, 4, 3, 4, 1, 4, 1, 3], [3, 3])
    integer, parameter :: odd_order(3, 3) = reshape( &
        [3, 1, 4, 1, 4, 3, 4, 3, 1], [3, 3])
    real(dp), allocatable :: xi(:), eta(:), zeta(:), w(:)
    real(dp) :: coordinates(3, 5), worst_jump, worst_residual
    integer :: degree, ia, ib, ncases, status

    coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    coordinates(:, 2) = [1.1_dp, 0.1_dp, -0.1_dp]
    coordinates(:, 3) = [0.1_dp, 0.9_dp, 0.0_dp]
    coordinates(:, 4) = [-0.1_dp, 0.2_dp, 1.2_dp]
    coordinates(:, 5) = [-0.8_dp, 0.3_dp, 0.2_dp]
    call tetra_duffy_quadrature(8, xi, eta, zeta, w, status)
    call check_condition(status == 0, "tetrahedral quadrature")

    do degree = 1, 3
        worst_jump = 0.0_dp
        worst_residual = 0.0_dp
        ncases = 0
        do ia = 1, 3
            do ib = 1, 3
                call run_case(degree, even_order(:, ia), odd_order(:, ib))
            end do
        end do
        write (*, '(A,I0,A,I0,A,2ES10.2)') "RT_", degree, " cases ", ncases, &
            ": worst normal jump, L2 residual ", worst_jump, worst_residual
        call check_condition(ncases == 9, "all face-1 vertex orderings run")
        call check_condition(worst_jump < tolerance, &
            "RT global dofs give continuous normal trace across shared face 1")
        call check_condition(worst_residual < tolerance, &
            "global RT projection reproduces a smooth P_k field")
    end do
    call check_summary("Tetrahedral RT face-1 conformity")

contains

    subroutine run_case(k, order_a, order_b)
        integer, intent(in) :: k, order_a(3), order_b(3)
        type(tetra_rt_t) :: basis
        integer, allocatable :: faces(:, :), global_dofs(:, :)
        integer, allocatable :: orientations(:, :), permutations(:, :, :)
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: transform(:, :, :), mass(:, :), rhs(:, :)
        real(dp) :: jac(3, 3), det(2), vertices(3, 4, 2), field(3)
        real(dp) :: phi(3, (k + 1)*(k + 2)*(k + 4)/2), dv(size(phi, 2))
        real(dp) :: ref_values(3, size(phi, 2)), ref_div(size(phi, 2))
        real(dp) :: psi(3, size(phi, 2)), x(3), residual_sq, trace(2), normal(3)
        real(dp) :: lambda(3), point(3), ref_point(3), jinv(3, 3), uh(3)
        integer :: tetrahedra(4, 2), nrt, ng, t, q, i, j, info, p
        integer :: ref_status

        tetrahedra(:, 1) = [order_a(1), 2, order_a(2), order_a(3)]
        tetrahedra(:, 2) = [order_b(1), 5, order_b(2), order_b(3)]
        do t = 1, 2
            vertices(:, :, t) = coordinates(:, tetrahedra(:, t))
            do i = 1, 3
                jac(:, i) = vertices(:, i + 1, t) - vertices(:, 1, t)
            end do
            det(t) = determinant(jac)
            if (.not. (det(t) > 0.0_dp)) then
                call check_condition(.false., "tetrahedron is positively oriented")
                return
            end if
        end do
        ncases = ncases + 1
        call initialize_tetra_rt(k, basis, status)
        nrt = basis%dof_count
        call build_tetra_rt_dof_map(k, tetrahedra, faces, global_dofs, &
            orientations, permutations, status)
        call check_condition(status == 0 .and. size(faces, 2) == 7, &
            "two tetrahedra share exactly one face")
        allocate (transform(nrt, nrt, 2))
        do t = 1, 2
            call build_tetra_rt_basis_transform(k, orientations(:, t), &
                permutations(:, :, t), transform(:, :, t), status)
        end do
        ng = maxval(global_dofs)
        allocate (mass(ng, ng), rhs(ng, 1), ipiv(ng))
        mass = 0.0_dp
        rhs = 0.0_dp
        do t = 1, 2
            do i = 1, 3
                jac(:, i) = vertices(:, i + 1, t) - vertices(:, 1, t)
            end do
            do q = 1, size(w)
                call evaluate_tetra_rt(basis, [xi(q), eta(q), zeta(q)], &
                    ref_values, ref_div, status)
                call map_tetra_rt_contravariant(jac, ref_values, ref_div, phi, dv, &
                    ref_status)
                psi = matmul(phi, transform(:, :, t))
                x = vertices(:, 1, t) + matmul(jac, [xi(q), eta(q), zeta(q)])
                field = smooth_field(k, x)
                do j = 1, nrt
                    rhs(global_dofs(j, t), 1) = rhs(global_dofs(j, t), 1) &
                        + w(q)*det(t)*dot_product(psi(:, j), field)
                    do i = 1, nrt
                        mass(global_dofs(i, t), global_dofs(j, t)) = &
                            mass(global_dofs(i, t), global_dofs(j, t)) &
                            + w(q)*det(t)*dot_product(psi(:, i), psi(:, j))
                    end do
                end do
            end do
        end do
        call dgesv(ng, 1, mass, ng, ipiv, rhs, ng, info)
        call check_condition(info == 0, "global RT mass matrix is solvable")

        residual_sq = 0.0_dp
        do t = 1, 2
            do i = 1, 3
                jac(:, i) = vertices(:, i + 1, t) - vertices(:, 1, t)
            end do
            do q = 1, size(w)
                call evaluate_tetra_rt(basis, [xi(q), eta(q), zeta(q)], &
                    ref_values, ref_div, status)
                call map_tetra_rt_contravariant(jac, ref_values, ref_div, phi, dv, &
                    ref_status)
                psi = matmul(phi, transform(:, :, t))
                x = vertices(:, 1, t) + matmul(jac, [xi(q), eta(q), zeta(q)])
                uh = 0.0_dp
                do j = 1, nrt
                    uh = uh + rhs(global_dofs(j, t), 1)*psi(:, j)
                end do
                residual_sq = residual_sq + w(q)*det(t) &
                    *sum((uh - smooth_field(k, x))**2)
            end do
        end do
        worst_residual = max(worst_residual, sqrt(residual_sq))

        ! Normal trace of both cells at points of the shared face (1, 3, 4).
        normal = cross(coordinates(:, 3) - coordinates(:, 1), &
            coordinates(:, 4) - coordinates(:, 1))
        normal = normal/norm2(normal)
        do p = 1, 4
            lambda = [0.2_dp*p, 0.1_dp + 0.15_dp*p, 0.05_dp]
            lambda = lambda/sum(lambda)
            point = lambda(1)*coordinates(:, 1) + lambda(2)*coordinates(:, 3) &
                + lambda(3)*coordinates(:, 4)
            do t = 1, 2
                do i = 1, 3
                    jac(:, i) = vertices(:, i + 1, t) - vertices(:, 1, t)
                end do
                call invert3(jac, jinv)
                ref_point = matmul(jinv, point - vertices(:, 1, t))
                call evaluate_tetra_rt(basis, ref_point, ref_values, ref_div, &
                    status)
                call map_tetra_rt_contravariant(jac, ref_values, ref_div, phi, dv, &
                    ref_status)
                psi = matmul(phi, transform(:, :, t))
                uh = 0.0_dp
                do j = 1, nrt
                    uh = uh + rhs(global_dofs(j, t), 1)*psi(:, j)
                end do
                trace(t) = dot_product(uh, normal)
            end do
            worst_jump = max(worst_jump, abs(trace(1) - trace(2)))
        end do
    end subroutine run_case

    function smooth_field(k, x) result(f)
        integer, intent(in) :: k
        real(dp), intent(in) :: x(3)
        real(dp) :: f(3)
        integer :: a, b, c, n

        f = 0.0_dp
        n = 0
        do a = 0, k
            do b = 0, k - a
                do c = 0, k - a - b
                    n = n + 1
                    f = f + x(1)**a*x(2)**b*x(3)**c &
                        *[sin(1.3_dp*n), cos(0.7_dp*n + 1.0_dp), sin(2.1_dp*n + 0.5_dp)]
                end do
            end do
        end do
    end function smooth_field

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross

    pure function determinant(m) result(d)
        real(dp), intent(in) :: m(3, 3)
        real(dp) :: d

        d = dot_product(m(:, 1), cross(m(:, 2), m(:, 3)))
    end function determinant

    pure subroutine invert3(m, minv)
        real(dp), intent(in) :: m(3, 3)
        real(dp), intent(out) :: minv(3, 3)

        minv(1, :) = cross(m(:, 2), m(:, 3))
        minv(2, :) = cross(m(:, 3), m(:, 1))
        minv(3, :) = cross(m(:, 1), m(:, 2))
        minv = minv/determinant(m)
    end subroutine invert3

end program test_tetra_rt_face1_conformity
