program test_elasticity_p1_csc_2d
    !! Sparse vector P1 elasticity assembly: equal to the dense assembly of
    !! elasticity_p1_2d_core on a jittered mesh with mixed orientations, JVP
    !! against central differences, VJP against JVP by the dot-product test, and
    !! thread-count invariance of values and vertex cotangents.
    use check, only: check_condition, check_summary
    use fortfem_block_csc, only: block_csc_new_matrix, block_csc_pattern_t
    use fortfem_elasticity_p1_2d_core, only: assemble_elasticity_dynamic_p1_2d
    use fortfem_elasticity_p1_csc_2d, only: elasticity_p1_assemble, &
        elasticity_p1_assemble_jvp, elasticity_p1_assemble_vjp, &
        elasticity_p1_pattern
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_is_valid, csc_z_t
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)

    call jittered_mesh(7, vertices, triangles)
    call check_against_dense()
    call check_derivatives()
    call check_thread_invariance()
    call check_summary("elasticity_p1_csc_2d")

contains

    subroutine check_against_dense()
        type(block_csc_pattern_t) :: pattern
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: dense(:, :), sparse(:, :)
        complex(dp) :: young, lambda, mu
        real(dp) :: omega, nu, rho
        integer :: cell_group, column, entry, n, nt, status, dense_status

        nt = size(triangles, 2)
        n = 2*size(vertices, 2)
        omega = 1.7_dp
        nu = 0.31_dp
        rho = 1.3_dp
        young = (2.0_dp, 0.15_dp)
        mu = young/(2.0_dp*(1.0_dp + nu))
        lambda = young*nu/((1.0_dp + nu)*(1.0_dp - 2.0_dp*nu))
        call elasticity_p1_pattern(size(vertices, 2), triangles, pattern, cell_group)
        call block_csc_new_matrix(pattern, matrix)
        call elasticity_p1_assemble(pattern, cell_group, vertices, triangles, &
            spread(lambda, 1, nt), spread(mu, 1, nt), &
            spread(cmplx(rho*omega**2, 0.0_dp, dp), 1, nt), matrix%val, status)
        allocate(dense(n, n), sparse(n, n))
        call assemble_elasticity_dynamic_p1_2d(vertices, triangles, omega, young, &
            nu, rho, dense, dense_status)
        sparse = (0.0_dp, 0.0_dp)
        do column = 1, n
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                sparse(matrix%row_idx(entry), column) = matrix%val(entry)
            end do
        end do
        call check_condition(status == 0 .and. dense_status == 0 .and. &
            csc_is_valid(matrix), "assembly succeeds with a valid CSC matrix")
        call check_condition(maxval(abs(sparse - dense)) <= &
            1.0e-12_dp*maxval(abs(dense)), "sparse assembly equals dense core")
        call check_condition(matrix%nnz >= count(abs(dense) > 0.0_dp) .and. &
            matrix%nnz <= 4*7*size(vertices, 2), &
            "pattern holds the nonzeros with O(vertices) entries")
    end subroutine check_against_dense

    subroutine check_derivatives()
        real(dp), parameter :: h = 1.0e-6_dp
        type(block_csc_pattern_t) :: pattern
        complex(dp), allocatable :: lambda(:), mu(:), inertia(:)
        complex(dp), allocatable :: lambda_dot(:), mu_dot(:), inertia_dot(:)
        complex(dp), allocatable :: lambda_bar(:), mu_bar(:), inertia_bar(:)
        complex(dp), allocatable :: tangent(:), plus(:), minus(:), values_bar(:)
        real(dp), allocatable :: vertices_dot(:, :), vertices_bar(:, :), r(:)
        real(dp) :: left, right
        integer :: cell_group, nt, status(5)

        nt = size(triangles, 2)
        call elasticity_p1_pattern(size(vertices, 2), triangles, pattern, cell_group)
        call materials(nt, 1, lambda, mu, inertia)
        call materials(nt, 2, lambda_dot, mu_dot, inertia_dot)
        allocate(vertices_dot, mold=vertices)
        call random_number(vertices_dot)
        vertices_dot = vertices_dot - 0.5_dp
        allocate(tangent(pattern%nnz), plus(pattern%nnz), minus(pattern%nnz))
        tangent = (0.0_dp, 0.0_dp)
        plus = (0.0_dp, 0.0_dp)
        minus = (0.0_dp, 0.0_dp)
        call elasticity_p1_assemble_jvp(pattern, cell_group, vertices, triangles, &
            lambda, mu, inertia, vertices_dot, lambda_dot, mu_dot, inertia_dot, &
            tangent, status(1))
        call elasticity_p1_assemble(pattern, cell_group, vertices + h*vertices_dot, &
            triangles, lambda + h*lambda_dot, mu + h*mu_dot, &
            inertia + h*inertia_dot, plus, status(2))
        call elasticity_p1_assemble(pattern, cell_group, vertices - h*vertices_dot, &
            triangles, lambda - h*lambda_dot, mu - h*mu_dot, &
            inertia - h*inertia_dot, minus, status(3))
        call check_condition(all(status(1:3) == 0) .and. &
            maxval(abs((plus - minus)/(2.0_dp*h) - tangent)) <= &
            1.0e-7_dp*maxval(abs(tangent)), "assembly JVP matches central differences")

        allocate(values_bar(pattern%nnz), r(2*pattern%nnz))
        call random_number(r)
        values_bar = cmplx(r(:pattern%nnz) - 0.5_dp, r(pattern%nnz + 1:) - 0.5_dp, dp)
        allocate(vertices_bar, mold=vertices)
        allocate(lambda_bar(nt), mu_bar(nt), inertia_bar(nt))
        call elasticity_p1_assemble_vjp(pattern, cell_group, vertices, triangles, &
            lambda, mu, inertia, values_bar, vertices_bar, lambda_bar, mu_bar, &
            inertia_bar, status(4))
        left = real(sum(conjg(values_bar)*tangent), dp)
        right = sum(vertices_bar*vertices_dot) + real(sum(conjg(lambda_bar)* &
            lambda_dot) + sum(conjg(mu_bar)*mu_dot) + &
            sum(conjg(inertia_bar)*inertia_dot), dp)
        call check_condition(status(4) == 0 .and. abs(left - right) <= &
            1.0e-12_dp*sum(abs(values_bar)*abs(tangent)), &
            "assembly VJP is the adjoint of the JVP")
        ! The material part is complex linear: the full complex identity holds.
        tangent = (0.0_dp, 0.0_dp)
        call elasticity_p1_assemble_jvp(pattern, cell_group, vertices, triangles, &
            lambda, mu, inertia, 0.0_dp*vertices_dot, lambda_dot, mu_dot, &
            inertia_dot, tangent, status(5))
        call check_condition(abs(sum(conjg(values_bar)*tangent) - &
            (sum(conjg(lambda_bar)*lambda_dot) + sum(conjg(mu_bar)*mu_dot) + &
            sum(conjg(inertia_bar)*inertia_dot))) <= &
            1.0e-12_dp*sum(abs(values_bar)*abs(tangent)), &
            "material cotangents are conjugate-linear")
    end subroutine check_derivatives

    subroutine check_thread_invariance()
        type(block_csc_pattern_t) :: pattern
        complex(dp), allocatable :: lambda(:), mu(:), inertia(:), values_bar(:)
        complex(dp), allocatable :: one(:), many(:), lb(:), mb(:), ib(:)
        real(dp), allocatable :: big_vertices(:, :), bar_one(:, :), bar_many(:, :)
        integer, allocatable :: big_triangles(:, :)
        integer :: cell_group, nt, status(4), threads

        threads = 1
        !$ threads = omp_get_max_threads()
        call jittered_mesh(80, big_vertices, big_triangles)
        nt = size(big_triangles, 2)
        call materials(nt, 3, lambda, mu, inertia)
        call elasticity_p1_pattern(size(big_vertices, 2), big_triangles, pattern, &
            cell_group)
        allocate(one(pattern%nnz), many(pattern%nnz), values_bar(pattern%nnz))
        allocate(bar_one, bar_many, mold=big_vertices)
        allocate(lb(nt), mb(nt), ib(nt))
        one = (0.0_dp, 0.0_dp)
        many = (0.0_dp, 0.0_dp)
        values_bar = (1.0_dp, -0.5_dp)
        !$ call omp_set_num_threads(1)
        call elasticity_p1_assemble(pattern, cell_group, big_vertices, &
            big_triangles, lambda, mu, inertia, one, status(1))
        call elasticity_p1_assemble_vjp(pattern, cell_group, big_vertices, &
            big_triangles, lambda, mu, inertia, values_bar, bar_one, lb, mb, ib, &
            status(2))
        !$ call omp_set_num_threads(max(4, threads))
        call elasticity_p1_assemble(pattern, cell_group, big_vertices, &
            big_triangles, lambda, mu, inertia, many, status(3))
        call elasticity_p1_assemble_vjp(pattern, cell_group, big_vertices, &
            big_triangles, lambda, mu, inertia, values_bar, bar_many, lb, mb, ib, &
            status(4))
        !$ call omp_set_num_threads(threads)
        call check_condition(all(status == 0) .and. all(one == many) .and. &
            all(bar_one == bar_many), "values and cotangents independent of threads")
    end subroutine check_thread_invariance

    subroutine materials(nt, seed, lambda, mu, inertia)
        integer, intent(in) :: nt, seed
        complex(dp), allocatable, intent(out) :: lambda(:), mu(:), inertia(:)
        real(dp), allocatable :: r(:, :)
        integer :: e, k

        allocate(r(6, nt), lambda(nt), mu(nt), inertia(nt))
        do e = 1, nt
            r(:, e) = [(modulo(0.6180339887_dp*(e*6 + k + 7*seed), 1.0_dp), k=1, 6)]
        end do
        lambda = cmplx(1.0_dp + r(1, :), 0.1_dp*r(2, :), dp)
        mu = cmplx(0.5_dp + r(3, :), -0.1_dp*r(4, :), dp)
        inertia = cmplx(2.0_dp + r(5, :), 0.2_dp*r(6, :), dp)
    end subroutine materials

    subroutine jittered_mesh(nx, vertices, triangles)
        !! Unit square, nx x nx vertices moved by up to 20% of the spacing, with
        !! the vertex order of every third triangle reversed (clockwise).
        integer, intent(in) :: nx
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: triangles(:, :)
        integer :: i, j, k, node(4)
        real(dp) :: hx

        hx = 1.0_dp/real(nx - 1, dp)
        allocate(vertices(2, nx*nx), triangles(3, 2*(nx - 1)**2))
        do j = 1, nx
            do i = 1, nx
                k = i + (j - 1)*nx
                vertices(:, k) = hx*[real(i - 1, dp), real(j - 1, dp)]
                if (i > 1 .and. i < nx .and. j > 1 .and. j < nx) &
                    vertices(:, k) = vertices(:, k) + 0.2_dp*hx*[ &
                    sin(12.9898_dp*k), cos(78.233_dp*k)]
            end do
        end do
        k = 0
        do j = 1, nx - 1
            do i = 1, nx - 1
                node = [i + (j - 1)*nx, i + 1 + (j - 1)*nx, i + 1 + j*nx, i + j*nx]
                triangles(:, k + 1) = [node(1), node(2), node(3)]
                triangles(:, k + 2) = [node(1), node(3), node(4)]
                k = k + 2
            end do
        end do
        triangles(2:3, ::3) = triangles(3:2:-1, ::3)
    end subroutine jittered_mesh

end program test_elasticity_p1_csc_2d
