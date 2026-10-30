program test_scalar_helmholtz_p1_csc
    !! Sparse P1 scalar Helmholtz cells on intervals, triangles and tetrahedra:
    !! the reference-cell matrices sum_d g_d (d_d v, d_d u) - k**2 m (v, u) in
    !! closed form, constants in the kernel of the stiffness part, the mass part
    !! summing to the measure, the JVP against central differences and the VJP
    !! against the JVP by the dot-product test, on jittered meshes with complex
    !! coefficients.
    use check, only: check_condition, check_summary
    use fortfem_block_csc, only: block_csc_add_group, block_csc_finalize, &
        block_csc_init, block_csc_new_matrix, block_csc_pattern_t
    use fortfem_kinds, only: dp
    use fortfem_scalar_helmholtz_p1_csc, only: helmholtz_p1_assemble, &
        helmholtz_p1_assemble_jvp, helmholtz_p1_assemble_vjp
    use fortsparse, only: csc_z_t
    implicit none

    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: cells(:, :)
    integer :: dim

    call check_reference_cells()
    do dim = 1, 3
        call jittered_mesh(dim, vertices, cells)
        call check_invariants(dim)
        call check_derivatives(dim)
    end do
    call check_summary("scalar_helmholtz_p1_csc")

contains

    subroutine check_reference_cells()
        complex(dp), parameter :: g(3) = [(1.3_dp, 0.2_dp), (0.7_dp, -0.4_dp), &
            (1.1_dp, 0.5_dp)]
        complex(dp), parameter :: m = (0.9_dp, 0.3_dp)
        real(dp), parameter :: k = 1.7_dp, h = 0.6_dp
        complex(dp) :: expected(4, 4)
        real(dp) :: grads(3, 4), ones(4, 4), eye(4, 4)
        integer :: i, j

        ones = 1.0_dp
        eye = 0.0_dp
        do i = 1, 4
            eye(i, i) = 1.0_dp
        end do
        ! Interval of length h.
        expected(1:2, 1:2) = g(1)/h*reshape([1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp], &
            [2, 2]) - k**2*m*h/6.0_dp*(ones(1:2, 1:2) + eye(1:2, 1:2))
        call compare_single_cell(reshape([0.4_dp, 0.4_dp + h], [1, 2]), g(1:1), &
            m, k, expected(1:2, 1:2), "interval matrix in closed form")
        ! Unit right triangle (area 1/2), clockwise as well.
        grads(1:2, 1:3) = reshape([-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            1.0_dp], [2, 3])
        do j = 1, 3
            do i = 1, 3
                expected(i, j) = 0.5_dp*sum(g(1:2)*grads(1:2, i)*grads(1:2, j))
            end do
        end do
        expected(1:3, 1:3) = expected(1:3, 1:3) - k**2*m/24.0_dp* &
            (ones(1:3, 1:3) + eye(1:3, 1:3))
        call compare_single_cell(reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp], [2, 3]), g(1:2), m, k, expected(1:3, 1:3), &
            "triangle matrix in closed form")
        call compare_single_cell(reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
            1.0_dp, 0.0_dp], [2, 3]), g(1:2), m, k, &
            expected([1, 3, 2], [1, 3, 2]), "clockwise triangle matrix")
        ! Unit right tetrahedron (volume 1/6).
        grads = reshape([-1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
        do j = 1, 4
            do i = 1, 4
                expected(i, j) = sum(g*grads(:, i)*grads(:, j))/6.0_dp - &
                    k**2*m*(1.0_dp + eye(i, j))/120.0_dp
            end do
        end do
        call compare_single_cell(reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
            [3, 4]), g, m, k, expected, "tetrahedron matrix in closed form")
    end subroutine check_reference_cells

    subroutine compare_single_cell(corners, g, m, k, expected, label)
        real(dp), intent(in) :: corners(:, :), k
        complex(dp), intent(in) :: g(:), m, expected(:, :)
        character(*), intent(in) :: label
        type(block_csc_pattern_t) :: pattern
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: dense(:, :)
        integer :: group, column, entry, n, status

        n = size(corners, 2)
        call build_pattern(n, reshape([(column, column=1, n)], [n, 1]), pattern, &
            group)
        call block_csc_new_matrix(pattern, matrix)
        call helmholtz_p1_assemble(pattern, group, corners, &
            reshape([(column, column=1, n)], [n, 1]), reshape(g, [size(g), 1]), &
            [m], k, matrix%val, status)
        allocate(dense(n, n))
        dense = (0.0_dp, 0.0_dp)
        do column = 1, n
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                dense(matrix%row_idx(entry), column) = matrix%val(entry)
            end do
        end do
        call check_condition(status == 0 .and. maxval(abs(dense - expected)) < &
            1.0e-14_dp*maxval(abs(expected)), label)
    end subroutine compare_single_cell

    subroutine check_invariants(dim)
        !! g = 1, k = 0: rows sum to zero; g = 0, m = 1, k = 1: all entries sum
        !! to minus the measure of the mesh (the jittered meshes fill [0, 1]^dim).
        integer, intent(in) :: dim
        type(block_csc_pattern_t) :: pattern
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: row_sum(:)
        integer :: column, entry, group, nc, status(2)

        nc = size(cells, 2)
        call build_pattern(size(vertices, 2), cells, pattern, group)
        call block_csc_new_matrix(pattern, matrix)
        call helmholtz_p1_assemble(pattern, group, vertices, cells, &
            spread([((1.0_dp, 0.0_dp), column=1, dim)], 2, nc), &
            spread((1.0_dp, 0.0_dp), 1, nc), 0.0_dp, matrix%val, status(1))
        allocate(row_sum(size(vertices, 2)))
        row_sum = (0.0_dp, 0.0_dp)
        do column = 1, matrix%ncol
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row_sum(matrix%row_idx(entry)) = row_sum(matrix%row_idx(entry)) + &
                    matrix%val(entry)
            end do
        end do
        call check_condition(status(1) == 0 .and. maxval(abs(row_sum)) < &
            1.0e-12_dp*maxval(abs(matrix%val)), &
            label(dim, "stiffness annihilates constants"))
        matrix%val = (0.0_dp, 0.0_dp)
        call helmholtz_p1_assemble(pattern, group, vertices, cells, &
            spread([((0.0_dp, 0.0_dp), column=1, dim)], 2, nc), &
            spread((1.0_dp, 0.0_dp), 1, nc), 1.0_dp, matrix%val, status(2))
        call check_condition(status(2) == 0 .and. &
            abs(sum(matrix%val) + 1.0_dp) < 1.0e-13_dp, &
            label(dim, "mass sums to the measure"))
    end subroutine check_invariants

    subroutine check_derivatives(dim)
        integer, intent(in) :: dim
        real(dp), parameter :: h = 1.0e-6_dp, k = 1.3_dp, k_dot = -0.4_dp
        type(block_csc_pattern_t) :: pattern
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: g(:, :), g_dot(:, :), g_bar(:, :)
        complex(dp), allocatable :: m(:), m_dot(:), m_bar(:)
        complex(dp), allocatable :: tangent(:), plus(:), minus(:), values_bar(:)
        real(dp), allocatable :: x_dot(:, :), x_bar(:, :)
        real(dp) :: k_bar, left, right
        integer :: group, i, j, nc, status(4)

        nc = size(cells, 2)
        call build_pattern(size(vertices, 2), cells, pattern, group)
        call block_csc_new_matrix(pattern, matrix)
        allocate(g(dim, nc), g_dot(dim, nc), g_bar(dim, nc))
        allocate(m(nc), m_dot(nc), m_bar(nc))
        allocate(x_dot, mold=vertices)
        allocate(x_bar, mold=vertices)
        do i = 1, nc
            do j = 1, dim
                g(j, i) = cmplx(1.0_dp + 0.1_dp*sin(real(i*(dim + j), dp)), &
                    0.3_dp*cos(real(i + j, dp)), dp)
                g_dot(j, i) = cmplx(cos(0.7_dp*i + 0.1_dp*j), &
                    sin(0.3_dp*i - 0.1_dp*(j + 3)), dp)
            end do
            m(i) = cmplx(1.0_dp + 0.2_dp*cos(1.1_dp*i), 0.4_dp*sin(0.9_dp*i), dp)
            m_dot(i) = cmplx(sin(1.3_dp*i), cos(0.5_dp*i), dp)
        end do
        x_dot = reshape([(0.1_dp*sin(1.7_dp*i), i=1, size(vertices))], &
            shape(vertices))
        allocate(tangent(pattern%nnz), plus(pattern%nnz), minus(pattern%nnz))
        allocate(values_bar(pattern%nnz))
        tangent = (0.0_dp, 0.0_dp)
        plus = (0.0_dp, 0.0_dp)
        minus = (0.0_dp, 0.0_dp)
        call helmholtz_p1_assemble_jvp(pattern, group, vertices, cells, g, m, k, &
            x_dot, g_dot, m_dot, k_dot, tangent, status(1))
        call helmholtz_p1_assemble(pattern, group, vertices + h*x_dot, cells, &
            g + h*g_dot, m + h*m_dot, k + h*k_dot, plus, status(2))
        call helmholtz_p1_assemble(pattern, group, vertices - h*x_dot, cells, &
            g - h*g_dot, m - h*m_dot, k - h*k_dot, minus, status(3))
        call check_condition(all(status(1:3) == 0) .and. &
            maxval(abs(tangent - (plus - minus)/(2.0_dp*h))) < &
            1.0e-7_dp*maxval(abs(tangent)), label(dim, "JVP matches differences"))
        values_bar = [(cmplx(cos(0.37_dp*i), sin(0.61_dp*i), dp), &
            i=1, pattern%nnz)]
        call helmholtz_p1_assemble_vjp(pattern, group, vertices, cells, g, m, k, &
            values_bar, x_bar, g_bar, m_bar, k_bar, status(4))
        left = real(sum(conjg(values_bar)*tangent), dp)
        right = sum(x_bar*x_dot) + real(sum(conjg(g_bar)*g_dot), dp) + &
            real(sum(conjg(m_bar)*m_dot), dp) + k_bar*k_dot
        call check_condition(status(4) == 0 .and. abs(left - right) < &
            1.0e-12_dp*max(1.0_dp, abs(left)), label(dim, "VJP dot-product test"))
    end subroutine check_derivatives

    subroutine build_pattern(n, conn, pattern, group)
        integer, intent(in) :: n, conn(:, :)
        type(block_csc_pattern_t), intent(out) :: pattern
        integer, intent(out) :: group

        call block_csc_init(pattern, n)
        call block_csc_add_group(pattern, conn, group)
        call block_csc_finalize(pattern)
    end subroutine build_pattern

    subroutine jittered_mesh(dim, points, conn)
        !! Unit interval, square or cube, 4 cells per direction, interior
        !! vertices jittered; triangles of both orientations, positive tetrahedra.
        integer, intent(in) :: dim
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: conn(:, :)
        integer, parameter :: n = 5
        integer :: i, j, l, c, v, t, corner(8), tets(4, 6), swap, index(3)
        real(dp) :: e(3, 3), det, jitter(3)

        allocate(points(dim, n**dim))
        do v = 1, n**dim
            i = mod(v - 1, n)
            j = mod((v - 1)/n, n)
            l = (v - 1)/n**2
            index = [i, j, l]
            jitter = 0.05_dp*[sin(1.3_dp*v), cos(2.1_dp*v), sin(0.7_dp*v)]
            points(:, v) = real(index(1:dim), dp)/real(n - 1, dp)
            if (all(index(1:dim) > 0 .and. index(1:dim) < n - 1)) &
                points(:, v) = points(:, v) + jitter(1:dim)
        end do
        select case (dim)
        case (1)
            conn = reshape([((c + j, j=0, 1), c=1, n - 1)], [2, n - 1])
        case (2)
            allocate(conn(3, 2*(n - 1)**2))
            t = 0
            do j = 0, n - 2
                do i = 1, n - 1
                    c = i + j*n
                    t = t + 1
                    conn(:, t) = [c, c + 1, c + n + 1]
                    t = t + 1
                    ! Clockwise.
                    conn(:, t) = [c, c + n, c + n + 1]
                end do
            end do
        case default
            allocate(conn(4, 6*(n - 1)**3))
            tets = reshape([1, 2, 4, 8, 1, 4, 3, 8, 1, 3, 7, 8, 1, 7, 5, 8, &
                1, 5, 6, 8, 1, 6, 2, 8], [4, 6])
            t = 0
            do l = 0, n - 2
                do j = 0, n - 2
                    do i = 1, n - 1
                        c = i + n*(j + n*l)
                        corner = [c, c + 1, c + n, c + n + 1, c + n*n, &
                            c + n*n + 1, c + n*n + n, c + n*n + n + 1]
                        do v = 1, 6
                            t = t + 1
                            conn(:, t) = corner(tets(:, v))
                            do c = 1, 3
                                e(:, c) = points(:, conn(c + 1, t)) - &
                                    points(:, conn(1, t))
                            end do
                            det = dot_product(e(:, 1), [e(2, 2)*e(3, 3) - &
                                e(3, 2)*e(2, 3), e(3, 2)*e(1, 3) - e(1, 2)* &
                                e(3, 3), e(1, 2)*e(2, 3) - e(2, 2)*e(1, 3)])
                            if (det < 0.0_dp) then
                                swap = conn(3, t)
                                conn(3, t) = conn(4, t)
                                conn(4, t) = swap
                            end if
                            c = i + n*(j + n*l)
                        end do
                    end do
                end do
            end do
        end select
    end subroutine jittered_mesh

    function label(dim, text) result(full)
        integer, intent(in) :: dim
        character(*), intent(in) :: text
        character(:), allocatable :: full
        character(4) :: prefix

        write (prefix, "(i0,'D: ')") dim
        full = trim(prefix)//" "//text
    end function label

end program test_scalar_helmholtz_p1_csc
