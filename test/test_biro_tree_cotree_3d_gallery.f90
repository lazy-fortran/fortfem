program test_biro_tree_cotree_3d_gallery
    !! Independent topology and manufactured-field oracle for the 3-D gallery.
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_tree_cotree_prolongation, &
        build_tree_cotree_gauge, reduce_tree_cotree_dense_system, &
        tree_cotree_gauge_edges, tree_cotree_gauge_t
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nv = 8, ne = 12, nf = 6
    integer, parameter :: incidence(nv, ne) = reshape([ &
        -1, 1, 0, 0, 0, 0, 0, 0, &
        0, -1, 1, 0, 0, 0, 0, 0, &
        0, 0, -1, 1, 0, 0, 0, 0, &
        1, 0, 0, -1, 0, 0, 0, 0, &
        0, 0, 0, 0, -1, 1, 0, 0, &
        0, 0, 0, 0, 0, -1, 1, 0, &
        0, 0, 0, 0, 0, 0, -1, 1, &
        0, 0, 0, 0, 1, 0, 0, -1, &
        -1, 0, 0, 0, 1, 0, 0, 0, &
        0, -1, 0, 0, 0, 1, 0, 0, &
        0, 0, -1, 0, 0, 0, 1, 0, &
        0, 0, 0, -1, 0, 0, 0, 1], &
        [nv, ne])
    integer, parameter :: face_edge(nf, ne) = reshape([ &
        1, 0, 1, 0, 0, 0, &
        1, 0, 0, 1, 0, 0, &
        1, 0, 0, 0, 1, 0, &
        1, 0, 0, 0, 0, 1, &
        0, 1, -1, 0, 0, 0, &
        0, 1, 0, -1, 0, 0, &
        0, 1, 0, 0, -1, 0, &
        0, 1, 0, 0, 0, -1, &
        0, 0, -1, 0, 0, 1, &
        0, 0, 1, -1, 0, 0, &
        0, 0, 0, 1, -1, 0, &
        0, 0, 0, 0, 1, -1], [nf, ne])
    real(dp), parameter :: x(nv) = &
        [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: y(nv) = &
        [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: z(nv) = &
        [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: vertex_potential(nv) = &
        [-0.30_dp, 0.18_dp, 0.65_dp, -0.22_dp, 0.44_dp, &
        0.72_dp, 0.31_dp, -0.11_dp]
    integer, parameter :: edge_nodes(2, ne) = reshape([ &
        1, 2, 2, 3, 3, 4, 4, 1, 5, 6, 6, 7, 7, 8, 8, 5, &
        1, 5, 2, 6, 3, 7, 4, 8], [2, ne])
    real(dp) :: curl_incidence(nf, ne), gram(ne, ne), edge_target(ne)
    real(dp) :: rhs(ne), midpoint(3), tangent(3), potential(3)
    real(dp) :: discrete_gradient(ne)
    real(dp), allocatable :: reduced_matrix(:, :), reduced_rhs(:)
    real(dp), allocatable :: reduced(:), prolonged(:)
    integer, allocatable :: tree(:), cotree(:)
    type(tree_cotree_gauge_t) :: gauge
    integer :: edge, first, second, status

    curl_incidence = real(face_edge, dp)
    gram = matmul(transpose(curl_incidence), curl_incidence)
    do edge = 1, ne
        first = edge_nodes(1, edge)
        second = edge_nodes(2, edge)
        midpoint = 0.5_dp*[x(first) + x(second), y(first) + y(second), &
            z(first) + z(second)]
        tangent = [x(second) - x(first), y(second) - y(first), &
            z(second) - z(first)]
        potential = [-0.30_dp*midpoint(2) + 0.10_dp*midpoint(3), &
            0.22_dp*midpoint(1) + 0.05_dp*midpoint(3), &
            0.16_dp*midpoint(2) - 0.08_dp*midpoint(1)]
        edge_target(edge) = dot_product(potential, tangent)
    end do
    rhs = matmul(gram, edge_target)
    discrete_gradient = matmul(transpose(real(incidence, dp)), vertex_potential)
    call build_tree_cotree_gauge(incidence, gauge, status)
    call check_condition(status == 0, "3-D cube tree-cotree gauge builds")
    if (status /= 0) then
        call check_summary("Biro 3-D tree-cotree gallery")
        return
    end if
    call tree_cotree_gauge_edges(gauge, tree, cotree, status)
    call check_condition(status == 0 .and. size(tree) == 7 .and. &
        size(cotree) == 5, "cube has seven tree and five cotree edges")
    call check_condition(maxval(abs(matmul(gram, discrete_gradient))) < &
        3.0e-14_dp, "face curl Gram has its expected five-dimensional kernel")
    call reduce_tree_cotree_dense_system(gauge, gram, rhs, &
        reduced_matrix, reduced_rhs, status)
    call check_condition(status == 0, "3-D reduced block assembles")
    if (status == 0) then
        call check_condition(minval([(reduced_matrix(edge, edge), &
            edge=1, size(reduced_matrix, 1))]) > 0.0_dp, &
            "cotree block has positive diagonal energy")
        call solve_reference(reduced_matrix, reduced_rhs, reduced, status)
        call check_condition(status == 0, "independent reference solve succeeds")
        call apply_tree_cotree_prolongation(gauge, reduced, prolonged, status)
        call check_condition(status == 0, "3-D cotree prolongation succeeds")
        if (status == 0) then
            call check_condition(maxval(abs(matmul(gram, prolonged) - rhs)) < &
                3.0e-11_dp, "prolonged edge solution satisfies manufactured curl-curl")
            call check_condition(maxval(abs(prolonged(tree))) < 1.0e-13_dp, &
                "tree edge potential is fixed by the direct gauge")
        end if
    end if
    call check_summary("Biro 3-D tree-cotree gallery")

contains

    subroutine solve_reference(matrix, rhs, solution, status)
        real(dp), intent(in) :: matrix(:, :), rhs(:)
        real(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status
        real(dp), allocatable :: work(:, :), load(:)
        real(dp) :: factor
        integer :: n, row, column

        n = size(rhs)
        status = 1
        allocate(solution(0))
        if (size(matrix, 1) /= n .or. size(matrix, 2) /= n) return
        allocate(work(n, n), load(n))
        work = matrix
        load = rhs
        do column = 1, n - 1
            if (abs(work(column, column)) < 1.0e-13_dp) return
            do row = column + 1, n
                factor = work(row, column)/work(column, column)
                work(row, column:n) = work(row, column:n) - &
                    factor*work(column, column:n)
                load(row) = load(row) - factor*load(column)
            end do
        end do
        if (abs(work(n, n)) < 1.0e-13_dp) return
        if (allocated(solution)) deallocate(solution)
        allocate(solution(n))
        solution(n) = load(n)/work(n, n)
        do row = n - 1, 1, -1
            solution(row) = (load(row) - &
                dot_product(work(row, row + 1:n), solution(row + 1:n)))/ &
                work(row, row)
        end do
        status = 0
    end subroutine solve_reference

end program test_biro_tree_cotree_3d_gallery
