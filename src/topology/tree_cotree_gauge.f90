module fortfem_tree_cotree_gauge
    !! Fixed-topology tree-cotree reduction metadata.
    !!
    !! A curl-conforming edge potential has a discrete gradient kernel. A
    !! spanning forest of the vertex-edge incidence graph identifies tree
    !! edge degrees of freedom that can be fixed to zero; the remaining
    !! cotree degrees of freedom form the direct-solve system. This module
    !! owns only graph topology and restriction/prolongation algebra. Metric
    !! matrices, high-order basis maps, and physical boundary laws remain in
    !! the consuming H(curl) layer.
    use fortfem_kinds, only: dp
    implicit none
    private

    type, public :: tree_cotree_gauge_t
        private
        integer :: vertex_count = 0
        integer :: edge_count = 0
        integer, allocatable :: tree_edges(:)
        integer, allocatable :: cotree_edges(:)
    end type tree_cotree_gauge_t

    public :: build_tree_cotree_gauge
    public :: validate_tree_cotree_gauge
    public :: tree_cotree_gauge_edges
    public :: apply_tree_cotree_restriction
    public :: apply_tree_cotree_prolongation
    public :: reduce_tree_cotree_dense_system
    public :: reduce_tree_cotree_dense_system_jvp
    public :: reduce_tree_cotree_dense_system_vjp

contains

    subroutine build_tree_cotree_gauge(incidence, gauge, status)
        !! Build a spanning-forest tree/cotree split from an oriented graph.
        integer, intent(in) :: incidence(:, :)
        type(tree_cotree_gauge_t), intent(inout) :: gauge
        integer, intent(out) :: status
        integer, allocatable :: parent(:), tree_buffer(:), cotree_buffer(:)
        integer :: vertex_count, edge_count, edge, vertex
        integer :: plus_vertex, minus_vertex, root_plus, root_minus
        integer :: tree_count, cotree_count

        call clear_tree_cotree_gauge(gauge)
        status = 1
        vertex_count = size(incidence, 1)
        edge_count = size(incidence, 2)
        if (vertex_count < 1) return
        if (edge_count < 0) return
        allocate(parent(vertex_count))
        do vertex = 1, vertex_count
            parent(vertex) = vertex
        end do
        allocate(tree_buffer(edge_count), cotree_buffer(edge_count))
        tree_count = 0
        cotree_count = 0
        do edge = 1, edge_count
            call incidence_endpoints( &
                incidence(:, edge), plus_vertex, minus_vertex, status)
            if (status /= 0) then
                call clear_tree_cotree_gauge(gauge)
                return
            end if
            root_plus = find_root(parent, plus_vertex)
            root_minus = find_root(parent, minus_vertex)
            if (root_plus /= root_minus) then
                tree_count = tree_count + 1
                tree_buffer(tree_count) = edge
                parent(root_minus) = root_plus
            else
                cotree_count = cotree_count + 1
                cotree_buffer(cotree_count) = edge
            end if
        end do

        gauge%vertex_count = vertex_count
        gauge%edge_count = edge_count
        allocate(gauge%tree_edges(tree_count))
        allocate(gauge%cotree_edges(cotree_count))
        if (tree_count > 0) gauge%tree_edges = tree_buffer(:tree_count)
        if (cotree_count > 0) then
            gauge%cotree_edges = cotree_buffer(:cotree_count)
        end if
        call validate_tree_cotree_gauge(gauge, status)
    end subroutine build_tree_cotree_gauge

    subroutine validate_tree_cotree_gauge(gauge, status)
        type(tree_cotree_gauge_t), intent(in) :: gauge
        integer, intent(out) :: status
        integer :: edge

        status = 1
        if (gauge%vertex_count < 1 .or. gauge%edge_count < 0) return
        if (.not. allocated(gauge%tree_edges) .or. &
            .not. allocated(gauge%cotree_edges)) return
        if (size(gauge%tree_edges) + size(gauge%cotree_edges) /= &
            gauge%edge_count) then
            status = 2
            return
        end if
        if (any(gauge%tree_edges < 1) .or. &
            any(gauge%tree_edges > gauge%edge_count) .or. &
            any(gauge%cotree_edges < 1) .or. &
            any(gauge%cotree_edges > gauge%edge_count)) then
            status = 3
            return
        end if
        do edge = 1, gauge%edge_count
            if (count(gauge%tree_edges == edge) + &
                count(gauge%cotree_edges == edge) /= 1) then
                status = 4
                return
            end if
        end do
        if (size(gauge%tree_edges) > gauge%vertex_count - 1) then
            status = 5
            return
        end if
        status = 0
    end subroutine validate_tree_cotree_gauge

    subroutine tree_cotree_gauge_edges( &
            gauge, tree_edges, cotree_edges, status)
        type(tree_cotree_gauge_t), intent(in) :: gauge
        integer, allocatable, intent(out) :: tree_edges(:), cotree_edges(:)
        integer, intent(out) :: status

        if (allocated(tree_edges)) deallocate(tree_edges)
        if (allocated(cotree_edges)) deallocate(cotree_edges)
        call validate_tree_cotree_gauge(gauge, status)
        if (status /= 0) then
            allocate(tree_edges(0), cotree_edges(0))
            return
        end if
        allocate(tree_edges(size(gauge%tree_edges)))
        allocate(cotree_edges(size(gauge%cotree_edges)))
        if (size(tree_edges) > 0) tree_edges = gauge%tree_edges
        if (size(cotree_edges) > 0) cotree_edges = gauge%cotree_edges
        status = 0
    end subroutine tree_cotree_gauge_edges

    subroutine apply_tree_cotree_restriction( &
            gauge, full_vector, reduced_vector, status)
        type(tree_cotree_gauge_t), intent(in) :: gauge
        real(dp), intent(in) :: full_vector(:)
        real(dp), allocatable, intent(out) :: reduced_vector(:)
        integer, intent(out) :: status
        integer :: reduced

        if (allocated(reduced_vector)) deallocate(reduced_vector)
        call validate_tree_cotree_gauge(gauge, status)
        if (status /= 0 .or. size(full_vector) /= gauge%edge_count) then
            status = 2
            allocate(reduced_vector(0))
            return
        end if
        allocate(reduced_vector(size(gauge%cotree_edges)))
        do reduced = 1, size(reduced_vector)
            reduced_vector(reduced) = full_vector(gauge%cotree_edges(reduced))
        end do
        status = 0
    end subroutine apply_tree_cotree_restriction

    subroutine apply_tree_cotree_prolongation( &
            gauge, reduced_vector, full_vector, status)
        type(tree_cotree_gauge_t), intent(in) :: gauge
        real(dp), intent(in) :: reduced_vector(:)
        real(dp), allocatable, intent(out) :: full_vector(:)
        integer, intent(out) :: status
        integer :: reduced

        if (allocated(full_vector)) deallocate(full_vector)
        call validate_tree_cotree_gauge(gauge, status)
        if (status /= 0 .or. size(reduced_vector) /= &
            size(gauge%cotree_edges)) then
            status = 2
            allocate(full_vector(0))
            return
        end if
        allocate(full_vector(gauge%edge_count))
        full_vector = 0.0_dp
        do reduced = 1, size(reduced_vector)
            full_vector(gauge%cotree_edges(reduced)) = reduced_vector(reduced)
        end do
        status = 0
    end subroutine apply_tree_cotree_prolongation

    subroutine reduce_tree_cotree_dense_system( &
            gauge, matrix, right_hand_side, reduced_matrix, reduced_rhs, status)
        !! Extract R^T matrix R and R^T right_hand_side for fixed R.
        type(tree_cotree_gauge_t), intent(in) :: gauge
        real(dp), intent(in) :: matrix(:, :), right_hand_side(:)
        real(dp), allocatable, intent(out) :: reduced_matrix(:, :)
        real(dp), allocatable, intent(out) :: reduced_rhs(:)
        integer, intent(out) :: status
        integer :: row, column

        if (allocated(reduced_matrix)) deallocate(reduced_matrix)
        if (allocated(reduced_rhs)) deallocate(reduced_rhs)
        call validate_tree_cotree_gauge(gauge, status)
        if (status /= 0 .or. size(matrix, 1) /= gauge%edge_count .or. &
            size(matrix, 2) /= gauge%edge_count .or. &
            size(right_hand_side) /= gauge%edge_count) then
            status = 2
            allocate(reduced_matrix(0, 0), reduced_rhs(0))
            return
        end if
        allocate(reduced_matrix( &
            size(gauge%cotree_edges), size(gauge%cotree_edges)))
        allocate(reduced_rhs(size(gauge%cotree_edges)))
        do row = 1, size(gauge%cotree_edges)
            reduced_rhs(row) = right_hand_side(gauge%cotree_edges(row))
            do column = 1, size(gauge%cotree_edges)
                reduced_matrix(row, column) = matrix( &
                    gauge%cotree_edges(row), gauge%cotree_edges(column))
            end do
        end do
        status = 0
    end subroutine reduce_tree_cotree_dense_system

    subroutine reduce_tree_cotree_dense_system_jvp( &
            gauge, matrix_dot, right_hand_side_dot, reduced_matrix_dot, &
            reduced_rhs_dot, status)
        !! Apply the fixed-topology JVP of the reduced dense system.
        type(tree_cotree_gauge_t), intent(in) :: gauge
        real(dp), intent(in) :: matrix_dot(:, :), right_hand_side_dot(:)
        real(dp), allocatable, intent(out) :: reduced_matrix_dot(:, :)
        real(dp), allocatable, intent(out) :: reduced_rhs_dot(:)
        integer, intent(out) :: status

        call reduce_tree_cotree_dense_system( &
            gauge, matrix_dot, right_hand_side_dot, reduced_matrix_dot, &
            reduced_rhs_dot, status)
    end subroutine reduce_tree_cotree_dense_system_jvp

    subroutine reduce_tree_cotree_dense_system_vjp( &
            gauge, reduced_matrix_bar, reduced_rhs_bar, matrix_bar, rhs_bar, &
            status)
        !! Apply the real VJP of the fixed-topology reduced system.
        type(tree_cotree_gauge_t), intent(in) :: gauge
        real(dp), intent(in) :: reduced_matrix_bar(:, :), reduced_rhs_bar(:)
        real(dp), allocatable, intent(out) :: matrix_bar(:, :), rhs_bar(:)
        integer, intent(out) :: status
        integer :: row, column

        if (allocated(matrix_bar)) deallocate(matrix_bar)
        if (allocated(rhs_bar)) deallocate(rhs_bar)
        call validate_tree_cotree_gauge(gauge, status)
        if (status /= 0 .or. &
            size(reduced_matrix_bar, 1) /= size(gauge%cotree_edges) .or. &
            size(reduced_matrix_bar, 2) /= size(gauge%cotree_edges) .or. &
            size(reduced_rhs_bar) /= size(gauge%cotree_edges)) then
            status = 2
            allocate(matrix_bar(0, 0), rhs_bar(0))
            return
        end if
        allocate(matrix_bar(gauge%edge_count, gauge%edge_count))
        allocate(rhs_bar(gauge%edge_count))
        matrix_bar = 0.0_dp
        rhs_bar = 0.0_dp
        do row = 1, size(gauge%cotree_edges)
            rhs_bar(gauge%cotree_edges(row)) = reduced_rhs_bar(row)
            do column = 1, size(gauge%cotree_edges)
                matrix_bar(gauge%cotree_edges(row), &
                    gauge%cotree_edges(column)) = reduced_matrix_bar(row, column)
            end do
        end do
        status = 0
    end subroutine reduce_tree_cotree_dense_system_vjp

    subroutine incidence_endpoints(column, plus_vertex, minus_vertex, status)
        integer, intent(in) :: column(:)
        integer, intent(out) :: plus_vertex, minus_vertex
        integer, intent(out) :: status
        integer :: vertex

        plus_vertex = 0
        minus_vertex = 0
        status = 1
        do vertex = 1, size(column)
            if (column(vertex) == 1) then
                if (plus_vertex /= 0) return
                plus_vertex = vertex
            else if (column(vertex) == -1) then
                if (minus_vertex /= 0) return
                minus_vertex = vertex
            else if (column(vertex) /= 0) then
                return
            end if
        end do
        if (plus_vertex == 0 .or. minus_vertex == 0) return
        status = 0
    end subroutine incidence_endpoints

    integer function find_root(parent, node) result(root)
        integer, intent(inout) :: parent(:)
        integer, intent(in) :: node
        integer :: current, next

        root = node
        do while (parent(root) /= root)
            root = parent(root)
        end do
        current = node
        do while (parent(current) /= current)
            next = parent(current)
            parent(current) = root
            current = next
        end do
    end function find_root

    subroutine clear_tree_cotree_gauge(gauge)
        type(tree_cotree_gauge_t), intent(inout) :: gauge

        if (allocated(gauge%tree_edges)) deallocate(gauge%tree_edges)
        if (allocated(gauge%cotree_edges)) deallocate(gauge%cotree_edges)
        gauge%vertex_count = 0
        gauge%edge_count = 0
    end subroutine clear_tree_cotree_gauge

end module fortfem_tree_cotree_gauge
