program test_lab_unstructured_mesh
    ! Curved unstructured quad/hex input (Gmsh MSH 4.1 fixtures, island-chain
    ! generator, quad x phi extrusion) on the matrix-free gather-map operator.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh, &
                                   build_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                    apply_tensor_operator, assemble_tensor_operator, tensor_load, &
                                       tensor_load_flux
    use fortfem_tensor_mesh_tools, only: tensor_mesh_quality_t, tensor_mesh_check, &
                                         orient_tensor_cells, extrude_tensor_mesh, &
                                         mark_boundary_faces
    use fortfem_tensor_mesh_gmsh, only: read_gmsh_tensor_mesh
    use fortfem_island_quad_mesh, only: island_chain_t, island_mesh_counts_t, &
                                        island_chain_quad_mesh
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), delta = 2.e-3_dp
    character(len=*), parameter :: data = 'test/data/gmsh/'
    type(tensor_mesh_t) :: mesh, other
    type(tensor_mesh_quality_t) :: quality
    type(island_chain_t) :: chain
    type(island_mesh_counts_t) :: counts
    character(len=:), allocatable :: message
    real(dp) :: area(3), measure
    real(dp), allocatable :: geom(:, :, :)
    integer, allocatable :: corners(:, :)
    integer :: order, p, status, flipped, k, nrow, ncell, n
    ! Gmsh quadrilaterals of geometry order 1..3 on the unit disk.
    do order = 1, 3
        do p = 1, 6
            call read_gmsh_tensor_mesh(data//'disk_quad_o'//char(48 + order)// &
                                       '.msh', p, mesh, status, message=message)
            if (status /= 0) then
                print *, message
                error stop 'read disk'
            end if
            if (mesh%geometry_degree /= order) error stop 'disk order'
            call check_mesh(mesh, 'disk', measure)
            if (p == order) area(order) = measure
        end do
    end do
    print '(a,3es10.2)', 'disk area errors o1..o3:', abs(area - pi)
    ! Optimised boundary nodes stay on the circle but are not equispaced, so
    ! only bounds are asserted; the energy oracle above is exact.
    if (abs(area(1) - pi) > 0.1_dp .or. abs(area(2) - pi) > 2.e-3_dp .or. &
        abs(area(3) - pi) > 2.e-3_dp) error stop 'isoparametric disk area'
    ! Gmsh's own minSJ is negative in three cells of this unoptimised order-3
    ! disk (make_fixtures.py): the check must report exactly those.
    call read_gmsh_tensor_mesh(data//'disk_quad_o3_inverted.msh', 3, mesh, status)
    if (status /= 0) error stop 'read inverted disk'
    call tensor_mesh_check(mesh, quality)
    if (quality%inverted_cells /= 3) error stop 'Gmsh-confirmed inverted cells'
    ! Mirrored cells are re-oriented, not rejected.
    call read_gmsh_tensor_mesh(data//'disk_quad_o2.msh', 2, mesh, status)
    geom = mesh%geometry
    corners = mesh%vertices
    geom(1, :, :) = -geom(1, :, :)
    call orient_tensor_cells(2, 2, corners, geom, flipped)
    if (flipped /= mesh%nelem) error stop 'all mirrored cells flipped'
    call build_tensor_mesh(2, 2, 2, corners, geom, other, status)
    if (status /= 0) error stop 'mirrored mesh'
    call mark_boundary_faces(other)
    call check_mesh(other, 'mirrored', measure)
    if (abs(measure - area(2)) > 1.e-13_dp) error stop 'mirrored area'
    ! Curved extruded hexahedra (hex27) and a z-periodic curved box (hex27/64).
    do p = 1, 4
        call read_gmsh_tensor_mesh(data//'cylinder_hex_o2.msh', p, mesh, status)
        if (status /= 0) error stop 'read cylinder'
        call check_mesh(mesh, 'cylinder', measure)
        if (abs(measure - pi/2) > 2.e-2_dp) error stop 'cylinder volume'
    end do
    do order = 2, 3
        do p = 1, 4
            call read_gmsh_tensor_mesh(data//'box_periodic_hex_o'// &
                                       char(48 + order)//'.msh', p, mesh, status)
            if (status /= 0) error stop 'read periodic box'
            ! Two cells in x and y, three periodic layers in z.
            if (mesh%ndof /= (2*p + 1)**2*3*p) error stop 'periodic box dofs'
            if (count(mesh%boundary) /= 24*p*p) error stop 'periodic box boundary'
            call check_mesh(mesh, 'box', measure, periodic=.true.)
            if (abs(measure - 1) > 1.e-3_dp) error stop 'box volume'
        end do
    end do
    ! Island-chain generator for the level function -(r-1/2)^2 + delta cos(4u):
    ! O-points at u = k pi/2 and X-points at u = pi/4 + k pi/2, radius 1/2.
    allocate (chain%opoint(2, 4), chain%xpoint(2, 4))
    do k = 1, 4
        chain%opoint(:, k) = 0.5_dp*[cos((k - 1)*pi/2), sin((k - 1)*pi/2)]
        chain%xpoint(:, k) = 0.5_dp*[cos((k - 0.5_dp)*pi/2), &
                                     sin((k - 0.5_dp)*pi/2)]
    end do
    chain%r_core = 0.15_dp
    chain%r_inner = 0.3_dp
    chain%r_outer = 0.7_dp
    chain%r_wall = 1
    counts%half = 2
    do p = 1, 6
        call island_chain_quad_mesh(level, chain, counts, p, 4, mesh, status)
        if (status /= 0) error stop 'island mesh'
        call check_mesh(mesh, 'island', measure)
        if (abs(measure - pi) > 1.e-6_dp) error stop 'island disk area'
    end do
    n = 2*4*counts%half
    nrow = counts%core_rings + counts%inner_rows + 1 + counts%outer_rows + &
           counts%wall_rows
    ncell = (n/4)**2 + (nrow - 1)*n + 4*(4*counts%half*(counts%island_rings + 1) + &
                                         counts%half**2)
    if (mesh%nelem /= ncell) error stop 'island cell count'
    call check_island_topology(mesh)
    call extrude_tensor_mesh(mesh, 3, 2, other, status, map=torus, periodic=.true.)
    if (status /= 0) error stop 'extrusion'
    ! Wall of n edges at p = 2 (2n nodes) times 3 periodic layers (6 planes).
    if (count(other%boundary) /= 2*n*3*2) error stop 'periodic extrusion boundary'
    call check_mesh(other, 'island x phi', measure, periodic=.true.)
    if (abs(measure/(6*pi**2) - 1) > 1.e-3_dp) error stop 'torus volume'
    call check_flux_load()
    print '(a)', 'unstructured mesh tests passed'
contains
    function level(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value, r
        r = norm2(x)
        value = -(r - 0.5_dp)**2 + delta*(x(1)**4 - 6*x(1)**2*x(2)**2 + x(2)**4)/r**4
    end function

    subroutine torus(xi, x)
        ! Quad x phi: poloidal plane at major radius 3, toroidal angle 2 pi s.
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        x(1) = (3 + xi(1))*cos(2*pi*xi(3))
        x(2) = (3 + xi(1))*sin(2*pi*xi(3))
        x(3) = xi(2)
    end subroutine

    subroutine check_mesh(mesh, name, measure, periodic)
        type(tensor_mesh_t), intent(in) :: mesh
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: measure
        logical, intent(in), optional :: periodic
        type(tensor_operator_t) :: op, plain
        type(csr_matrix_t) :: a
        type(tensor_mesh_quality_t) :: quality
        real(dp), allocatable :: x(:), y(:), z(:), one(:)
        real(dp) :: err, slope(3), energy
        integer :: i, j
        call tensor_mesh_check(mesh, quality)
        if (quality%inverted_cells /= 0 .or. quality%min_jacobian <= 0) then
            print *, name, quality%min_jacobian
            error stop 'positive Jacobian'
        end if
        if (.not. present(periodic)) then
            if (quality%max_gap > 1.e-13_dp) then
                print *, name, quality%max_gap
                error stop 'conforming shared nodes'
            end if
        end if
        call build_tensor_operator(mesh, coefficient, op, status)
        if (status /= 0) error stop 'operator'
        allocate (x(mesh%ndof), y(mesh%ndof), z(mesh%ndof))
        do i = 1, mesh%ndof
            x(i) = sin(real(13*i, dp))
        end do
        call apply_tensor_operator(op, x, y)
        call assemble_tensor_operator(op, a)
        z = 0
        do i = 1, a%n
            do j = a%rowptr(i), a%rowptr(i + 1) - 1
                z(i) = z(i) + a%val(j)*x(a%colidx(j))
            end do
        end do
        err = norm2(y - z)/max(1.0_dp, norm2(y))
        if (err > 1.e-12_dp) then
            print '(a,a,i2,es10.2)', name, ' p=', mesh%degree, err
            error stop 'MF versus assembled'
        end if
        ! Isoparametric reproduction of a linear field: energy = |a|^2 * measure.
        ! Enough points to integrate the geometry Jacobian exactly.
        call build_tensor_operator(mesh, identity, plain, status, &
                                   nq=max(mesh%degree, mesh%geometry_degree) + 2)
        if (status /= 0) error stop 'unit operator'
        measure = sum(plain%volume)
        if (mesh%degree >= mesh%geometry_degree) then
            slope = [0.7_dp, -0.4_dp, 0.0_dp]
            if (mesh%dim == 3 .and. .not. present(periodic)) slope(3) = 0.3_dp
            do i = 1, mesh%ndof
                x(i) = dot_product(slope(1:mesh%dim), mesh%coordinates(:, i))
            end do
            if (present(periodic)) then
                ! Only fields that are single-valued across the period.
                if (mesh%dim == 3 .and. name == 'island x phi') then
                    do i = 1, mesh%ndof
                        x(i) = 0.4_dp*mesh%coordinates(3, i)
                    end do
                    slope = [0.0_dp, 0.0_dp, 0.4_dp]
                end if
            end if
            call apply_tensor_operator(plain, x, y, eliminate=.false.)
            energy = dot_product(x, y)
            if (abs(energy - sum(slope**2)*measure) > &
                1.e-12_dp*sum(slope**2)*measure) then
                print *, name, mesh%degree, energy, sum(slope**2)*measure
                error stop 'linear energy oracle'
            end if
        end if
        allocate (one(mesh%ndof))
        one = 1
        call apply_tensor_operator(plain, one, y, eliminate=.false.)
        if (maxval(abs(y)) > 1.e-12_dp) error stop 'constants in the kernel'
    end subroutine

    subroutine check_island_topology(mesh)
        type(tensor_mesh_t), intent(in) :: mesh
        integer, allocatable :: valence(:)
        integer :: e, c, k, v, nv, nedge, i, j, pairs(2, 4)
        integer, allocatable :: keys(:, :)
        real(dp) :: corner(2), sep
        logical :: found
        logical, allocatable :: onsep(:)
        nv = maxval(mesh%vertices)
        allocate (valence(nv), keys(2, 4*mesh%nelem))
        valence = 0
        nedge = 0
        pairs = reshape([1, 2, 2, 4, 4, 3, 3, 1], [2, 4])
        do e = 1, mesh%nelem
            do c = 1, 4
                valence(mesh%vertices(c, e)) = valence(mesh%vertices(c, e)) + 1
                i = minval(mesh%vertices(pairs(:, c), e))
                j = maxval(mesh%vertices(pairs(:, c), e))
                found = .false.
                do k = 1, nedge
                    if (keys(1, k) == i .and. keys(2, k) == j) found = .true.
                end do
                if (.not. found) then
                    nedge = nedge + 1
                    keys(:, nedge) = [i, j]
                end if
            end do
        end do
        if (nv - nedge + mesh%nelem /= 1) error stop 'disk Euler characteristic'
        ! Each X-point is one vertex; its four separatrix sectors (inner band,
        ! outer band, two islands) are each split by one edge: eight cells.
        sep = level(chain%xpoint(:, 1))
        do k = 1, 4
            found = .false.
            do e = 1, mesh%nelem
                do c = 1, 4
                    corner = mesh%geometry(:, 1 + merge(mesh%geometry_degree, 0, &
                                     modulo(c - 1, 2) == 1) + (mesh%geometry_degree + 1)* &
                                     merge(mesh%geometry_degree, 0, c > 2), e)
                    if (norm2(corner - chain%xpoint(:, k)) > 1.e-14_dp) cycle
                    v = mesh%vertices(c, e)
                    found = .true.
                end do
            end do
            if (.not. found) error stop 'X-point is a vertex'
            if (valence(v) /= 8) error stop 'X-point valence eight'
        end do
        ! Separatrix alignment: both branches, 2n - m distinct vertices, lie on
        ! the level set through the X-points.
        allocate (onsep(nv))
        onsep = .false.
        do e = 1, mesh%nelem
            do c = 1, 4
                corner = mesh%geometry(:, 1 + merge(mesh%geometry_degree, 0, &
                                 modulo(c - 1, 2) == 1) + (mesh%geometry_degree + 1)* &
                                 merge(mesh%geometry_degree, 0, c > 2), e)
                if (abs(level(corner) - sep) < 1.e-12_dp) &
                    onsep(mesh%vertices(c, e)) = .true.
            end do
        end do
        if (count(onsep) /= 2*(8*counts%half) - 4) then
            print *, count(onsep)
            error stop 'separatrix vertices'
        end if
    end subroutine

    subroutine check_flux_load()
        ! Weak flux load equals the strong load for a quadratic field on affine
        ! cells, where both quadratures are exact.
        type(tensor_mesh_t) :: box
        type(tensor_operator_t) :: op
        real(dp), allocatable :: strong(:), weak(:)
        integer :: dim, p
        do dim = 2, 3
            do p = 2, 4
                call structured_tensor_mesh(dim, [3, 2, 2], p, box, status)
                call build_tensor_operator(box, identity, op, status)
                call tensor_load(op, laplacian, strong)
                call tensor_load_flux(op, gradient, weak)
                if (norm2(strong - weak) > 1.e-12_dp*norm2(strong)) &
                    error stop 'flux load'
            end do
        end do
    end subroutine

    function laplacian(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = -(4 - 2 + merge(6.0_dp, 0.0_dp, size(x) == 3)) + 0*x(1)
    end function

    subroutine gradient(x, f)
        ! T = 2x^2 + 3xy - y^2 (+ 3 z^2): grad T.
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        f(1) = 4*x(1) + 3*x(2)
        f(2) = 3*x(1) - 2*x(2)
        if (size(x) == 3) f(3) = 6*x(3)
    end subroutine

    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1 + 0.1_dp*x(i)**2
        end do
        k(1, 2) = 0.17_dp
        k(2, 1) = 0.17_dp
        c = 0.3_dp
    end subroutine

    subroutine identity(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1
        end do
        c = 0
    end subroutine
end program
