module assembly_2d_module
    use fortfem_kinds
    use basis_p1_2d_module
    use fortfem_mesh_2d
    use fortfem_sparse_matrix, only: triplet_matrix_t
    implicit none
    private
    
    public :: assembly_2d_t, assemble_laplacian, assemble_mass_matrix, assemble_mass_rhs
    
    type :: assembly_2d_t
        integer :: n_dofs = 0
        integer :: max_nnz = 0
        type(basis_p1_2d_t) :: basis
    contains
        procedure :: init
        procedure :: destroy
        procedure :: element_mass_matrix
        procedure :: element_stiffness_matrix
        procedure :: element_load_vector
        procedure :: assemble_global
    end type assembly_2d_t
    
    ! 3-point edge midpoint rule for triangles (exact for degree 2 polynomials)
    integer, parameter :: n_gauss = 3
    real(dp), parameter :: gauss_points(2,3) = reshape([ &
        0.5_dp, 0.0_dp, &     ! Edge midpoint (1/2, 0)
        0.5_dp, 0.5_dp, &     ! Edge midpoint (1/2, 1/2)  
        0.0_dp, 0.5_dp  &     ! Edge midpoint (0, 1/2)
    ], [2, 3])
    real(dp), parameter :: gauss_weights(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
    
contains

    subroutine init(this, n_dofs, max_nnz)
        class(assembly_2d_t), intent(inout) :: this
        integer, intent(in) :: n_dofs, max_nnz
        
        this%n_dofs = n_dofs
        this%max_nnz = max_nnz
        
    end subroutine init
    
    subroutine destroy(this)
        class(assembly_2d_t), intent(inout) :: this
        
        this%n_dofs = 0
        this%max_nnz = 0
        
    end subroutine destroy
    
    pure subroutine element_mass_matrix(this, vertices, mass)
        class(assembly_2d_t), intent(in) :: this
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: mass(3,3)
        
        real(dp) :: jac(2,2), det_j
        real(dp) :: phi_i, phi_j
        integer :: i, j, q
        
        ! Get Jacobian (constant for linear elements)
        call this%basis%compute_jacobian(vertices, jac, det_j)
        
        ! Initialize mass matrix
        mass = 0.0_dp
        
        ! Quadrature loop
        do q = 1, n_gauss
            ! Basis function values at quadrature point
            do i = 1, 3
                phi_i = this%basis%eval(i, gauss_points(1,q), gauss_points(2,q))
                do j = 1, 3
                    phi_j = this%basis%eval(j, gauss_points(1,q), gauss_points(2,q))
                    mass(i,j) = mass(i,j) + gauss_weights(q) * phi_i * phi_j * abs(det_j)
                end do
            end do
        end do
        
    end subroutine element_mass_matrix
    
    pure subroutine element_stiffness_matrix(this, vertices, stiff)
        class(assembly_2d_t), intent(in) :: this
        real(dp), intent(in) :: vertices(2,3)
        real(dp), intent(out) :: stiff(3,3)
        
        real(dp) :: jac(2,2), det_j, jac_inv(2,2)
        real(dp) :: grad_ref(2), grad_phys(2)
        real(dp) :: grad_i(2), grad_j(2)
        integer :: i, j, q
        
        ! Get Jacobian and its inverse
        call this%basis%compute_jacobian(vertices, jac, det_j)
        
        ! Compute inverse Jacobian
        jac_inv(1,1) = jac(2,2) / det_j
        jac_inv(1,2) = -jac(1,2) / det_j
        jac_inv(2,1) = -jac(2,1) / det_j
        jac_inv(2,2) = jac(1,1) / det_j
        
        ! Initialize stiffness matrix
        stiff = 0.0_dp
        
        ! For P1 elements, gradients are constant, so we only need one quadrature point
        do i = 1, 3
            ! Get gradient in reference coordinates
            grad_ref = this%basis%grad(i, 0.0_dp, 0.0_dp)
            
            ! Transform to physical coordinates: grad_phys = J^{-T} * grad_ref
            grad_i(1) = jac_inv(1,1) * grad_ref(1) + jac_inv(2,1) * grad_ref(2)
            grad_i(2) = jac_inv(1,2) * grad_ref(1) + jac_inv(2,2) * grad_ref(2)
            
            do j = 1, 3
                grad_ref = this%basis%grad(j, 0.0_dp, 0.0_dp)
                grad_j(1) = jac_inv(1,1) * grad_ref(1) + jac_inv(2,1) * grad_ref(2)
                grad_j(2) = jac_inv(1,2) * grad_ref(1) + jac_inv(2,2) * grad_ref(2)
                
                ! Stiffness entry: integral of grad(phi_i) . grad(phi_j)
                stiff(i,j) = abs(det_j) * (grad_i(1) * grad_j(1) + grad_i(2) * grad_j(2))
            end do
        end do
        
    end subroutine element_stiffness_matrix
    
    pure subroutine element_load_vector(this, vertices, source_func, load)
        class(assembly_2d_t), intent(in) :: this
        real(dp), intent(in) :: vertices(2,3)
        interface
            pure function source_func(x, y) result(val)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp) :: val
            end function source_func
        end interface
        real(dp), intent(out) :: load(3)
        
        real(dp) :: jac(2,2), det_j
        real(dp) :: x, y, f_val, phi_i
        integer :: i, q
        
        ! Get Jacobian
        call this%basis%compute_jacobian(vertices, jac, det_j)
        
        ! Initialize load vector
        load = 0.0_dp
        
        ! Quadrature loop
        do q = 1, n_gauss
            ! Map quadrature point to physical coordinates
            call this%basis%transform_to_physical(gauss_points(1,q), gauss_points(2,q), &
                                                 vertices, x, y)
            
            ! Evaluate source function
            f_val = source_func(x, y)
            
            ! Add contribution from each basis function
            do i = 1, 3
                phi_i = this%basis%eval(i, gauss_points(1,q), gauss_points(2,q))
                load(i) = load(i) + 2.0_dp * gauss_weights(q) * f_val * phi_i * abs(det_j)
            end do
        end do
        
    end subroutine element_load_vector
    
    subroutine assemble_global(this, mesh, source_func, matrix, rhs)
        class(assembly_2d_t), intent(in) :: this
        type(mesh_2d_t), intent(in) :: mesh
        interface
            pure function source_func(x, y) result(val)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp) :: val
            end function source_func
        end interface
        type(triplet_matrix_t), intent(inout) :: matrix
        real(dp), intent(out) :: rhs(:)
        
        real(dp) :: vertices(2,3), elem_stiff(3,3), elem_load(3)
        integer :: elem, i, j, gi, gj
        
        ! Initialize
        call matrix%init(mesh%n_vertices, this%max_nnz)
        rhs = 0.0_dp
        
        ! Loop over elements
        do elem = 1, mesh%n_triangles
            ! Get element vertices
            do i = 1, 3
                vertices(:,i) = mesh%vertices(:, mesh%triangles(i,elem))
            end do
            
            ! Compute element matrices
            call this%element_stiffness_matrix(vertices, elem_stiff)
            call this%element_load_vector(vertices, source_func, elem_load)
            
            ! Assemble into global system
            do i = 1, 3
                gi = mesh%triangles(i, elem)
                
                ! Add to RHS
                rhs(gi) = rhs(gi) + elem_load(i)
                
                ! Add to matrix
                do j = 1, 3
                    gj = mesh%triangles(j, elem)
                    call matrix%add(gi, gj, elem_stiff(i,j))
                end do
            end do
        end do
        
    end subroutine assemble_global

    ! Standalone assembly routines for backward compatibility
    subroutine assemble_laplacian(space, matrix)
        use function_space_module
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: matrix(:,:)
        
        integer :: i, j, n
        
        n = space%n_dofs
        matrix = 0.0_dp
        
        ! Simple diagonal dominant matrix for now
        do i = 1, n
            matrix(i, i) = 4.0_dp
            if (i > 1) matrix(i, i-1) = -1.0_dp
            if (i < n) matrix(i, i+1) = -1.0_dp
        end do
    end subroutine assemble_laplacian

    subroutine assemble_mass_matrix(space, matrix)
        use function_space_module
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: matrix(:,:)
        
        integer :: i, n
        
        n = space%n_dofs
        matrix = 0.0_dp
        
        ! Simple mass matrix (identity for now)
        do i = 1, n
            matrix(i, i) = 1.0_dp
        end do
    end subroutine assemble_mass_matrix

    subroutine assemble_mass_rhs(space, rhs, source_func)
        use function_space_module
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: rhs(:)
        interface
            function source_func(x, y) result(f)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp) :: f
            end function source_func
        end interface
        
        integer :: i, n
        
        n = space%n_dofs
        rhs = 0.0_dp
        
        ! Simple unit RHS vector
        do i = 1, n
            rhs(i) = source_func(real(i, dp)/real(n, dp), 0.5_dp)
        end do
    end subroutine assemble_mass_rhs

end module assembly_2d_module