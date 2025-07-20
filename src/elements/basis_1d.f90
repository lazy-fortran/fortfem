module fortfem_basis_1d
    use fortfem_kinds
    implicit none
    private
    
    public :: p1_basis, p1_basis_derivative, basis_1d_t
    
    type :: basis_1d_t
        integer :: n_basis = 2  ! P1 elements have 2 basis functions
        integer :: element_type = 1  ! 1=P1
    contains
        procedure :: init => init_basis_1d
        procedure :: n_dofs => n_dofs_basis_1d
        procedure :: evaluate => evaluate_basis_1d
        procedure :: evaluate_derivative => evaluate_derivative_1d
    end type basis_1d_t
    
contains

    function p1_basis(i, xi) result(phi)
        integer, intent(in) :: i      ! Basis function index (1 or 2)
        real(dp), intent(in) :: xi    ! Local coordinate in [0,1]
        real(dp) :: phi
        
        select case(i)
        case(1)
            phi = 1.0_dp - xi
        case(2)
            phi = xi
        case default
            error stop "Invalid basis function index"
        end select
        
    end function p1_basis
    
    function p1_basis_derivative(i, xi) result(dphi)
        integer, intent(in) :: i      ! Basis function index (1 or 2)
        real(dp), intent(in) :: xi    ! Local coordinate (not used for P1)
        real(dp) :: dphi
        
        associate(dummy => xi)
        end associate
        
        select case(i)
        case(1)
            dphi = -1.0_dp
        case(2)
            dphi = 1.0_dp
        case default
            error stop "Invalid basis function index"
        end select
        
    end function p1_basis_derivative

    subroutine init_basis_1d(this, mesh)
        use mesh_1d, only: mesh_1d_t
        class(basis_1d_t), intent(inout) :: this
        type(mesh_1d_t), intent(in) :: mesh
        
        associate(dummy => mesh)
        end associate
        
        this%n_basis = 2
        this%element_type = 1
    end subroutine init_basis_1d

    function n_dofs_basis_1d(this) result(n)
        class(basis_1d_t), intent(in) :: this
        integer :: n
        n = this%n_basis
    end function n_dofs_basis_1d

    function evaluate_basis_1d(this, i, xi) result(phi)
        class(basis_1d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi
        real(dp) :: phi
        
        associate(dummy => this)
        end associate
        
        phi = p1_basis(i, xi)
    end function evaluate_basis_1d

    function evaluate_derivative_1d(this, i, xi) result(dphi)
        class(basis_1d_t), intent(in) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: xi
        real(dp) :: dphi
        
        associate(dummy => this)
        end associate
        
        dphi = p1_basis_derivative(i, xi)
    end function evaluate_derivative_1d

end module fortfem_basis_1d