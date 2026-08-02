module fortfem_boundary
    !! Direct boundary-facing facade.
    !!
    !! Mesh boundary geometry remains available here for compatibility.  The
    !! open-boundary primitives are re-exported from their canonical modules
    !! so clients can depend on this small boundary surface without importing
    !! the umbrella ``fortfem_api`` (or any parity/comparison implementation).
    use fortfem_boundary_operator_contract, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_BACKEND_NESTOR, &
        BOUNDARY_OPERATOR_BACKEND_BIEST, &
        BOUNDARY_OPERATOR_BACKEND_PML, &
        BOUNDARY_OPERATOR_BACKEND_USER, &
        BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    use fortfem_circular_dtn_2d, only: &
        apply_circular_helmholtz_dtn, circular_helmholtz_dtn_eigenvalue
    use fortfem_kinds, only: dp
    use fortfem_laplace_boundary_operators_2d, only: &
        assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant
    use fortfem_helmholtz_exterior_2d, only: &
        evaluate_helmholtz_combined_potential_constant, &
        solve_helmholtz_cfie_constant
    use fortfem_laplace_symmetric_coupling_2d, only: &
        solve_laplace_symmetric_coupling_p1_p0
    use fortfem_planar_helmholtz_dtn, only: &
        apply_planar_helmholtz_dtn, &
        apply_planar_helmholtz_dtn_jvp, &
        apply_planar_helmholtz_dtn_vjp, &
        assemble_planar_helmholtz_dtn_form, &
        assemble_planar_helmholtz_dtn_form_jvp, &
        assemble_planar_helmholtz_dtn_form_vjp
    implicit none

    private
    public :: boundary_t
    public :: BOUNDARY_OPERATOR_BACKEND_BEM
    public :: BOUNDARY_OPERATOR_BACKEND_DTN
    public :: BOUNDARY_OPERATOR_BACKEND_FEM
    public :: BOUNDARY_OPERATOR_BACKEND_NESTOR
    public :: BOUNDARY_OPERATOR_BACKEND_BIEST
    public :: BOUNDARY_OPERATOR_BACKEND_PML
    public :: BOUNDARY_OPERATOR_BACKEND_USER
    public :: BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL
    public :: boundary_operator_contract_t
    public :: initialize_boundary_operator_contract
    public :: initialize_boundary_operator_trace_metadata
    public :: validate_boundary_operator_contract
    public :: apply_planar_helmholtz_dtn
    public :: apply_planar_helmholtz_dtn_jvp
    public :: apply_planar_helmholtz_dtn_vjp
    public :: assemble_planar_helmholtz_dtn_form
    public :: assemble_planar_helmholtz_dtn_form_jvp
    public :: assemble_planar_helmholtz_dtn_form_vjp
    public :: apply_circular_helmholtz_dtn
    public :: circular_helmholtz_dtn_eigenvalue
    public :: assemble_laplace_single_layer_constant
    public :: assemble_laplace_hypersingular_linear
    public :: assemble_helmholtz_double_layer_constant
    public :: assemble_helmholtz_hypersingular_linear
    public :: assemble_helmholtz_single_layer_constant
    public :: evaluate_helmholtz_combined_potential_constant
    public :: solve_helmholtz_cfie_constant
    public :: solve_laplace_symmetric_coupling_p1_p0

    ! Boundary type for defining domains
    type :: boundary_t
        integer :: n_points = 0
        real(dp), allocatable :: points(:,:) ! (2, n_points)
        integer, allocatable :: labels(:) ! (n_points-1) segment labels
        logical :: is_closed = .false.
    end type boundary_t

end module fortfem_boundary
