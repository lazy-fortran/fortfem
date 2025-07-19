module fortfem_boundary_assembly_2d
    ! Interface module for 2D boundary assembly
    use boundary_assembly_2d_module, only: boundary_assembly_2d_t, assemble_neumann_bc
    implicit none
    public :: boundary_assembly_2d_t, assemble_neumann_bc
end module fortfem_boundary_assembly_2d