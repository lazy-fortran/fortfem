module fortfem
    ! Main module that exports the clean FEniCS-style API
    use fortfem_api
    use fortfem_simple_solvers
    implicit none
    
    ! Re-export only the clean API
    public
    
end module fortfem
