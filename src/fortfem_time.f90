module fortfem_time
    !! Canonical direct-import facade for structure-preserving time updates.
    !!
    !! The facade owns no numerical implementation.  Mixed Hamiltonian wave
    !! steps, the symplectic-map diagnostic, and the dissipative Cayley
    !! contract remain in their domain modules so generated and legacy code
    !! have one implementation while new clients avoid the legacy umbrella.
    use fortfem_mixed_wave_time, only: &
        advance_mixed_wave_midpoint, &
        advance_mixed_wave_midpoint_jvp, &
        advance_mixed_wave_midpoint_vjp, &
        advance_mixed_wave_symplectic_euler, &
        advance_mixed_wave_symplectic_euler_jvp, &
        advance_mixed_wave_symplectic_euler_vjp, &
        advance_mixed_wave_strang, &
        advance_mixed_wave_strang_jvp, &
        advance_mixed_wave_strang_vjp
    use fortfem_symplectic_map_defect, only: &
        assemble_symplectic_map_defect, &
        assemble_symplectic_map_defect_jvp, &
        assemble_symplectic_map_defect_vjp
    use fortfem_dissipative_cayley, only: &
        advance_dissipative_cayley, &
        advance_dissipative_cayley_jvp, &
        advance_dissipative_cayley_vjp
    implicit none
    private

    public :: advance_dissipative_cayley
    public :: advance_dissipative_cayley_jvp
    public :: advance_dissipative_cayley_vjp
    public :: advance_mixed_wave_midpoint
    public :: advance_mixed_wave_midpoint_jvp
    public :: advance_mixed_wave_midpoint_vjp
    public :: advance_mixed_wave_strang
    public :: advance_mixed_wave_strang_jvp
    public :: advance_mixed_wave_strang_vjp
    public :: advance_mixed_wave_symplectic_euler
    public :: advance_mixed_wave_symplectic_euler_jvp
    public :: advance_mixed_wave_symplectic_euler_vjp
    public :: assemble_symplectic_map_defect
    public :: assemble_symplectic_map_defect_jvp
    public :: assemble_symplectic_map_defect_vjp

end module fortfem_time
