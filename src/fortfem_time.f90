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
    use fortfem_assembly_bspline_2d, only: &
        advance_bspline_jorek_poloidal_flux_midpoint_steps
    use fortfem_mixed_wave_wall_time, only: &
        advance_mixed_wave_wall_midpoint, &
        advance_mixed_wave_wall_midpoint_jvp, &
        advance_mixed_wave_wall_midpoint_vjp, &
        evaluate_mixed_wave_wall_energy_balance
    use fortfem_continuation_event, only: &
        CONTINUATION_EVENT_NEAR_ZERO, CONTINUATION_EVENT_SIGN_CROSSING, &
        classify_continuation_event
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
    public :: advance_bspline_jorek_poloidal_flux_midpoint_steps
    public :: advance_mixed_wave_wall_midpoint
    public :: advance_mixed_wave_wall_midpoint_jvp
    public :: advance_mixed_wave_wall_midpoint_vjp
    public :: assemble_symplectic_map_defect
    public :: assemble_symplectic_map_defect_jvp
    public :: assemble_symplectic_map_defect_vjp
    public :: evaluate_mixed_wave_wall_energy_balance
    public :: CONTINUATION_EVENT_NEAR_ZERO
    public :: CONTINUATION_EVENT_SIGN_CROSSING
    public :: classify_continuation_event

end module fortfem_time
