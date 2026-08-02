module fortfem_fourier
    !! Canonical facade for Fourier-FEM mode metadata and synthesis.
    !!
    !! This module is deliberately a re-export layer: mode registry and
    !! expansion implementations remain in their domain modules, while the
    !! pinned half-integer toroidal harmonics remain owned by FortNum through
    !! ``fortfem_toroidal_harmonics``.  New consumers can import this coherent
    !! Fourier surface without depending on the compatibility umbrella
    !! ``fortfem_api``.
    use fortfem_fourier_mode_registry, only: &
        fourier_mode_registry_t, initialize_fourier_mode_registry, &
        validate_fourier_mode_registry, evaluate_fourier_mode, &
        evaluate_fourier_mode_jvp, evaluate_fourier_mode_vjp, &
        find_fourier_mode, fourier_mode_conjugate_index, &
        fourier_mode_triad_closed, build_fourier_mode_triad_map, &
        build_fourier_mode_padded_registry, build_fourier_mode_closure_registry
    use fortfem_fourier_mode_expansion, only: &
        evaluate_fourier_mode_expansion, &
        evaluate_fourier_mode_expansion_jvp, &
        evaluate_fourier_mode_expansion_vjp, &
        evaluate_fourier_mode_expansion_hessian, &
        evaluate_fourier_mode_expansion_hvp
    use fortfem_toroidal_harmonics, only: &
        toroidal_p, toroidal_q, toroidal_p_derivative, toroidal_q_derivative, &
        toroidal_p_second_derivative, toroidal_q_second_derivative
    implicit none
    private

    public :: build_fourier_mode_closure_registry
    public :: build_fourier_mode_padded_registry
    public :: build_fourier_mode_triad_map
    public :: evaluate_fourier_mode
    public :: evaluate_fourier_mode_expansion
    public :: evaluate_fourier_mode_expansion_hessian
    public :: evaluate_fourier_mode_expansion_hvp
    public :: evaluate_fourier_mode_expansion_jvp
    public :: evaluate_fourier_mode_expansion_vjp
    public :: evaluate_fourier_mode_jvp
    public :: evaluate_fourier_mode_vjp
    public :: find_fourier_mode
    public :: fourier_mode_conjugate_index
    public :: fourier_mode_registry_t
    public :: fourier_mode_triad_closed
    public :: initialize_fourier_mode_registry
    public :: toroidal_p
    public :: toroidal_p_derivative
    public :: toroidal_p_second_derivative
    public :: toroidal_q
    public :: toroidal_q_derivative
    public :: toroidal_q_second_derivative
    public :: validate_fourier_mode_registry

end module fortfem_fourier
