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
    use fortfem_fourier_vector_product, only: &
        apply_fourier_bilinear_product, &
        apply_fourier_bilinear_product_jvp, &
        apply_fourier_bilinear_product_vjp, &
        assemble_fourier_vector_product, &
        assemble_fourier_vector_product_jvp, &
        assemble_fourier_vector_product_vjp
    use fortfem_fourier_mode_linear_operator, only: &
        apply_fourier_mode_linear_operator, &
        apply_fourier_mode_linear_operator_jvp, &
        apply_fourier_mode_linear_operator_vjp
    use fortfem_fourier_mode_energy, only: &
        assemble_fourier_mode_energy, &
        assemble_fourier_mode_energy_jvp, &
        assemble_fourier_mode_energy_vjp, &
        fourier_coefficients_conjugate_symmetric
    use fortfem_axis_regular_fourier_modes, only: &
        AXIS_RADIAL_PARITY_EVEN, AXIS_RADIAL_PARITY_ODD, &
        axis_regular_mode_record_t, axis_regular_mode_table_t, &
        axis_regular_mode_requirements, build_axis_regular_mode_table, &
        validate_axis_regular_mode_table
    use fortfem_toroidal_harmonics, only: &
        toroidal_p, toroidal_q, toroidal_p_derivative, toroidal_q_derivative, &
        toroidal_p_second_derivative, toroidal_q_second_derivative
    use fortfem_magnetic_curvilinear_coefficients_2d, only: &
        scalar_reluctivity_curvilinear_fourier_coefficients
    use fortfem_toroidal_poisson_dtn, only: &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_ampere_field_p, &
        toroidal_poisson_exterior_dtn_p
    use fortfem_toroidal_poisson_dtn_ad, only: &
        evaluate_toroidal_harmonic_p_jvp, evaluate_toroidal_harmonic_p_vjp, &
        evaluate_toroidal_ampere_field_p_jvp, &
        evaluate_toroidal_ampere_field_p_vjp, &
        toroidal_poisson_exterior_dtn_p_jvp, &
        toroidal_poisson_exterior_dtn_p_vjp
    use fortfem_toroidal_spectral_trace, only: &
        evaluate_toroidal_spectral_trace, &
        evaluate_toroidal_spectral_trace_jvp, &
        evaluate_toroidal_spectral_trace_vjp, &
        evaluate_toroidal_spectral_trace_grid, &
        evaluate_toroidal_spectral_trace_grid_jvp, &
        evaluate_toroidal_spectral_trace_grid_vjp, &
        solve_toroidal_spectral_neumann, &
        solve_toroidal_spectral_neumann_jvp, &
        solve_toroidal_spectral_neumann_vjp
    use fortfem_toroidal_modal_convolution, only: &
        apply_toroidal_modal_convolution, &
        apply_toroidal_modal_convolution_jvp, &
        apply_toroidal_modal_convolution_vjp
    use fortfem_assembly_bspline_2d, only: apply_toroidal_fourier_derivative
    use fortfem_harmonic_period_normalization, only: &
        normalize_harmonic_one_forms, normalize_harmonic_one_forms_jvp, &
        normalize_harmonic_one_forms_vjp
    implicit none
    private

    public :: build_fourier_mode_closure_registry
    public :: build_axis_regular_mode_table
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
    public :: apply_fourier_bilinear_product
    public :: apply_fourier_bilinear_product_jvp
    public :: apply_fourier_bilinear_product_vjp
    public :: assemble_fourier_vector_product
    public :: assemble_fourier_vector_product_jvp
    public :: assemble_fourier_vector_product_vjp
    public :: apply_fourier_mode_linear_operator
    public :: apply_fourier_mode_linear_operator_jvp
    public :: apply_fourier_mode_linear_operator_vjp
    public :: assemble_fourier_mode_energy
    public :: assemble_fourier_mode_energy_jvp
    public :: assemble_fourier_mode_energy_vjp
    public :: fourier_coefficients_conjugate_symmetric
    public :: axis_regular_mode_requirements
    public :: axis_regular_mode_record_t
    public :: axis_regular_mode_table_t
    public :: find_fourier_mode
    public :: fourier_mode_conjugate_index
    public :: fourier_mode_registry_t
    public :: AXIS_RADIAL_PARITY_EVEN
    public :: AXIS_RADIAL_PARITY_ODD
    public :: fourier_mode_triad_closed
    public :: initialize_fourier_mode_registry
    public :: toroidal_p
    public :: toroidal_p_derivative
    public :: toroidal_p_second_derivative
    public :: toroidal_q
    public :: toroidal_q_derivative
    public :: toroidal_q_second_derivative
    public :: scalar_reluctivity_curvilinear_fourier_coefficients
    public :: evaluate_toroidal_harmonic_p
    public :: evaluate_toroidal_harmonic_p_jvp
    public :: evaluate_toroidal_harmonic_p_vjp
    public :: evaluate_toroidal_ampere_field_p
    public :: evaluate_toroidal_ampere_field_p_jvp
    public :: evaluate_toroidal_ampere_field_p_vjp
    public :: toroidal_poisson_exterior_dtn_p
    public :: toroidal_poisson_exterior_dtn_p_jvp
    public :: toroidal_poisson_exterior_dtn_p_vjp
    public :: evaluate_toroidal_spectral_trace
    public :: evaluate_toroidal_spectral_trace_jvp
    public :: evaluate_toroidal_spectral_trace_vjp
    public :: evaluate_toroidal_spectral_trace_grid
    public :: evaluate_toroidal_spectral_trace_grid_jvp
    public :: evaluate_toroidal_spectral_trace_grid_vjp
    public :: solve_toroidal_spectral_neumann
    public :: solve_toroidal_spectral_neumann_jvp
    public :: solve_toroidal_spectral_neumann_vjp
    public :: apply_toroidal_modal_convolution
    public :: apply_toroidal_modal_convolution_jvp
    public :: apply_toroidal_modal_convolution_vjp
    public :: apply_toroidal_fourier_derivative
    public :: normalize_harmonic_one_forms
    public :: normalize_harmonic_one_forms_jvp
    public :: normalize_harmonic_one_forms_vjp
    public :: validate_fourier_mode_registry
    public :: validate_axis_regular_mode_table

end module fortfem_fourier
