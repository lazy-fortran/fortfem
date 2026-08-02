module fortfem_interop
    !! Canonical facade for neutral external-code interoperability.
    !!
    !! This surface groups the metadata-only interchange contracts used by
    !! FEM, BEM, DtN, PML, and independent oracle clients.  It imports the
    !! implementation modules directly; the compatibility umbrella
    !! `fortfem_api` is intentionally not a dependency.  Keeping this facade
    !! narrow makes it safe to migrate downstream clients before broader API
    !! namespaces are reorganized.
    use fortfem_interchange_samples, only: &
        compare_interchange_samples, compare_interchange_samples_jvp, &
        compare_interchange_samples_vjp, initialize_interchange_samples, &
        interchange_sample_set_t, validate_interchange_samples
    use fortfem_complex_interchange_samples, only: &
        compare_complex_interchange_samples, &
        compare_complex_interchange_samples_jvp, &
        compare_complex_interchange_samples_vjp, &
        complex_interchange_sample_set_t, &
        initialize_complex_interchange_samples, &
        validate_complex_interchange_samples
    use fortfem_oracle_manifest, only: &
        oracle_manifest_schema_magic, oracle_manifest_schema_version, &
        oracle_manifest_t, oracle_normalization_t, oracle_timing_t, &
        oracle_tolerance_t, initialize_oracle_manifest, &
        validate_oracle_manifest, read_oracle_manifest, write_oracle_manifest
    use fortfem_boundary_operator_contract, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, BOUNDARY_OPERATOR_BACKEND_BIEST, &
        BOUNDARY_OPERATOR_BACKEND_DTN, BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_BACKEND_NESTOR, BOUNDARY_OPERATOR_BACKEND_PML, &
        BOUNDARY_OPERATOR_BACKEND_USER, &
        BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
        boundary_operator_contract_t, initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    use fortfem_boundary_operator_parity, only: &
        boundary_operator_parity_t, compare_boundary_operator_parity, &
        compare_boundary_operator_parity_jvp, &
        compare_boundary_operator_parity_vjp, &
        validate_boundary_operator_parity
    use fortfem_sheet_current_surface_parity, only: &
        compare_sheet_current_surface_representations, &
        compare_sheet_current_surface_representations_jvp
    use fortfem_sheet_current_parity, only: &
        compare_sheet_current_representations
    use fortfem_regularized_surface_current_layer, only: &
        evaluate_regularized_surface_current_layer, &
        evaluate_regularized_surface_current_layer_jvp, &
        evaluate_regularized_surface_current_layer_vjp, &
        evaluate_regularized_surface_current_integral
    use fortfem_linear_response_interchange, only: &
        assemble_linear_response_operator, &
        assemble_linear_response_operator_jvp, &
        assemble_linear_response_operator_vjp, &
        assemble_linear_response_residual, &
        assemble_linear_response_residual_jvp, &
        assemble_linear_response_residual_vjp, &
        evaluate_linear_response_diagnostics, &
        initialize_linear_response_interchange, &
        linear_response_interchange_t, validate_linear_response_interchange
    use fortfem_linear_response_schema, only: &
        linear_response_schema_magic, read_linear_response_interchange, &
        write_linear_response_interchange
    use fortfem_linear_response_cross, only: &
        assemble_linear_response_eigen_cross_residual, &
        assemble_linear_response_eigen_cross_residual_jvp, &
        assemble_linear_response_eigen_cross_residual_vjp, &
        initialize_linear_response_cross_metadata, &
        linear_response_cross_metadata_t, validate_linear_response_cross_metadata
    use fortfem_linear_perturbation_composition, only: &
        LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE, LINEAR_PASSIVITY_UNDECLARED, &
        LINEAR_RECIPROCITY_TRANSPOSE, LINEAR_RECIPROCITY_UNDECLARED, &
        assemble_linear_perturbation_operator, &
        assemble_linear_perturbation_operator_jvp, &
        assemble_linear_perturbation_operator_vjp, &
        initialize_linear_perturbation_metadata, &
        linear_perturbation_metadata_t, validate_linear_perturbation_metadata
    use fortfem_nonlinear_resistive_mhd_composition, only: &
        RESISTIVE_MHD_AMPERE, RESISTIVE_MHD_ANISOTROPIC_TRANSPORT, &
        RESISTIVE_MHD_BLOCK_COUNT, RESISTIVE_MHD_FARADAY, &
        RESISTIVE_MHD_FREE_BOUNDARY, RESISTIVE_MHD_MOMENTUM, &
        RESISTIVE_MHD_PRESSURE, RESISTIVE_MHD_TENSOR, RESISTIVE_MHD_WALL, &
        assemble_nonlinear_resistive_mhd_residual, &
        assemble_nonlinear_resistive_mhd_residual_jvp, &
        assemble_nonlinear_resistive_mhd_residual_vjp, &
        nonlinear_resistive_mhd_energy_ledger_t
    use fortfem_resistive_mhd_branch_history, only: &
        compare_resistive_mhd_branch_histories, &
        evaluate_resistive_mhd_branch_diagnostics, &
        evaluate_resistive_mhd_branch_path_metric, &
        evaluate_resistive_mhd_branch_path_metric_jvp, &
        evaluate_resistive_mhd_branch_path_metric_vjp, &
        initialize_resistive_mhd_branch_history, &
        resistive_mhd_branch_diagnostics_t, resistive_mhd_branch_history_t, &
        validate_resistive_mhd_branch_history
    implicit none
    private

    public :: BOUNDARY_OPERATOR_BACKEND_BEM
    public :: BOUNDARY_OPERATOR_BACKEND_BIEST
    public :: BOUNDARY_OPERATOR_BACKEND_DTN
    public :: BOUNDARY_OPERATOR_BACKEND_FEM
    public :: BOUNDARY_OPERATOR_BACKEND_NESTOR
    public :: BOUNDARY_OPERATOR_BACKEND_PML
    public :: BOUNDARY_OPERATOR_BACKEND_USER
    public :: BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL
    public :: boundary_operator_contract_t
    public :: boundary_operator_parity_t
    public :: compare_boundary_operator_parity
    public :: compare_boundary_operator_parity_jvp
    public :: compare_boundary_operator_parity_vjp
    public :: compare_complex_interchange_samples
    public :: compare_complex_interchange_samples_jvp
    public :: compare_complex_interchange_samples_vjp
    public :: compare_interchange_samples
    public :: compare_interchange_samples_jvp
    public :: compare_interchange_samples_vjp
    public :: compare_sheet_current_surface_representations
    public :: compare_sheet_current_surface_representations_jvp
    public :: compare_sheet_current_representations
    public :: evaluate_regularized_surface_current_layer
    public :: evaluate_regularized_surface_current_layer_jvp
    public :: evaluate_regularized_surface_current_layer_vjp
    public :: evaluate_regularized_surface_current_integral
    public :: assemble_linear_response_operator
    public :: assemble_linear_response_operator_jvp
    public :: assemble_linear_response_operator_vjp
    public :: assemble_linear_response_residual
    public :: assemble_linear_response_residual_jvp
    public :: assemble_linear_response_residual_vjp
    public :: evaluate_linear_response_diagnostics
    public :: initialize_linear_response_interchange
    public :: linear_response_interchange_t
    public :: validate_linear_response_interchange
    public :: linear_response_schema_magic
    public :: read_linear_response_interchange
    public :: write_linear_response_interchange
    public :: assemble_linear_response_eigen_cross_residual
    public :: assemble_linear_response_eigen_cross_residual_jvp
    public :: assemble_linear_response_eigen_cross_residual_vjp
    public :: initialize_linear_response_cross_metadata
    public :: linear_response_cross_metadata_t
    public :: validate_linear_response_cross_metadata
    public :: LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE
    public :: LINEAR_PASSIVITY_UNDECLARED
    public :: LINEAR_RECIPROCITY_TRANSPOSE
    public :: LINEAR_RECIPROCITY_UNDECLARED
    public :: assemble_linear_perturbation_operator
    public :: assemble_linear_perturbation_operator_jvp
    public :: assemble_linear_perturbation_operator_vjp
    public :: initialize_linear_perturbation_metadata
    public :: linear_perturbation_metadata_t
    public :: validate_linear_perturbation_metadata
    public :: RESISTIVE_MHD_AMPERE
    public :: RESISTIVE_MHD_ANISOTROPIC_TRANSPORT
    public :: RESISTIVE_MHD_BLOCK_COUNT
    public :: RESISTIVE_MHD_FARADAY
    public :: RESISTIVE_MHD_FREE_BOUNDARY
    public :: RESISTIVE_MHD_MOMENTUM
    public :: RESISTIVE_MHD_PRESSURE
    public :: RESISTIVE_MHD_TENSOR
    public :: RESISTIVE_MHD_WALL
    public :: assemble_nonlinear_resistive_mhd_residual
    public :: assemble_nonlinear_resistive_mhd_residual_jvp
    public :: assemble_nonlinear_resistive_mhd_residual_vjp
    public :: nonlinear_resistive_mhd_energy_ledger_t
    public :: compare_resistive_mhd_branch_histories
    public :: evaluate_resistive_mhd_branch_diagnostics
    public :: evaluate_resistive_mhd_branch_path_metric
    public :: evaluate_resistive_mhd_branch_path_metric_jvp
    public :: evaluate_resistive_mhd_branch_path_metric_vjp
    public :: initialize_resistive_mhd_branch_history
    public :: resistive_mhd_branch_diagnostics_t
    public :: resistive_mhd_branch_history_t
    public :: validate_resistive_mhd_branch_history
    public :: complex_interchange_sample_set_t
    public :: initialize_boundary_operator_contract
    public :: initialize_boundary_operator_trace_metadata
    public :: initialize_complex_interchange_samples
    public :: initialize_interchange_samples
    public :: initialize_oracle_manifest
    public :: interchange_sample_set_t
    public :: oracle_manifest_schema_magic
    public :: oracle_manifest_schema_version
    public :: oracle_manifest_t
    public :: oracle_normalization_t
    public :: oracle_timing_t
    public :: oracle_tolerance_t
    public :: read_oracle_manifest
    public :: validate_boundary_operator_contract
    public :: validate_boundary_operator_parity
    public :: validate_complex_interchange_samples
    public :: validate_interchange_samples
    public :: validate_oracle_manifest
    public :: write_oracle_manifest

end module fortfem_interop
