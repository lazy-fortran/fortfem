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
        compare_sheet_current_surface_representations
    use fortfem_sheet_current_parity, only: &
        compare_sheet_current_representations
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
    public :: compare_sheet_current_representations
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
