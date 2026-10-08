
The reference-basis generator also emits
`fortfem_reference_scalar_intervals.f90`: rigorous value/physical-gradient
reconstruction for a scalar P2 field on an affine triangle. It uses the same
canonical basis and centered nodal differences as the floating products.
Exact degree-zero through degree-two polynomial reproduction is proved
natively before emission. The runtime requires the existing FortNum interval
provider; callers must reject singular geometry and validate mesh conformity.
The independent native test checks a physical quadratic, reversed winding,
translation, representable gauge shifts and incorrect edge ordering. Existing
CAPI/CMake interfaces do not consume this interval module; a consumer must
provide FortNum when compiling this additional source.
