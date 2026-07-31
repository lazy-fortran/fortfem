# Interface traces, jumps, and averages

`compute_interface_scalar_jump_average` and
`compute_interface_vector_traces` establish the common orientation convention
for fitted, cut, DG, BEM, and plasma-interface formulations:

- the supplied unit normal points from the minus side to the plus side;
- every jump is plus minus minus;
- averages are arithmetic averages;
- normal traces are (n\cdot u);
- tangential traces are (u-(n\cdot u)n);
- the rotated tangential jump is (n\times(u^+-u^-)).

The last quantity is the natural Ampere surface-current trace
(\mathbf K=\mathbf n\times[\![\mathbf H]\!]).  The routines validate vector
shapes, finite values, and unit normals, and return a status rather than
silently changing orientation.

The focused test uses an independent vector oracle, checks tangential
orthogonality and the rotated-current sign, and rejects a non-unit normal.
Broken finite-element spaces, surface measures, distributional source terms,
and Nitsche/mortar coupling will consume this contract in later slices.  The
public `assemble_interface_jump_penalty` routine now assembles the symmetric
positive-semidefinite \([T^+,-T^-]^T[T^+,-T^-]\) block used by SIPG and the
penalty portion of Nitsche coupling. Its JVP propagates plus/minus traces,
surface weights, and penalty directions; its VJP returns all corresponding
cotangents while retaining the fixed orientation. The public
`assemble_symmetric_nitsche_interface` routine adds the common-normal average
flux consistency terms and is checked against an independent block oracle. Its
product-rule JVP propagates trace, flux, surface-weight, and penalty
directions; its VJP returns the corresponding cotangents for arbitrary real
matrix seeds, so fitted and cut interface geometry can remain differentiable.
`assemble_mortar_trace_coupling` supplies the weighted cross-mass block when
the test and trial trace spaces have different degree-of-freedom counts. Its
JVP propagates trace and surface-weight directions, while its VJP returns the
test-trace, trial-trace, and quadrature-weight cotangents. The focused AD test
checks the real dot-product identity and rejects incompatible directions.
`assemble_scalar_sipg_interface` generalizes the same contract to independent
test and trial traces. The consistency selector `theta=1,0,-1` gives symmetric,
incomplete, and nonsymmetric interior-penalty blocks, with value/JVP/VJP
actions for all trace, flux, weight, and penalty inputs.
`assemble_surface_triangle_measures_3d` supplies the matching oriented area
and unit-normal arrays for a triangular fitted surface. Its geometry JVP and
VJP accumulate shared-vertex derivatives and are checked by independent
finite-difference and real dot-product oracles, so current-sheet and delta-load
primitives can consume one orientation-preserving measure contract.
`evaluate_level_set_triangle_interface_2d` supplies the corresponding first
unfitted-manifold primitive: two edge intersections, physical segment length,
and the level-set-gradient normal for a linear triangle cut. Vertex crossings
are treated as topology events for future derivative and cut-cell paths.
