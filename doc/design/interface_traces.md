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
penalty portion of Nitsche coupling.  The public
`assemble_symmetric_nitsche_interface` routine adds the common-normal average
flux consistency terms and is checked against an independent block oracle.
