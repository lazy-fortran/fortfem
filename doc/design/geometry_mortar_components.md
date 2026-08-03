# Component-valued IGA mortar traces

`assemble_geometry_mortar_component_coupling` is the neutral trace
contraction used when two spline patches carry vector or flattened tensor
unknowns.  It does not choose a Piola convention, material law, or physical
coordinate system: those remain caller-owned providers.

For a fixed reference quadrature row `q`, test degree of freedom `i`, and
trial degree of freedom `j`, the assembled cross-mass is

\[
 M_{ij}=\sum_q w_q J_q
       \sum_{a,b} T_{a q i}\,G_{a b q}\,R_{b q j}.
\]

The trace arrays have shape `(component, quadrature, dof)`.  `G` has shape
`(test_component, trial_component, quadrature)` and can be a pulled-back
metric, a constitutive tensor, or a caller-provided component transform.
`w_q` must be positive reference weights and `J_q` must be a positive physical
surface Jacobian.  The routine returns the physical weights `w_q J_q` as a
separate diagnostic, which lets a caller audit orientation and geometry
without inspecting the assembled matrix.

The JVP applies the complete product rule to both traces, the component
metric, and the physical measure.  The VJP uses the real dot-product
convention and returns cotangents for every input.  Topology, patch
identification, signed H(curl)/H(div) maps, and Piola transformations are
deliberately outside this contraction; their outputs can be supplied as the
trace and component arrays.

`test_geometry_mortar_component_coupling` compares the value against an
independent nested-loop tensor oracle, checks central differences for the
JVP, and checks the real adjoint identity for the VJP.  It exercises both a
two-component vector trace and a four-component flattened tensor trace.
