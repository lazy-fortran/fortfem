# Sheet-current surface gallery

This physical-first gallery samples the same caller-owned Ampère trace on a
slab, cylinder, sphere, and torus.  Each view shows the oriented surface and
short arrows for the tangential surface current
\(\mathbf K=\mathbf n\times[\![\mathbf H]\!]\).  The numerical layer is a
normalized Gaussian \(\mathbf K\,\delta_\epsilon\) in signed normal distance;
the fitted and resolved integrated ledgers are compared with FortFEM's
canonical `fortfem_interop` surface-parity contract.

The fields and geometries are manufactured, closure-neutral fixtures.  No
equilibrium, coil, plasma profile, BEM, or material model is inferred.  The
independent gate checks positive quadrature measures, unit outward normals,
tangentiality, all four geometry labels, integrated-current parity, and the
fixed-topology physical-before-diagnostics sequence.
