This example interpolates the manufactured curl-free field

\[
\boldsymbol{E}=\nabla\!\left(x^4+\frac{2}{3}y^3+z^3\right)
 = (4x^3,\,2y^2,\,3z^2)
\]

with first-kind Nédélec spaces on the reference tetrahedron.  The gallery
starts with the computed solution: the 3-D view shows the tetrahedron edges,
field magnitude, and direction arrows; a triangular slice shows the same
vector field with quiver arrows.  The p-convergence curve is reported after
the physical solution plots.  Each run writes `gallery_sequence.txt`, whose
two records are an execution oracle for that ordering.
