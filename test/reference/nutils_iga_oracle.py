"""Independent MIT-licensed Nutils oracle for quadratic IGA operators.

Nutils documentation: https://docs.nutils.org/
Pinned oracle version: 9.2
"""

import json

import numpy
from nutils import function, mesh


topology, geometry = mesh.rectilinear(
    [numpy.linspace(0.0, 1.0, 4), numpy.linspace(0.0, 1.0, 3)]
)
basis = topology.basis("spline", degree=2)
namespace = function.Namespace()
namespace.x = geometry
namespace.basis = basis

mass = topology.integral(
    namespace.eval_nm("basis_n basis_m d:x"), degree=6
).eval(legacy=True)
stiffness = topology.integral(
    namespace.eval_nm("basis_n,i basis_m,i d:x"), degree=6
).eval(legacy=True)
coordinate_load = topology.integral(
    namespace.eval_n("basis_n x_0 d:x"), degree=6
).eval()

constant_coefficients = numpy.ones(basis.shape[0])
coordinate_coefficients = mass.solve(coordinate_load)
mass_energy = float(constant_coefficients @ (mass @ constant_coefficients))
coordinate_energy = float(
    coordinate_coefficients @ (stiffness @ coordinate_coefficients)
)
constant_residual = float(
    numpy.linalg.norm(stiffness @ constant_coefficients, ord=numpy.inf)
)

assert basis.shape == (20,)
assert abs(mass_energy - 1.0) < 2.0e-13
assert abs(coordinate_energy - 1.0) < 2.0e-12
assert constant_residual < 2.0e-13

print(
    json.dumps(
        {
            "oracle": "nutils",
            "version": "9.2",
            "basis_dofs": basis.shape[0],
            "constant_mass_energy": mass_energy,
            "coordinate_stiffness_energy": coordinate_energy,
            "constant_stiffness_residual": constant_residual,
        },
        sort_keys=True,
    )
)
