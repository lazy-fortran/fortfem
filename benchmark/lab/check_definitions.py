#!/usr/bin/env python3
"""Independent exact-algebra oracle for the documented MMS (SymPy, not FortSym)."""
import sympy as sp

x, y, z, r, phi = sp.symbols('x y z r phi', real=True)
U = sp.sin(x)+sp.cos(y)+sp.sin(z)+(x+y+z)/10
f = sp.sin(x)+sp.cos(y)+sp.sin(z)
assert sp.simplify(-sum(sp.diff(U, q, 2) for q in (x, y, z))-f) == 0
A = sp.Matrix([sp.sin(y), sp.sin(z), sp.sin(x)])


def curl(a):
    return sp.Matrix([sp.diff(a[2], y)-sp.diff(a[1], z),
                      sp.diff(a[0], z)-sp.diff(a[2], x),
                      sp.diff(a[1], x)-sp.diff(a[0], y)])


assert curl(curl(A))+A == 2*A
psi = sp.cos(sp.pi*x)*sp.cos(sp.pi*y)
b = sp.Matrix([-sp.diff(psi, y), sp.diff(psi, x)])
assert sp.expand(b.dot(sp.Matrix([sp.diff(psi, x), sp.diff(psi, y)]))) == 0
assert sp.simplify(-sp.diff(psi, x, 2)-sp.diff(psi, y, 2)-2*sp.pi**2*psi) == 0
# Prove tangency on the full disk without a polar chart or normalisation.
a = 1+sp.cos(5*phi)/10
rc = 3+sp.cos(5*phi)/10
zc = sp.sin(5*phi)/10
u, v = (r-rc)/a, (z-zc)*sp.Rational(10, 13)
s = u*u+v*v
iota = (2+s)/5
t = sp.Matrix([-sp.sin(5*phi)*(1+u)/2-a*iota*v,
               r, sp.cos(5*phi)/2+sp.Rational(13, 10)*iota*u])
grad_s = sp.Matrix([sp.diff(s, r), sp.diff(s, phi)/r, sp.diff(s, z)])
assert sp.trigsimp(sp.expand(t.dot(grad_s))) == 0
# Positive coupled model, and its exact physical forcing.
C = sp.Matrix([[2, 1, 0], [1, 2, 0], [0, 0, 1]])
assert C.eigenvals() == {3: 1, 1: 2}
w = sp.Matrix([U, 2*U, 3*U])
F = C*sp.Matrix([f, 2*f, 3*f])+w
for j in range(3):
    assert sp.simplify(-sum(C[j,k]*sum(sp.diff(w[k], q, 2)
                       for q in (x,y,z)) for k in range(3))+w[j]-F[j]) == 0
print('PASS: exact L1, L2, L3a, L3b tangency and L4 forcing identities')
