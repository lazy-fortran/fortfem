# Manufactured linear perturbation response

This closure-neutral fixture composes the seven public linear-perturbation
blocks

\[
A(\omega)=L+P+V+W+S-\omega^2M+i\omega R
\]

and solves `A state = source` for five coupled Fourier coefficients under an
asymmetric complex drive. The first plot reconstructs the complex state over
physical poloidal angle. It is a
manufactured displacement-like response, not a plasma equilibrium, stability
model, singular-layer closure, or plasma-code file reader.

The frequency sweep forms the dense response matrix `A^-1` and passes it to
FortFEM's public reciprocity/passivity diagnostic. It also compares the exact
frequency JVP of the seven-block composition with an independent centered
difference and checks the forced residual after each solve.

Outputs:

- `linear_response_state_1d.png`: physical complex state first, with magnitude
  and real/imaginary components over poloidal angle;
- `linear_response_phase_1d.png`: unwrapped-free principal phase of that state;
- `linear_response_diagnostics_1d.png`: reciprocity, passivity, derivative,
  and forced-residual diagnostics across frequency;
- `linear_response_state.csv` and `linear_response_diagnostics.csv`:
  reproducible plot sources;
- `benchmark.txt`: small runner-local composition, solve, and diagnostic timing
  record.

The timings are regression context for this tiny manufactured algebraic case,
not production performance claims.
