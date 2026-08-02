# Total-pressure interface balance

`assemble_total_pressure_jump` is the closure-neutral scalar interface block

\[
  p^+ + |B^+|^2/(2\mu)-p^- - |B^-|^2/(2\mu)-t=0.
\]

The caller supplies both pressure and field traces, permeability, and the
target jump.  Exact real JVP/VJP products make it usable with CGL, elastic,
Maxwell, or multi-region clients without selecting a plasma law or geometry.
