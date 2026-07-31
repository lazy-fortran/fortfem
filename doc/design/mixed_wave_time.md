# Mixed first-order wave time stepping

`advance_mixed_wave_midpoint` advances the compatible semidiscrete state

\[
 M_q\dot q + C^T v = 0, \qquad
 M_v\dot v - Cq = 0.
\]

The state may represent pressure/particle velocity, displacement/momentum,
electric/magnetic variables, or another port-Hamiltonian wave pair. The
implicit-midpoint block is

\[
\begin{bmatrix}
 M_q & \frac{\Delta t}{2}C^T\\
 -\frac{\Delta t}{2}C & M_v
\end{bmatrix}
\begin{bmatrix}q_{n+1}\\v_{n+1}\end{bmatrix}
=
\begin{bmatrix}
 M_q & -\frac{\Delta t}{2}C^T\\
 \frac{\Delta t}{2}C & M_v
\end{bmatrix}
\begin{bmatrix}q_n\\v_n\end{bmatrix}.
\]

For symmetric positive mass blocks, this Cayley map preserves

\[
 E=\tfrac12q^TM_qq+\tfrac12v^TM_vv
\]

and reversing the signed step reverses the state to roundoff. The focused
oscillator test checks the independent Cayley formula, energy preservation,
and reversibility. Dissipative terms are intentionally not hidden in this
routine; they belong in a declared split or metriplectic substep.

The same contract is suitable for tensor-valued pressure, mixed elasticity,
Maxwell, and acoustic-elastic coupling once their compatible mass and
interconnection blocks are assembled.
