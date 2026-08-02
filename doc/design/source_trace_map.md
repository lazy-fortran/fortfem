# Source-to-trace map

`fortfem_source_trace_map` is a small, neutral dense-map contract. Given a
caller-owned complex matrix `A` and source coefficients `x`,

```text
trace = A x.
```

The value routine returns an allocatable trace and validates dimensions and
finite complex values. Its JVP is the exact product rule

```text
trace_dot = A_dot x + A x_dot.
```

The VJP uses the real-part complex pairing
`Re(sum(conjg(trace_bar) * trace_dot))`. Consequently,

```text
A_bar = trace_bar conjg(x)^T,
x_bar = A^H trace_bar.
```

All reverse outputs are allocatable and are reset to zero-sized arrays on a
rejected call, making repeated calls safe without caller-side deallocation.

`evaluate_weighted_source_trace_reciprocity_defect` compares two supplied
source/target experiments. For `trace_k = A source_k`, positive target-side
weights `W`, and the transpose (not Hermitian) reciprocity pairing, it returns

```text
|target_1^T W trace_2 - target_2^T W trace_1|
------------------------------------------------
max(1, |target_1^T W trace_2|, |target_2^T W trace_1|).
```

No geometry, Green kernel, constitutive law, source normalization, or physical
interpretation is inferred. Such conventions remain with the caller.
