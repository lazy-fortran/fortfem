# FreeFEM Ampere interoperability runner

This runner independently solves the shared unit-square Ampere problem:

\[
\operatorname{curl}\operatorname{curl}\mathbf E+\mathbf E=\mathbf f,
\qquad
\mathbf E\times\mathbf n=0,
\]

with

\[
\mathbf E=(\sin(\pi y),\sin(\pi x)),\qquad
\mathbf f=(\pi^2+1)\mathbf E.
\]

The field has exactly zero tangential trace on all four sides and
\(\operatorname{curl}\operatorname{curl}\mathbf E=\pi^2\mathbf E\). These
identities are generated and tested independently in FortFEM by
`gen_interoperability_oracles`.

Run it with a pinned FreeFEM installation:

```sh
FreeFem++ curl_curl_2d.edp
```

The runner writes `ampere.csv`. Generated data and images are ignored; they
belong in CI artifacts and, after review, the external benchmark-data
repository. It exits nonzero if either analytical error fails to decrease.

`RT0Ortho` is FreeFEM's lowest-order triangular H(curl) element. The boundary
assignment acts on its tangential edge degrees of freedom.
