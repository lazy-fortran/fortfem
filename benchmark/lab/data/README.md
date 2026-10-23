# Round-4 measured inputs of the scoreboard

Copies of the measured sweep outputs that `../scoreboard.py` reads. All are synthetic
manufactured or published-oracle problems; nothing here is non-public data. The files
are byte-identical to the job outputs, except three (the SEAM1 line endings are described below). `fb1-sweep-3356163.csv` holds the
`FB1_TRANSMISSION` rows of the sweep log with its column header and without the
log prefix. `uns1-sweep-3356137.csv` joins the three array-task logs (one
`# HOST/JOB/REVISION` header line each) into one CSV with the host, job and
revision as leading columns. `SHA256SUMS` covers every CSV in
this directory (`sha256sum -c SHA256SUMS`); `scoreboard.py --check` verifies it.
`seam1-l3-3356493.csv` has its CRLF line endings converted to LF (the original
file has the SHA256 listed below); its command columns keep the original cluster
paths as recorded provenance.

| file | job | host | revision | what |
|---|---|---|---|---|
| igaax-scan-3356269.csv | Slurm 3356269 | sCluster node7, exclusive, 16 threads, GCC 12.2 | 3bf52b95 | IGA-AX L1 scan, four matrix-free families |
| igaax-threads-3356269.csv | Slurm 3356269 | same | same | IGA-AX 1/4/16 thread series |
| igaax-disk-3356269.csv | Slurm 3356269 | same, 4 threads | same | IGA-AX disk near-axis table (GLL, B-spline, polar) |
| fb1-sweep-3356163.csv | Slurm 3356163 | sCluster node6, exclusive, GCC 12.2 | 46fefb79 | FB1 FEM+BEM threads 1..16, trace 384..6144 (the sweep log repeats some thread counts) |
| seam1-l3-3356493.csv | Slurm 3356493 | sCluster node17, 4 threads, GCC 12.2 | 817dad00 | SEAM1 L3 disk, aligned mortar seam vs unaligned |
| uns1-sweep-3356137.csv | Slurm array 3356137_0..2 | sCluster node6/13/4, exclusive, 16 threads, GCC 12.2 | 3c2fe688 | UNS1 island vs O-grid sweep, 288 rows |

Source SHA256 of the original job outputs (the same digests appear in `SHA256SUMS`
for the byte-identical copies):

| source | sha256 |
|---|---|
| IGA-AX `scan.csv` | 525da0119088fb9728d7b5cd4ddbb031912b7a2d15699ad87a2d9cf85bf4b681 |
| IGA-AX `threads.csv` | bb40884bffe3b9ab7ce9b92d367f7ab10c563e3022237b9ff290ecf744a25c37 |
| IGA-AX `disk.csv` | 8e15be456aba13d2394ccb8b046a3f679bb0c07a13347b51a1afd98c947b3d01 |
| FB1 `sweep-3356163.log` (rows extracted verbatim) | 8e5599000c9d6b4bbb643866bcf08755666afecaa058ea3f38a68e0b7864db2b |
| SEAM1 `seam-l3.csv` (CRLF original) | 2eafcbcfe4f842a1f8b22e6ec3797978338f4d5d4703c8826e838cd7ca2dcf41 |
| UNS1 `sweep-3356137_0.log` | 9c67f31915d754980b2c6813f765ad950d5ad7cf39d47e133dd4a6fee4a9c824 |
| UNS1 `sweep-3356137_1.log` | 8066e3f629cb0e46aca54b99acf17fafe60433c774932ca8422ee8a09243105f |
| UNS1 `sweep-3356137_2.log` | a7c8387ee1dcf5859b87a03c4aad1ef6b3a5e84311fd5eb279b39de3b0bd8e8b |

The NEQ1 rows are not files: they are the literal rows of the NEQ1 handoff (case A and
case B, node13, 8 threads, jobs 3355912 and 3355913, revision 95ca9ddc) coded in
`scoreboard.py` with their source line.
