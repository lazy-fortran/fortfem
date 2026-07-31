---
title: Internal-manifold graph contract
---

# Internal-manifold graph contract

`fortfem_internal_manifold_graph` is the neutral topology layer for interfaces
that lie inside a volume. Each oriented manifold records its plus and minus
region and an optional pair of boundary/junction IDs. It deliberately owns no
coordinates, mesh, trace space, pressure law, current model, or input format.

An endpoint ID of zero means that the manifold is boundaryless. Equal nonzero
endpoints represent a periodic self-identification. Distinct nonzero endpoints
represent an open manifold whose two boundary components meet the corresponding
junctions. The junction incidence convention is

\[
    B_{j,m}=+1\quad\text{at the end junction},\qquad
    B_{j,m}=-1\quad\text{at the start junction}.
\]

## API

```fortran
type(internal_manifold_graph_t) :: graph
call initialize_internal_manifold_graph( &
    graph, region_count, junction_count, plus_region, minus_region, &
    start_junction, end_junction, status)
call internal_manifold_graph_region_incidence(graph, G, status)
call internal_manifold_graph_junction_incidence(graph, B, status)
call internal_manifold_graph_closed(graph, closed, status)
call internal_manifold_graph_components(graph, components, count, status)
```

`G` is the oriented region incidence used by conservative interface laws.
`B` is the boundary ledger topology used by later surface-current divergence
and junction-balance operators. Closed manifolds have a zero column in `B`.
Components are labels of the manifold patches connected through a common
junction; boundaryless manifolds with no junction are separate components.

## Independent fixtures

`test_internal_manifold_graph` covers an open slab chain, a periodic cylinder,
a boundaryless sphere, a toroidal periodic component, orientation signs,
component labels, and malformed endpoint data. These are topology-only tests;
surface geometry and physical jump/pressure laws remain client-owned.
