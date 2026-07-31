---
title: Region/interface graph contract
---

# Region/interface graph contract

`fortfem_region_interface_graph` is the neutral graph layer between an
oriented cell complex and application-owned interface laws. It records which
side of each internal interface is the plus or minus region. It does not own a
pressure law, a current model, a mesh, a plasma profile, or an input-file
format.

## Oriented incidence

For interface \(e\), the region incidence column is

\[
G_{:,e}=\mathbf e_{r_+}-\mathbf e_{r_-}.
\]

The public routines construct this integer matrix and validate that every
endpoint is a valid region ID. Reversing an interface reverses its column and
therefore gives a direct sign oracle for trace and sheet-current assembly.
Equal endpoints are allowed: they represent a periodic self-identification
and produce a zero region-boundary column without silently inventing a second
region.

## API

```fortran
type(region_interface_graph_t) :: graph
integer, allocatable :: incidence(:, :), components(:)
integer :: component_count, status

call initialize_region_interface_graph( &
    graph, region_count, plus_region, minus_region, status)
call validate_region_interface_graph(graph, status)
call region_interface_graph_incidence(graph, incidence, status)
call region_interface_graph_components( &
    graph, components, component_count, status)
```

`components` contains compact labels starting at one. The component routine
uses a union-find traversal over the undirected region adjacency; interface
orientation does not affect connectivity. The incidence matrix retains the
orientation for later conservative balance laws.

## Independent fixtures

`test_region_interface_graph` checks a three-region slab-like chain, two
disconnected region pairs, a periodic self-interface, reversed plus/minus
signs through the incidence oracle, and malformed endpoint tables. These are
topological tests only. Surface quadrature, trace bases, sheet currents,
pressure balance, and application-owned boundary laws remain separate
contracts.
