# FortFEM Mesh Benchmark

This example benchmarks FortFEM's mesh generation performance against FreeFEM.

## Description

Compares mesh generation speed and quality between FortFEM and FreeFEM implementations for various mesh sizes and geometries.

## Usage

```bash
fpm run --example fortfem_mesh_benchmark
```

## What it does

- Generates meshes of different sizes using FortFEM
- Measures timing and mesh quality metrics
- Outputs benchmark results for comparison with FreeFEM