"""Regenerate the curved quad/hex MSH 4.1 fixtures (synthetic, no external data).

    python -I make_fixtures.py   # needs the gmsh Python module (PyPI gmsh 4.15.2)
"""
import math
import gmsh


def start(name):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)


def disk(order, optimize=True, name=None):
    start("disk")
    gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.8)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.8)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(order)
    if optimize and order > 1:
        # Cells with two edges on the circle fold at order 3 unless optimised.
        gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
        gmsh.model.mesh.optimize("HighOrder")
    gmsh.write(name or "disk_quad_o%d.msh" % order)
    gmsh.finalize()


def cylinder(order):
    start("cylinder")
    gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.9)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.9)
    gmsh.model.occ.extrude([(2, 1)], 0, 0, 0.5, numElements=[2], recombine=True)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    gmsh.write("cylinder_hex_o%d.msh" % order)
    gmsh.finalize()


def periodic_box(order):
    start("box")
    # Extruded, not transfinite: Gmsh then records the top/bottom node pairs.
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.extrude([(2, 1)], 0, 0, 1, numElements=[3], recombine=True)
    gmsh.model.occ.synchronize()
    eps = 1e-6
    bottom = gmsh.model.getEntitiesInBoundingBox(-eps, -eps, -eps, 1 + eps, 1 + eps,
                                                  eps, 2)
    top = gmsh.model.getEntitiesInBoundingBox(-eps, -eps, 1 - eps, 1 + eps, 1 + eps,
                                               1 + eps, 2)
    shift = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [top[0][1]], [bottom[0][1]], shift)
    for dim, tag in gmsh.model.getEntities(1):
        box = gmsh.model.getBoundingBox(dim, tag)
        if box[5] - box[2] < 0.5:
            gmsh.model.mesh.setTransfiniteCurve(tag, 3)
    gmsh.model.mesh.setTransfiniteSurface(bottom[0][1])
    gmsh.model.mesh.setRecombine(2, bottom[0][1])
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    # Curve the interior so that the geometry order is exercised.
    tags, coords, _ = gmsh.model.mesh.getNodes()
    for i, tag in enumerate(tags):
        x, y, z = coords[3 * i:3 * i + 3]
        bump = 0.05 * math.sin(math.pi * x) * math.sin(math.pi * y)
        gmsh.model.mesh.setNode(tag, [x, y, z + bump * math.sin(2 * math.pi * z)],
                                [])
    gmsh.write("box_periodic_hex_o%d.msh" % order)
    gmsh.finalize()


if __name__ == "__main__":
    for order in (1, 2, 3):
        disk(order)
    # Gmsh's own minSJ is negative in three cells: a detection fixture.
    disk(3, optimize=False, name="disk_quad_o3_inverted.msh")
    cylinder(2)
    periodic_box(2)
    periodic_box(3)
