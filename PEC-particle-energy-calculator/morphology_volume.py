"""
Calculates the volume and surface of the crystal morphology defined from
(hkl index, normal distance) data. The unit of measurement of the results
depends on the input data. Indices of symmetry-equivalent facets need
to be explicitly included.
"""

import numpy as np
import itertools
from scipy.spatial import ConvexHull
from ccdc.io import EntryReader


def find_vertices(points, facet_equations, tol=1E-6):
    """
    Limit plane intersection points to morphology vertices by checking if
    they are on the same or opposite side of origin. Check is made by
    comparing sign.

    :param points: all point intersections of the morphology facets
    :param facet_equations: the equations of the facets as (a, b, c, d)
    :param tol: the tolerance to check if a point lies on plane
    :return: the unique vertices of the morphology as an array
    """

    def eval_points(points, plane):
        """Plug the point coordinates in the plane equation"""
        return np.matmul(points, plane[:-1]) + plane[-1]

    center_point = [0, 0, 0]

    points = np.array(points)
    morphology_vertices = []

    for i, plane in enumerate(facet_equations):
        # compare the signs of center/plane vs point/plane
        if i == 0:
            working_array = points
        else:
            working_array = morphology_vertices
        test_c = np.sign(eval_points(center_point, np.array(plane)))
        p_planes = eval_points(working_array, np.array(plane))
        # make sure to get rid of negative zeros
        test_p = np.sign(np.where(abs(p_planes) < 0 + tol, 0, p_planes))
        # we then keep point that are either on the plane (test_p = 0) or on the same side of the plane
        filtered_points = np.array(working_array)[np.where((test_p == 0) | (test_p == test_c))[0]]
        morphology_vertices = filtered_points

    return np.array(morphology_vertices)


def morphology_properties(growth_rates):
    """
    Calculate the volume and surface of the morphology.
    Symmetry-equivalent facets need to be explicit.
    The units of measurement of the outputs depend on the input distances.

    :param growth_rates: an iterable of pairs, ccdc.crystal.Crystal.MillerIndices and perpendicular distance
    :return: the morphology volume and surface
    """

    # the facet equations in the form "ax + by + cz + d = 0"
    facet_equations = []
    for item in growth_rates:
        a, b, c = item[0].plane.normal[0], item[0].plane.normal[1], item[0].plane.normal[2]
        facet_equations.append([a, b, c, item[1]])

    # make triplets of facets
    triplets = [combo for combo in itertools.combinations(facet_equations, 3)]

    intersect_points = []
    for i, item in enumerate(triplets):
        item = np.array(item)
        rank = np.linalg.matrix_rank(item[:, :3])
        # check if planes intersect in a point and solve system
        if rank == 3:
            # the first term is multiplied by -1 to express the system as "ax + by + cz = d"
            point = np.linalg.solve(-1 * item[:, :3], item[:, -1])
            intersect_points.append(point)

    # reduce points to vertices and calculate
    morphology_vertices = find_vertices(intersect_points, facet_equations)
    hull = ConvexHull(morphology_vertices)
    volume = hull.volume
    surface = hull.area

    return volume, surface


def growth_rates_from_file(crystal, distance_file):
    """
    Read h k l indices and facet distances from file specified by user and return in
    right format for morphology_properties function.

    :param crystal: an instance of ccdc.crystal.Crystal
    :param distance_file: the file with the distance data (its path if not in working dir)
    :return: an iterable of pairs, ccdc.crystal.Crystal.MillerIndices and perpendicular distance
    """

    growth_rates = []
    with open(distance_file, 'r') as f:
        for line in f.readlines():
            # read only non-blank lines
            if line.strip():
                split = line.split()
                h, k, l = int(split[0]), int(split[1]), int(split[2])
                growth_rates.append((crystal.MillerIndices(h, k, l, crystal), float(split[3])))
    return growth_rates


def demo():
    """
    Function to test this script using the HXACAN growth_rates from
    CSD Python API documentation >> Descriptive documentation >> Morphology.
    """

    reader = EntryReader('CSD')
    entry = reader.entry('HXACAN')
    crystal = entry.crystal

    growth_rates = [
        (crystal.MillerIndices(2, 1, 1, crystal), 0.2004225755740747),
        (crystal.MillerIndices(-1, 1, -1, crystal), 0.19703805522607432),
        (crystal.MillerIndices(1, -1, -1, crystal), 0.19683866569762418),
        (crystal.MillerIndices(2, 1, -1, crystal), 0.20042073834713292),
        (crystal.MillerIndices(-2, 2, 0, crystal), 0.21988965080867295),
        (crystal.MillerIndices(-2, -1, -1, crystal), 0.2007281995185276),
        (crystal.MillerIndices(-2, 1, -1, crystal), 0.20072647044205774),
        (crystal.MillerIndices(2, -1, -1, crystal), 0.20041901523173083),
        (crystal.MillerIndices(-2, 1, 1, crystal), 0.20072463001602556),
        (crystal.MillerIndices(0, 2, 0, crystal), 0.23260015100470804),
        (crystal.MillerIndices(1, 1, -1, crystal), 0.1968395609520126),
        (crystal.MillerIndices(-2, 0, 1, crystal), 0.19266228580032901),
        (crystal.MillerIndices(2, 0, 1, crystal), 0.19235971851036732),
        (crystal.MillerIndices(-2, 0, 0, crystal), 0.1956121055387738),
        (crystal.MillerIndices(1, 2, 0, crystal), 0.2273675403933891),
        (crystal.MillerIndices(2, 0, 0, crystal), 0.1952180443198133),
        (crystal.MillerIndices(2, 0, -1, crystal), 0.19235690444211515),
        (crystal.MillerIndices(-2, 0, -1, crystal), 0.1926651034625492),
        (crystal.MillerIndices(-2, -2, 0, crystal), 0.21988995158315175),
        (crystal.MillerIndices(0, -2, 0, crystal), 0.23260015804584852),
        (crystal.MillerIndices(2, 2, 0, crystal), 0.21953027052221527),
        (crystal.MillerIndices(1, 1, 1, crystal), 0.19683926430195395),
        (crystal.MillerIndices(2, -1, 1, crystal), 0.20042085241706384),
        (crystal.MillerIndices(-1, 1, 1, crystal), 0.19703835382502216),
        (crystal.MillerIndices(-1, -1, 1, crystal), 0.19703925375799888),
        (crystal.MillerIndices(-1, 2, 0, crystal), 0.2276236144133509),
        (crystal.MillerIndices(1, -2, 0, crystal), 0.22736517051933225),
        (crystal.MillerIndices(-1, -2, 0, crystal), 0.2276260079936761),
        (crystal.MillerIndices(-2, -1, 1, crystal), 0.20072635905074157),
        (crystal.MillerIndices(1, -1, 1, crystal), 0.19683836895638027),
        (crystal.MillerIndices(2, -2, 0, crystal), 0.21952998690923803),
        (crystal.MillerIndices(-1, -1, -1, crystal), 0.1970389552502907)]

    # volume using this method
    morph_volume, morph_surface = morphology_properties(growth_rates)

    print("With distances from example")
    print(f"volume of the morphology: {morph_volume}")

    # do the same but double the input distances
    growth_rates = [[item[0], item[1] * 2] for item in growth_rates]

    # recalculate the volume
    morph_volume, morph_surface = morphology_properties(growth_rates)

    print("Now the input distances are doubled")
    print(f"volume of the morphology: {morph_volume}")

    return morph_volume, morph_surface


if __name__ == '__main__':
    morphology_volume, morphology_surface = demo()
