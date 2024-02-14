import os
import re
import numpy as np
import itertools
from ccdc.io import EntryReader, CrystalReader


def parse_input_file(filename, filepath):
    """
    Read the lines from an input file using a standard format.
    """

    def parse_line(line):
        pattern = r"[\(\[\{\}]\s*(.*?)\s*[\)\]\}]\s*(.*)"
        match = re.search(pattern, line)
        if match:
            string_with_brackets = "{" + match.group(1) + "}"
            float_value = float(match.group(2))
            return string_with_brackets, float_value
        else:
            raise ValueError("Can't read your data, try changing the format.")

    absolute_path = os.path.join(filepath, filename)

    hkl_strings = []
    attachment_energy_values = []
    with open(absolute_path, "r") as f:
        lines = f.readlines()
        title = lines[0].rstrip("\n")
        structure = lines[1].rstrip("\n")
        lattice_energy = float(lines[3].rstrip("\n"))
        facet_lines = lines[4:]
        for line in facet_lines:
            split = parse_line(line.rstrip("\n"))
            hkl_strings.append(split[0])
            attachment_energy_values.append((split[1]))

    return lattice_energy, hkl_strings, attachment_energy_values


def read_structure(structure_name, path=os.getcwd()) -> object:
    """Read a crystal structure from a CIF file or from the CSD"""

    if "." not in structure_name:
        reader = EntryReader()
        crystal = reader.entry(structure_name).crystal

    else:
        cif_file = os.path.join(path, structure_name)
        crystal_reader = CrystalReader(cif_file)
        crystal = crystal_reader[0]

    return crystal


def morphology_vertices(morphology):
    """
    Returns the coordinates of the unique vertices in the morphology

    :param morphology: a ccdc.morphology Morphology object
    """

    coordinates = []
    for f in morphology.facets:
        for coord in f.coordinates:
            # values are rounded so they can be compared
            x, y, z = round(coord.x, 6), round(coord.y, 6), round(coord.z, 6)
            coordinates.append([x, y, z])

    # reduce to unique coordinates
    morph_vertices = [list(x) for x in set(tuple(x) for x in coordinates)]

    return morph_vertices


def distance(p1, p2) -> float:
    """Euclidean distance between two points"""
    dist = np.sqrt(((p2[0] - p1[0]) ** 2) + ((p2[1] - p1[1]) ** 2) + ((p2[2] - p1[2]) ** 2))
    return dist


def vector_from_points(p1, p2):
    """Vector between two points from p1 to p2"""
    return p2 - p1


def vector_length(vector: np.ndarray):
    return np.linalg.norm(vector)


def orthogonal_box_aspect_ratio(vertices):
    """
    Calculate the aspect ratio of a morphology given its vertices.
    This aspect ratio is not a minimum volume (oriented bounding box) aspect ratio,
    but just the ratio of the dimensions of the orthogonal box aligned with the principal axes.
    """

    box_vert = [max(point[0] for point in vertices), max(point[1] for point in vertices),
                max(point[2] for point in vertices)]

    default_vertices = np.array([
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [-1, -1, +1],
        [+1, -1, +1],
        [+1, +1, +1],
        [-1, +1, +1],
    ])
    vertices = box_vert * default_vertices
    vert_pairs = [x for x in itertools.combinations(vertices, 2)]
    x_sides, y_sides, z_sides = [], [], []
    for pair in vert_pairs:
        p1, p2 = np.array(pair[0]), np.array(pair[1])
        vector = vector_from_points(p1, p2)
        if vector[1] == 0 and vector[2] == 0:
            x_sides.append(vector_length(vector))
        if vector[0] == 0 and vector[2] == 0:
            y_sides.append(vector_length(vector))
        if vector[0] == 0 and vector[1] == 0:
            z_sides.append(vector_length(vector))
    x, y, z = np.mean(x_sides), np.mean(y_sides), np.mean(z_sides)
    ratios = [x / z, y / z, z / z]
    return sorted(ratios, reverse=True)


def calculate_aspect_ratio(morphology):
    """Calculate the aspect ratio of the morphology"""
    return orthogonal_box_aspect_ratio(morphology_vertices(morphology))


def relative_facet_distances(morphology):
    """Returns the relative distances of the morphology"""

    # These distances are scaled according to the CSD Python API method
    # where the smallest distance is assigned a value of 100
    relative_distances = []
    for f in morphology.facets:
        miller = f.miller_indices
        facet_distance = f.perpendicular_distance
        relative_distances.append((miller, facet_distance))
    return relative_distances


def sym(p1, operator, crystal, fractional_coordinates=False):
    """
    Apply symmetry operator on a point with fractional coordinates.
    If coordinates are cartesian, they are first converted and then converted back.

    :param p1: the point coordinates
    :param operator: the symmetry operator -> string
    :param crystal: Crystal object
    :param fractional_coordinates: boolean for supplied coordinates
    :return:
    """
    opx, opy, opz = operator.split(',')[0], operator.split(',')[1], operator.split(',')[2]

    if not fractional_coordinates:
        # convert to fractional coordinates
        p1 = cart2frac(p1, crystal)

    # following line is needed so that builtin eval can work
    x, y, z = p1[0], p1[1], p1[2]
    # apply symmetry operation
    p2 = [eval(opx), eval(opy), eval(opz)]

    if not fractional_coordinates:
        # convert back to cartesian coordinates
        p2 = frac2cart(p2, crystal)

    return p2


def sym_strip(operator_string) -> str:
    """
    Read the symmetry operator string and retain only the rotation/reflection part in case of glide planes and/or screw
    axes.
    This was written with triclinic, monoclinic and orthorombic SGs in mind.
    Will need to modify for less common symmetry operators and space groups.
    """

    # remove the 1/2 fraction for glide planes and screw axes from string, we only care about rotations and reflections
    filtered = filter(lambda ch: ch not in "1/2", operator_string)
    new_string = ""
    for x in filtered:
        new_string += x
    return new_string


def index_strip(index_string) -> list:
    """
    Read hkl indices from a string written in various formats and return as list of integers.
    Will not work for indices > abs(9).

    :param index_string: a string with the hkl indices of a facet or family of facets e.g. "( -1 0 0)" or "{-1 0 0}" etc
    :return: a list of integers representing the hkl indices
    """

    filtered = filter(lambda ch: ch not in ",{}() ", index_string)
    values = []
    integers = []
    new_string = ""
    for x in filtered:
        new_string += x
    check = 0
    for i, char in enumerate(new_string):
        if char == "-":
            # if minus sign then assign negative
            values.append(-1)
            check = 1
            continue
        if char != "-":
            integers.append(int(char))
            if check != 1:
                values.append(1)
        check = 0

    return [values[i] * integers[i] for i in range(len(integers))]


def apply_symmetry_hkl(crystal, h, k, l):
    """
    Apply symmetry to the facet indices to calculate symmetrical facets
    """
    symm_ops = crystal.symmetry_operators

    unique_indices = []

    for op in symm_ops:

        op_new = sym_strip(op)
        value = sym((h, k, l), op_new, crystal, fractional_coordinates=True)
        if value not in unique_indices:
            unique_indices.append(value)

    return unique_indices


def assign_growth_rates_from_file(crystal, hkl_indices, attachment_energy_values):
    """
    Returns the growth rates in the right format for volume calculation after calculating symmetry equivalent facets.

    :param crystal: a ccdc.Crystal object
    :param hkl_indices: a list of strings for the hkl indices, can contain () and {}, as well as spaces
    :param attachment_energy_values: a list of attachment energies (negative float values)
    :return: an iterable of (Crystal.MillerIndices, attachment energy) for each symmetry independent facet
    """

    growth_rates = []

    for i in range(len(hkl_indices)):
        hkl = index_strip(hkl_indices[i])
        # generate symmetry-equivalent facets and assign distances
        sym_indices = apply_symmetry_hkl(crystal, hkl[0], hkl[1], hkl[2])
        for item in sym_indices:
            growth_rates.append((crystal.miller_indices(item[0], item[1], item[2]), attachment_energy_values[i]))

    return growth_rates


def get_attachment_energies(crystal, filename=None, filepath=None):
    """
    Get the attachment energies for your structure.

    :param crystal: a Crystal object
    :param filename: the name of the file with the user energies (could be combined with the path)
    :param filepath: the path of the file with the user energies
    :return: the lattice energy and the attachment energies
    """

    lattice_energy, hkl_indices, attachment_energy_values = parse_input_file(filename,
                                                                             filepath)
    attachment_energies = assign_growth_rates_from_file(crystal, hkl_indices,
                                                        attachment_energy_values)

    return lattice_energy, attachment_energies


def read_morphology_file(crystal, morphology_file):
    """
    Reads a morphology file and returns the perpendicular distances in the right format

    :param crystal: a ccdc.Crystal object
    :param morphology_file: a filename or filepath containing only text. Each line has h, k, l, d values.
    """
    morphology_distances = []
    with open(morphology_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        split = line.split(" ")
        h, k, l = int(split[0]), int(split[1]), int(split[2])
        d = float(split[3])
        morphology_distances.append((crystal.miller_indices(h, k, l), d))
    return morphology_distances
