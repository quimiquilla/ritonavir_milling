# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.

"""

Particle Energy Calculator (PEC)

Main class to calculate particle energy.

"""

import numpy as np
import morphology_volume
from ccdc.morphology import MorphologyBase


class NanoParticle:
    """Creates a NanoParticle object representing the idealised nano-sized crystal."""

    def __init__(self, crystal, attachment_energies, lattice_energy, size, size_type="diameter", user_distances=None,
                 surface_multiplier=1):
        """
        :param obj crystal: a ccdc.Crystal object
        :param list attachment_energies: a list of (ccdc.Crystal.MillerIndices, attachment energy)
        :param float lattice_energy: the crystal lattice energy
        :param float size: a size in Angstrom, Angstrom^2 or Angstrom^3, depending on the selected size_type
        :param str size_type: use to control the particle size, can be "diameter", "volume" or "surface"
        :param list user_distances: optional, a list of (ccdc.Crystal.MillerIndices, perpendicular distances)
        :param int surface_multiplier: an integer to determine the thickness of the particle's surface.
                Higher values result in larger surface energy penalties.
        """
        self.crystal = crystal
        self.identifier = self.crystal.identifier
        self.bulk_energy = lattice_energy
        self.attachment_energies = attachment_energies
        self._size_type = size_type
        self._surface_multiplier = surface_multiplier

        if size_type == "diameter":
            self.size = size
            self.equivalent_diameter = size
        if size_type == "surface":
            self.size = size
            self.equivalent_diameter = self.surface_to_diameter(surface=size)
        if size_type == "volume":
            self.size = size
            self.equivalent_diameter = self.volume_to_diameter(volume=size)

        if user_distances:
            self._growth_rates = user_distances
        else:
            self._growth_rates = attachment_energies

        self.total_surface_molecules = None

        self._volume = None
        self._surface = None
        self._inner_hull_surface = None
        self._particle_energy = None
        self._surface_energy_penalty = None
        self._scaled_facet_distances = self._scale_facet_distances()
        self._inner_hull_volume = self.calculate_morphology_properties(self._inner_hull_distances)[0]

    @staticmethod
    def diameter_to_volume(diameter):
        """Converts a diameter to corresponding sphere volume"""
        return np.pi * diameter ** 3 / 6

    @staticmethod
    def volume_to_diameter(volume):
        """Converts a volume to corresponding sphere diameter"""
        return ((volume * 6) / np.pi) ** (1 / 3)

    @staticmethod
    def surface_to_diameter(surface):
        """Converts a surface area to corresponding sphere diameter"""
        return (surface / np.pi) ** (1 / 2)

    @staticmethod
    def calculate_morphology_properties(growth_rates):
        """Calculate the morphology properties after scaling to the desired size"""
        volume, surface_area = morphology_volume.morphology_properties(growth_rates)
        return volume, surface_area

    def _scale_facet_distances(self):
        """Scale the perpendicular facet distances to attain the specified particle size"""
        # start by using the relative distances to get scale factors
        volume_for_scaling, surface_for_scaling = morphology_volume.morphology_properties(self._relative_distances)
        if self._size_type == "volume":
            scale_factor = (self.size / volume_for_scaling) ** (1 / 3)
            scaled_distances = [(item[0], item[1] * scale_factor) for item in self._relative_distances]
        if self._size_type == "surface":
            scale_factor = (self.size / surface_for_scaling) ** (1 / 2)
            scaled_distances = [(item[0], item[1] * scale_factor) for item in self._relative_distances]
        if self._size_type == "diameter":
            # this is basically scaled like the volume
            desired_volume = self.diameter_to_volume(self.size)
            scale_factor = (desired_volume / volume_for_scaling) ** (1 / 3)
            scaled_distances = [(item[0], item[1] * scale_factor) for item in self._relative_distances]

        return scaled_distances

    def _particle_properties(self):
        """Get the volume and surface area of the particle"""
        self._volume, self._surface = self.calculate_morphology_properties(self._scaled_facet_distances)

    def _calc_particle_energy(self):
        """Calculate the particle energy"""
        self._particle_energy = self.bulk_energy + self.surface_energy_penalty

    def _calc_surface_energy_penalty(self):
        """Calculate the total surface energy penalty of the particle"""
        total_energy = 0
        for miller_indices, energy in self.attachment_energies:
            morphological_importance = self.morphology.relative_area(miller_indices)
            total_energy += self.surface_fraction * morphological_importance * energy
        self._surface_energy_penalty = -1 * total_energy / 2
        return self._surface_energy_penalty

    @property
    def _relative_distances(self):
        """Calculate the facets' relative distances"""
        if self._growth_rates[0][1] < 0:
            # if the input data are attachment energies (negative values), we use the least negative (max) as reference
            reference_value = max(self._growth_rates, key=lambda x: x[1])[1]
        else:
            # otherwise we use the minimum distance
            reference_value = min(self._growth_rates, key=lambda x: x[1])[1]
        relative_d = []
        for item in self._growth_rates:
            relative_d.append((item[0], item[1] / reference_value))
        return relative_d

    @property
    def morphology(self):
        return MorphologyBase.from_growth_rates(self.crystal, self._growth_rates)

    @property
    def total_particle_molecules(self):
        """Calculate the number of molecules in the particle"""
        return self.crystal.z_value * self.volume / self.crystal.cell_volume

    @property
    def volume(self):
        """Returns volume of scaled particle"""
        if self._volume:
            return self._volume
        else:
            self._particle_properties()
            return self._volume

    @property
    def total_surface_area(self):
        """Returns total surface area of scaled particle"""
        if self._surface:
            return self._surface
        else:
            self._particle_properties()
            return self._surface

    @property
    def particle_mass(self):
        """ Returns the crystal mass in grams"""
        # convert from Angstrom^3 to cm^3
        return self.volume * self.crystal.calculated_density * 1E-24

    @property
    def smallest_dimension(self):
        """ Returns the length of the smallest side of the particle bounding box """
        # size will be in Angstrom
        major, median, minor = self.morphology.oriented_bounding_box.lengths
        return (self.volume * minor ** 3 / self.morphology.volume) ** (1 / 3)

    @property
    def _inner_hull_distances(self):
        """Get the distances used to define the inner hull"""
        distances = []
        for i, item in enumerate(self._scaled_facet_distances):
            facet_distance = item[1]
            d_hkl = item[0].d_spacing
            hull_distance = facet_distance - self._surface_multiplier * d_hkl
            if facet_distance <= d_hkl:
                raise ValueError("The particle size you are using is too small! Try using a larger size.")
            distances.append((item[0], hull_distance))
        return distances

    @property
    def surface_fraction(self):
        """Calculate the total fraction of surface molecules"""
        # see https://pubs.acs.org/doi/10.1021/acs.jchemed.0c01247
        return (self.volume - self._inner_hull_volume) / self.volume

    @property
    def particle_energy(self):
        """Calculate the particle energy"""
        if self._particle_energy:
            return self._particle_energy
        else:
            self._calc_particle_energy()
            return self._particle_energy

    @property
    def surface_energy_penalty(self):
        """Return the total surface energy penalty of the nano-sized crystal"""
        if self._surface_energy_penalty:
            return self._surface_energy_penalty
        else:
            self._calc_surface_energy_penalty()
            return self._surface_energy_penalty
