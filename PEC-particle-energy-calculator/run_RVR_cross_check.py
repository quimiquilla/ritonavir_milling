"""
Run PEC to reproduce the results from Figure 2d
for the two polymorphs of Ritonavir (RVR), form 1 (RVR-I) and form 2 (RVR-II).

Particle energies are calculated for PED (Particle Equivalent Diameter) sizes between 20 nm and 100 nm.
"""

import os
import time

from particle_energy_calculator import NanoParticle
from pec_visualiser import CrossCheckVisualiser
import pec_utilities


def run_morphology_sweep(crystal, lattice_energy, attachment_energies, morphology_distances, equivalent_diameter=600):
    """ Calculate nanoparticles for a range of morphologies and a fixed size and return the NanoParticle objects """
    time_start = time.time()
    nanoparticles = []
    for user_distances in morphology_distances:
        # if the size is too small the PEC code will fail
        try:
            nanoparticle = NanoParticle(crystal, attachment_energies, lattice_energy, size=equivalent_diameter,
                                        size_type="diameter", user_distances=user_distances)
            nanoparticles.append(nanoparticle)
        except ValueError:
            continue
    print(f"time for {len(nanoparticles)} nanoparticles: {round(time.time() - time_start, 2)} s")
    return nanoparticles


def run_cross_check(nanoparticles_1, nanoparticles_2, energy_difference):
    """ Check all combinations of particles to see which result in switches for the given size"""
    count_total = 0
    count_switched = 0
    for part_1 in nanoparticles_1:
        for part_2 in nanoparticles_2:
            # the surface energy difference (changed of sign) needs to be larger than lattice energy difference
            en_diff = - (part_1.surface_energy_penalty - part_2.surface_energy_penalty)
            if en_diff >= energy_difference:
                count_switched += 1
            count_total += 1

    return count_switched, count_total


def read_morphologies(crystal, path):
    morphology_distances = []

    for f in os.listdir(path):
        morphology_distances.append(pec_utilities.read_morphology_file(crystal, os.path.join(path, f)))

    return morphology_distances


# read the .cif (DFT optimized structure) of RVR-I (CSD refcode YIGPIO02)
rvr_1_crystal = pec_utilities.read_structure(structure_name="YIGPIO02_opt_sym.cif",
                                             path=os.path.join(os.getcwd(), "input_files"))

# read the .cif (DFT optimized structure) of RVR-II (CSD refcode YIGPIO03)
rvr_2_crystal = pec_utilities.read_structure(structure_name="YIGPIO03_opt_sym.cif",
                                             path=os.path.join(os.getcwd(), "input_files"))

# get the energies of RVR-I and RVR-II
rvr_1_lattice_energy, rvr_1_attachment_energies = pec_utilities.get_attachment_energies(rvr_1_crystal,
                                                                                        filename="YIGPIO02_energies.txt",
                                                                                        filepath=os.path.join(
                                                                                            os.getcwd(), "input_files"))

rvr_2_lattice_energy, rvr_2_attachment_energies = pec_utilities.get_attachment_energies(rvr_2_crystal,
                                                                                        filename="YIGPIO03_energies.txt",
                                                                                        filepath=os.path.join(
                                                                                            os.getcwd(), "input_files"))

# we calculate the lattice energy difference to check if switch condition is satisfied
lattice_energy_difference = rvr_1_lattice_energy - rvr_2_lattice_energy

# get the morphologies
rvr_1_morphologies = read_morphologies(rvr_1_crystal, os.path.join("morphologies", "form_1"))
rvr_2_morphologies = read_morphologies(rvr_2_crystal, os.path.join("morphologies", "form_2"))

switch_fraction = []
# the size in nanometers
sizes = [20, 30, 40, 50, 60, 70, 80, 90, 100]

time_start_overall = time.time()
# then we loop
for size in sizes:
    # NanoParticle wants the size in Angstroms
    diameter = size * 10
    print(f"### Running calculation for size of {size} nm ###")
    nanoparticles_1 = run_morphology_sweep(rvr_1_crystal, rvr_1_lattice_energy, rvr_1_attachment_energies,
                                           rvr_1_morphologies, equivalent_diameter=diameter)
    nanoparticles_2 = run_morphology_sweep(rvr_2_crystal, rvr_2_lattice_energy, rvr_2_attachment_energies,
                                           rvr_2_morphologies, equivalent_diameter=diameter)
    time_start = time.time()
    count_switched, count_total = run_cross_check(nanoparticles_1, nanoparticles_2, lattice_energy_difference)
    switch_fraction.append(100 * count_switched / count_total)
    print(f"time to calculate {count_total} combinations: {round(time.time() - time_start, 2)}")

print(f"total time: {round(time.time() - time_start_overall, 2)} s")

# finally plot the results
CrossCheckVisualiser.generate_cross_check_plot(sizes, switch_fraction)
