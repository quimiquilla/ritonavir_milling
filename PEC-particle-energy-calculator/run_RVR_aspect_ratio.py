"""
Run PEC to reproduce the results from Figure 2b
for the two polymorphs of Ritonavir (RVR), form 1 (RVR-I) and form 2 (RVR-II).

Particle energies are calculated for PED (Particle Equivalent Diameter) size of 60 nm (can be changed).
"""

import os
import time

from particle_energy_calculator import NanoParticle
from pec_visualiser import AspectRatioVisualiser
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
    return nanoparticles


def read_morphologies(crystal, path):
    morphology_distances = []

    for f in os.listdir(path):
        morphology_distances.append(pec_utilities.read_morphology_file(crystal, os.path.join(path, f)))

    return morphology_distances


def plot_results(nanoparticles_1, nanoparticles_2):
    morphologies_1 = [nanoparticle.morphology for nanoparticle in nanoparticles_1]
    surface_penalties_1 = [nanoparticle.surface_energy_penalty for nanoparticle in nanoparticles_1]
    morphologies_2 = [nanoparticle.morphology for nanoparticle in nanoparticles_2]
    surface_penalties_2 = [nanoparticle.surface_energy_penalty for nanoparticle in nanoparticles_2]
    AspectRatioVisualiser.generate_AR_plot(morphologies_1, surface_penalties_1)
    AspectRatioVisualiser.generate_AR_plot(morphologies_2, surface_penalties_2)


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

# if needed, change the particle size
equivalent_diameter = 600
print(f"running calculation for particles with equivalent diameter of {equivalent_diameter / 10} nm")

# run for RVR-I
time_start = time.time()
rvr_1_morphologies = read_morphologies(rvr_1_crystal, os.path.join("morphologies", "form_1"))
rvr_1_nanoparticles = run_morphology_sweep(rvr_1_crystal, rvr_1_lattice_energy, rvr_1_attachment_energies,
                                           rvr_1_morphologies, equivalent_diameter=equivalent_diameter)
print(f"time for {len(rvr_1_nanoparticles)} nanoparticles of RVR-I: {round(time.time() - time_start, 2)} s")

# run for RVR-II
time_start = time.time()
rvr_2_morphologies = read_morphologies(rvr_2_crystal, os.path.join("morphologies", "form_2"))
rvr_2_nanoparticles = run_morphology_sweep(rvr_2_crystal, rvr_2_lattice_energy, rvr_2_attachment_energies,
                                           rvr_2_morphologies, equivalent_diameter=equivalent_diameter)
print(f"time for {len(rvr_2_nanoparticles)} nanoparticles of RVR-II: {round(time.time() - time_start, 2)} s")

# finally, plot the results
plot_results(nanoparticles_1=rvr_1_nanoparticles, nanoparticles_2=rvr_2_nanoparticles)
