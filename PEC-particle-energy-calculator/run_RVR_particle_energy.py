"""
Run PEC to reproduce the results from Figure 2c
for the two polymorphs of Ritonavir (RVR), form 1 (RVR-I) and form 2 (RVR-II).

Particle energies are calculated for PED (Particle Equivalent Diameter) sizes between 20 nm and 100 nm for two cases:
    - growth morphologies of RVR-I and RVR-II
    - needle morphologies of RVR-I and RVR-II

"""

import os
import numpy as np
import time

from particle_energy_calculator import NanoParticle
from pec_visualiser import ParticleEnergyVisualiser
import pec_utilities


def run_size_sweep(crystal, lattice_energy, attachment_energies, user_distances=None):
    """ Calculate a nanoparticle for a range of sizes and return the NanoParticle objects """

    nanoparticles = []

    # although the final values are reported in nm, we initially create a Nanoparticle using distances in Angstroms!
    ped_values = np.linspace(200, 1000, 50)

    for size in ped_values:
        # if a morphology is specified, use that morphology
        if user_distances:
            nanoparticle = NanoParticle(crystal, attachment_energies, lattice_energy, size=size, size_type="diameter",
                                        user_distances=user_distances)
        else:
            # otherwise, use the attachment energy morphology
            nanoparticle = NanoParticle(crystal, attachment_energies, lattice_energy, size=size, size_type="diameter")
        nanoparticles.append(nanoparticle)

    return nanoparticles


if __name__ == "__main__":
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
                                                                                                os.getcwd(),
                                                                                                "input_files"))

    rvr_2_lattice_energy, rvr_2_attachment_energies = pec_utilities.get_attachment_energies(rvr_2_crystal,
                                                                                            filename="YIGPIO03_energies.txt",
                                                                                            filepath=os.path.join(
                                                                                                os.getcwd(),
                                                                                                "input_files"))

    """ Run using the attachment energy morphologies """
    # run the calculation for RVR-I using the growth morphology (attachment energy morphology)
    time_start = time.time()
    rvr_1_nanoparticles = run_size_sweep(rvr_1_crystal, rvr_1_lattice_energy, rvr_1_attachment_energies)
    print(f"time for {len(rvr_1_nanoparticles)} particles of RVR-I: {round(time.time() - time_start, 2)} s")

    # run the calculation for RVR-II using the growth morphology (attachment energy morphology)
    time_start = time.time()
    rvr_2_nanoparticles = run_size_sweep(rvr_2_crystal, rvr_2_lattice_energy, rvr_2_attachment_energies)
    print(f"time for {len(rvr_2_nanoparticles)} particles of RVR-II: {round(time.time() - time_start, 2)} s")

    print(verify_cross_size(rvr_1_nanoparticles, rvr_2_nanoparticles))

    # finally, plot the energies
    ParticleEnergyVisualiser.generate_energy_plot(rvr_1_nanoparticles, rvr_2_nanoparticles)

    """ Run using needle morphologies """
    # run the calculation for RVR-I using a needle morphology (morphology file number 326)
    rvr_1_needle_distances = pec_utilities.read_morphology_file(rvr_1_crystal,
                                                                os.path.join("morphologies", "form_1",
                                                                             "form_1_morph_326"))
    time_start = time.time()
    rvr_1_nanoparticles = run_size_sweep(rvr_1_crystal, rvr_1_lattice_energy, rvr_1_attachment_energies,
                                         user_distances=rvr_1_needle_distances)
    print(f"time for {len(rvr_1_nanoparticles)} particles of RVR-I: {round(time.time() - time_start, 2)} s")

    # run the calculation for RVR-II using a needle morphology (morphology file number 270)
    rvr_2_needle_distances = pec_utilities.read_morphology_file(rvr_2_crystal,
                                                                os.path.join("morphologies", "form_2",
                                                                             "form_2_morph_270"))
    time_start = time.time()
    rvr_2_nanoparticles = run_size_sweep(rvr_2_crystal, rvr_2_lattice_energy, rvr_2_attachment_energies,
                                         user_distances=rvr_2_needle_distances)
    print(f"time for {len(rvr_2_nanoparticles)} particles of RVR-II: {round(time.time() - time_start, 2)} s")

    print(pec_utilities.verify_cross_size(rvr_1_nanoparticles, rvr_2_nanoparticles))

    # finally, plot the energies
    ParticleEnergyVisualiser.generate_energy_plot(rvr_1_nanoparticles, rvr_2_nanoparticles)
