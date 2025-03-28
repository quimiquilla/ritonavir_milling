# Particle Energy Calculator (PEC)

A Python code to calculate the total energy of idealized nanosized crystal particles.

For more details about the methodology
see [Crystal size, shape, and conformational changes drive both the disappearance and reappearance of ritonavir polymorphs in the mill](https://www.pnas.org/doi/abs/10.1073/pnas.2319127121).
Please cite this paper if you use PEC for your work.

## Dependencies

- [CSD Python API](https://www.ccdc.cam.ac.uk/solutions/software/csd-python/)
- Scipy
- itertools

## Contents

The principal scripts used by PEC are:

- `particle_energy_calculator.py`

- `morphology_volume.py`

- `pec_utilities`

- `pec_visualiser`

Additionally, the results presented in Figure 2
of [Crystal size, shape, and conformational changes drive both the disappearance and reappearance of ritonavir polymorphs in the mill](https://www.pnas.org/doi/abs/10.1073/pnas.2319127121)
can be reproduced by running the following scripts:

- Figure 2b: `run_RVR_aspect_ratio.py`
- Figure 2c: `run_RVR_particle_energy.py`
- Figure 2d: `run_RVR_cross_check.py`

The subdirectory **input_files** contains the *.cif* files of the DFT-d optimised crystal structures of Ritonavir form I and form II, *YIGPIO02* and *YIGPIO03*, respectively, as well as their corresponding energies.
The subdirectory **morphologies** contains the various crystal morphologies used in our work, which hare provided as a list of `h k l d` values, where *h*, *k* and *l* are the miller indices, and *d* is the perpendicular distance.

## Usage

This code requires the definition of a crystal structure by reading a `.cif` file or by reading
a [CSD](https://www.ccdc.cam.ac.uk/solutions/software/csd/) refcode.

Lattice energy and attachment energies of hkl planes are also required. These can be calculated with any preferred 
method. In the examples presented here, energies are read from input `.txt` files with the following format:

- line [1]: a title
- line [2]: the name of the `.cif` file (this doesn't need to be specified here, but it is useful to keep track of things)
- line [3]: further comments
- line [4]: the lattice energy
- lines [5:]: the {hkl} indices of the form and the relative attachment energy


```
RVR Form I
YIGPIO02_opt_sym.cif
user_input
-395.0
{  0  0  1}	-54.4
{  1  0  0}	-120.0
{  1  0 -1}	-121.6
{  0  1  1}	-326.5
{  0 -1  1}	-326.5
{  1  1  0}	-368.0
{  1 -1  0}	-368.0
{  1  1 -1}	-414.6
{  1 -1 -1}	-414.6
```

Have a look at `run_RVR_particle_energy.py` to see an example of a calculation with PEC.

## Author

Pietro Sacchi
