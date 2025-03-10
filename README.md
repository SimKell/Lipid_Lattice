Lipid Lattice
==============================

This program performs simple Monte-Carlo simulations of a lattice system of lipids and CHOL.
Based on interaction functions derived by MD simulations performed by Fabian Keller.
The first version of a squared lattice system for pure lipid mixtures was implementend by Davit Hokabyan & Andreas Heuer.
Example Input and Output Files can be found in "/Check" for comparison after changing the code.

File Structure
-------------------------------
* **input.txt**: Input for the simulations. Use only after compiling ``Lipid_Lattice.c``.
* **Kite_Plot.py**: Creates Kite Plots of the order parameter distributions.
* **Lattice_Plot.py**: Creates spatial maps of the lattice.
* **Lipid_Lattice.c**: The main code, responsible for the Monte-Carlo simulations.
* **run.sh**: A batch script, changing the input file for the simulation of different temperatures and CHOL-concentrations (and running the simulations).

Changelog
-------------------------------
* 0.0.1: 
   * Initial Release
* 0.0.2: 
   * Changed Output from ``.dat`` to ``.csv``
   * Added Example Files in ``/Check``
* 0.0.3:
   * Changed Function ``Phi_P`` in ``Lipid_Lattice.c`` (see there)
* 0.0.4:
   * Changed Function ``mc_move_order`` in ``Lipid_Lattice.c`` (see there)
