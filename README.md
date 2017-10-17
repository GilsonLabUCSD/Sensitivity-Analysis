# Sensitivity Analysis #

### Written by Jian (Jane) Yin (2017.10) ###

This code computes the analytical first-order derivatives of the standard chemical potential with respect to van der Waals parameters. To compute the derivative of the binding free 
energy, one may obtain the derivatives of the standard chemical potential from enthalpy calculations, i.e. simulations involving the initial and final states of the 
binding free energy calculations. For the detailed theoretical framework and applications of Sensitivity Analysis, please refer to the following publication:

Yin, J., Fenley, A. T., Henriksen, N. M., & Gilson, M. K. (2015). Toward improved force-field accuracy through sensitivity analysis of host-guest binding thermodynamics. 
The Journal of Physical Chemistry B, 119(32), 10145-10155. http://pubs.acs.org/doi/abs/10.1021/acs.jpcb.5b04262 

To run the Sensitivity Analysis scripts, you need to install Amber and parmEd.

## Usage ##

First of all, use the Python script provided here to generate molecular information (force field parameters and connectivity i.e. bonds, angles and torsions)  of your system
from the topology file.

For example:

	python parseTopoWithParmed.py full.topo

Only the Amber prmtop format is supported for the topology file (full.topo in this example).

Secondly, run the sensitivity executable with three flags: -crd, -cutoff and -gpu. 

Flags:

	-crd       MD trajectory file
	-cutoff    cutoff distance of nonbonded interactions
	-gpu       yes or no

Example:

	./sensitivity -crd traj.mdcrd -cutoff 9.0 -gpu yes

Only the Amber mdcrd format is supported for the coordinate file. If your trajectory uses the NetCDF or PDB format, you can use cpptraj or other programs to convert it to mdcrd. 
Enabling the -gpu option will increase the performance of the program significantly. 

The cutoff distance should match the abrupt cutoff of van ver Waals interaction used in simulations. The sensitivity code does not take into account the long-range corrections 
of the van der Waals interactions possibly used in MD simulations. However, it does not influence the predicting power of derivatives significantly. 

The derivatives of each frame and the mean values will be printed in radDerivatives.dat (derivatives to radius) and epsDerivatives.dat (derivatives to epsilon).
Note that in GAFF radius is used rather than sigma or R_min. So if you are interested in the derivatives to sigma or R_min, the derivatives to radius need to be converted accordingly.

If you would like to modify the script, issue

	make -f makefile
 
to compile the modified code. 

## Convergence ##

It is worth noting that the derivatives may not coverge for protein-ligand systems. We encourage you to write your own scripts to pull values from radDerivatives.dat
 and epsDerivatives.dat, and compute the uncertainties.   
