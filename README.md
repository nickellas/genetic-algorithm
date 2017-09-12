# genetic-algorithm
This directory contains the research I have been conducting since November 2016, I have assembled a previously constructed genetic algorithm and made changes to several aspects of the initial population step in an attemp to improve its speed and accuracy.

requirement to use this program:

*The ase package must be present on your system. It can be found here: https://wiki.fysik.dtu.dk/ase/install.html

*After ase has been installed clustering_v3.py and comparators_opt.py should be moved into the ase/data directory

*Several calculators can be used to optimaize and calculate energies. DFTB plus is the calculator used here but it is not recommended for production level work. DFT is recommended and the gpaw package should be used for production level work.
  https://wiki.fysik.dtu.dk/gpaw/

*If you wish to use Dftb plus you must install arpack on your system in addition to DFTB_plus.
  https://www.dftbplus.org/
  http://www.caam.rice.edu/software/ARPACK/
  
NOTE: These instructions worked on my schools cluster but any number of changes may be needed for them to work of your architecture. Follow all readme files on these programs.

___program guide___

The programs in this directory are designed to work together when genetic_algorithm.py is run. They each do the following:

genetic_algorithm.py (Nicholas Kellas)- The brains of the operation. To run, change the the stoichiometry desired (atom_numbers), a slab may be added if desired. The number of times the program will generate a new candidate is determined by n_to_test. By default the program
will run 5000 iterations with a 30% chance to mutate the created molecule rather then just cut-and-splice.

molemaker.py (Nicholas Kellas)- This module generate the molecules to form the initial population. It knows how to use Carbon, Hydrogen,
and Nitrogen currently, but can easily have other elements added.

Starting_pop_maker.py (Nicholas Kellas)- This program is a plug-in to use the two remaining programs to cluster the starting population.
With it, a large starting population can be generated without hugely time consuming calculations, and then the molecules can be grouped
based on intermolecular distances to create 20 distinct groups of molecules. A molecule is picked from each group at random to form the
true starting population.

clustering_v3.py (Mathias S Jorgensen)- Identifies similarity between molecules for grouping.

comparators_opt.py (Mathias S Jorgensen)- Chunky math stuff used to make the clustering program work.
