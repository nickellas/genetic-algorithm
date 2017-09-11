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
