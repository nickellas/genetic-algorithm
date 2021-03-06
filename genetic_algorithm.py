from ase.ga.data import PrepareDB, DataConnection
from molemaker import Starting_population
from random import random
from ase.io import write
from ase.optimize import BFGS
from ase.calculators.dftb import Dftb
from ase import Atoms
import numpy as np
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
from ase.constraints import FixAtoms
from starting_pop_maker import pop_maker
from time import time
import platform
from ase.visualize import view

# This program can be found on https://wiki.fysik.dtu.dk/ase/tutorials/ga/ga_optimize.html, modifications have been made by 
# Nicholas Kellas to improve how the staring population is created. The calculator used in this program, Dftb plus, does not come 
# installed with ase and must be installed manually along with arpack. Most other calculators may be used instead and for production
# level work DFT is recommended. Dftb plus is a fast calculator but is known to be inaccurate, it was used to test the improvement in 
# speed provided by alterations to the initial population phase of the GA and so accuracy was not required.

db_file = 'gadb.db'

# create the surface, In our testing we did not use a fixed surface but it was left in to allow quick addition if needed.
slab = Atoms('', positions = np.zeros((0,3)), cell = [30., 30., 30.])
slab.set_constraint(FixAtoms(mask=len(slab) * [True]))

# define the volume in which the adsorbed cluster is optimized
# the volume is defined by a corner position (p0)
# and three spanning vectors (v1, v2, v3)
pos = slab.get_positions()
cell = slab.get_cell()
p0 = np.array([0., 0., 12]) #max(pos[:, 2]) + 2.])
v1 = cell[0, :] * 0.8
v2 = cell[1, :] * 0.8
v3 = cell[2, :]
v3[2] = 3.

# Input the stoichiometry to be analyzed here the first number is the amount of the element and the bracketed number is the atomic number
# The stoichiometry below is for C9H7N
atom_numbers = 9 * [6] + 7 * [1] + 1 * [7]

# specifies how many molecules to generate for the initial population. This step is very fast so a large number is manageable if
# clustering will be used.
atcount = 500

start = time()
moleculegroup=[]
ga = Starting_population(atom_numbers, moleculegroup, cell)
ga.randomizer()
while atcount > len(moleculegroup):
    ga.randomizer()

write('candidates.traj', moleculegroup)
end = time()
print 'time to make candidate population', end-start

# pop_maker takes in the traj file made above and from the 500 candidates identifies 20 that are significantly different to use in the GA
start = time()
starting_pop = pop_maker()
end = time()
print 'time to make starting population', end-start

# creates the database file to hold the information generated
d = PrepareDB(db_file_name=db_file,
              simulation_cell = slab,
              stoichiometry=atom_numbers)

# The 20 molecules identified above become the starting population for the GA
for i in starting_pop:
    d.add_unrelaxed_candidate(i)

# n_to_test = the number of times the GA will generate a new molecule.
population_size = len(starting_pop)
mutation_probability = 0.3
n_to_test = 5000

da = DataConnection('gadb.db')
atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
blmin = closest_distances_generator(all_atom_types,
                                    ratio_of_covalent_radii=0.7)

comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                     pair_cor_cum_diff = 0.015,
                                     pair_cor_max=0.7,
                                     dE=0.02,
                                     mic=False)

pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
mutations = OperationSelector([1., 1., 1.],
                              [MirrorMutation(blmin, n_to_optimize),
                               RattleMutation(blmin, n_to_optimize),
                               PermutationMutation(n_to_optimize)])

# The starting population has its energy calculated and is optimized to a local minima or until it has tried for 100 iterations.
# Dftb_plus and BFGS can be substituted for other calculators and optimizers.
while da.get_number_of_unrelaxed_candidates() > 0:
    start_relax = time()
    a = da.get_an_unrelaxed_candidate()
    a.set_calculator(Dftb(label='C9H7N',
                          Hamiltonian_MaxAngularMomentum_='',
                          Hamiltonian_MaxAngularMomentum_C='"p"',
                          Hamiltonian_MaxAngularMomentum_H='"s"',
                          Hamiltonian_MaxAngularMomentum_N='"p"'
                          ))
    print('Relaxing starting candidate {0}'.format(a.info['confid']))
    dyn = BFGS(a, trajectory=None, logfile=None)
    dyn.run(fmax=0.05, steps=100)
    a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()
    da.add_relaxed_step(a)
    end_relax = time()
    print 'time to optimize:', end_relax-start_relax

population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

# test n_to_test new candidates
for i in range(n_to_test):
    print('Now starting configuration number {0}'.format(i))
    a1, a2 = population.get_two_candidates()
    a3, desc = pairing.get_new_individual([a1, a2])
    if a3 is None:
        continue
    da.add_unrelaxed_candidate(a3, description=desc)

    # Check if we want to do a mutation
    if random() < mutation_probability:
        a3_mut, desc = mutations.get_new_individual([a3])
        if a3_mut is not None:
            da.add_unrelaxed_step(a3_mut, desc)
            a3 = a3_mut

    # Relax the new candidate
    start_relax = time()
    a3.set_calculator(Dftb(label='C9H7N',
                           Hamiltonian_MaxAngularMomentum_='',
                           Hamiltonian_MaxAngularMomentum_C='"p"',
                           Hamiltonian_MaxAngularMomentum_H='"s"',
                           Hamiltonian_MaxAngularMomentum_N='"p"'
                           ))

    dyn = BFGS(a3, trajectory=None, logfile=None)
    dyn.run(fmax=0.05, steps=100)
    a3.info['key_value_pairs']['raw_score'] = -a3.get_potential_energy()
    da.add_relaxed_step(a3)
    population.update()
    end_relax = time()
    print 'time to optimize:', end_relax-start_relax
    # used in testing to identify if the known global minima had been found (Quinoline)
    #if a3.get_potential_energy() < -557.783851467:
        #break

write('all_candidates.traj', da.get_all_relaxed_candidates())
