from ase import Atoms, Atom
from random import uniform, choice, shuffle
from ase.io import write
from ase.data import covalent_radii, atomic_numbers
from numpy import pi, sin, cos
from atomplacer import molmaker, randomizer
from pickle import dump, load
from ase.calculators.emt import EMT


atoms = ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]

#atoms = []
#molstart = raw_input("what atom's are you testing. Enter them in the form 'NENE' where N is the number of atoms of one element followed by the element, denoted E:\n")
#molmaker(molstart, atoms)

#atcount = raw_input("How big do you want your population?")

moleinfo = {}

atcount = 10
while atcount > 0:
    molecule = Atoms(calculator = EMT())
    randomizer(atoms, molecule, moleinfo)
    atcount -= 1

#biginfo = moleinfo

with open('moleculedata.dat', 'rb') as f:
    biginfo = load(f)
    f.close()

biginfo.update(moleinfo)

with open("moleculedata.dat", "wb") as f:
    dump(biginfo, f)
    f.close()

print moleinfo
