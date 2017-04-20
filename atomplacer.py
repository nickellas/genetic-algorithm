from ase import Atoms, Atom
from random import uniform, choice, shuffle
from ase.io import write
from ase.data import covalent_radii, atomic_numbers
from numpy import pi, sin, cos
from ase.calculators.emt import EMT

#atoms = []
#molstart = raw_input("what atom's are you testing. Enter them in the form 'NENE' where N is the number of atoms of one element followed by the element, denoted E:\n")

def molmaker(molstart, atoms):
    molstart += "E"
    atomcount = 0
    element = []

    for i in list(molstart):
        try:
            check = int(i)
            if atomcount != 0:
                while atomcount != 0:
                    atoms.append(ele)
                    atomcount -= 1
                atomcount = int(i)
                continue
            if atomcount == 0:
                atomcount = int(i)
                continue
        except:
            if i == "E":
                while atomcount != 0:
                    atoms.append(ele)
                    atomcount -= 1
                break
            if element != []:
                if str.islower(i) == True:
                    element.append(i)
                    ele = ''.join(element)
                    while atomcount != 0:
                        atoms.append(ele)
                        atomcount -= 1
                    if atomcount == 0:
                        ele= ""
                if str.isupper(i) == True:
                    element = [i]
                    ele = i
                    continue
            if element == []:
                element = [i]
                ele = i
                continue

def randomizer(atoms, molecule, moleinfo):
    shuffle(atoms)
    
    for i in atoms:
        if len(molecule) != 0:
            atom = i
            achoose = choice(molecule)
            bondlength = covalent_radii[achoose.number] + covalent_radii[atomic_numbers[atom]]
            theta = uniform(0, 2*pi)
            fi = uniform(0, pi)
            i = Atom(atom, [achoose.position[0] + (bondlength*sin(fi)*cos(theta)), achoose.position[1] +
                            (bondlength*sin(fi)*sin(theta)), achoose.position[2] + (bondlength*cos(fi))])
            molecule.append(i)
            continue
        if len(molecule) == 0:
            atom = i
            i = Atom(atom, position = (uniform(-10.,10.), uniform(-10.,10.), uniform(-10.,10.)))
            molecule.append(i)
            continue
    mole = molecule.get_potential_energy()
    moleinfo[mole] = molecule
    return moleinfo
    return molecule

#calc = EMT()
#molecule.set_calculator(calc)

#opt = MinimaHopping(atoms = molecule)
#opt(totalsteps = 10)


#write("molecule.traj", molecule)


#number = 1
#plusone(number)
#print number

#atoms = ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]
#randomizer(atoms)
#print molecule

