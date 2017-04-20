from pickle import load
from ase.io import write


f = open("moleculedata.dat", 'rb')
biginfo = load(f)

atcount = 0
h_energy = 100
lowest_energy = 100000
for i in biginfo:
    print i
    if i < lowest_energy:
        lowest_energy = i
    if i > h_energy:
        h_energy = i
    atcount += 1

print "number of molecules", atcount
print 'lowest energy from this trial found'
print lowest_energy
print biginfo[lowest_energy]
print 'highest energy'
print h_energy

            
write('lowestenergy.traj', biginfo[lowest_energy])
write('highestenergy.traj', biginfo[h_energy])
