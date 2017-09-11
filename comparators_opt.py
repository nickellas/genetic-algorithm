from copy import copy
from scipy.special import erf
import numpy as np
from itertools import product
from ase.ga import get_raw_score, set_raw_score

def cosine_dist(f1,f2):
    norm1 = np.linalg.norm(f1)
    norm2 = np.linalg.norm(f2)
    distance = np.sum(np.array(f1)*np.array(f2))/(norm1*norm2)

    cos_dist = 0.5*(1-distance)

    return cos_dist

def norm2_dist(f1,f2):
    distance = np.linalg.norm(np.array(f1)-np.array(f2))

    return distance

def norm1_dist(f1,f2):
    distance = np.sum(abs(np.array(f1)-np.array(f2)))

    return distance

class BagOfBonds(object):

    def __init__(self, n_top=None, pair_cor_cum_diff=0.015,
                 pair_cor_max=0.7, dE=0.5, mic=True, excluded_types=[],
                 overwrite=True):
        self.pair_cor_cum_diff = pair_cor_cum_diff
        self.pair_cor_max = pair_cor_max
        self.dE = dE
        self.n_top = n_top or 0
        self.mic = mic
        self.excluded_types = excluded_types
        self.overwrite = overwrite

    def get_features(self,a):
        """ Utility method used to calculate interatomic distances
            returned as a dict sorted after atomic type
            (e.g. (6,6), (6,42), (42,42)).
        """
        if not self.overwrite:
            if 'features' in a.info:
                return a.info['features']

        atoms = a[-self.n_top:]

        unique_types = sorted(list(set(atoms.numbers)))
        unique_types = [u for u in unique_types if u not in self.excluded_types]
        pair_cor = {}
        for idx, u1 in enumerate(unique_types):
            i_u1 = [i for i in range(len(atoms)) if atoms[i].number == u1]
            for u2 in unique_types[idx:]:
                i_u2 = [i for i in range(len(atoms)) if atoms[i].number == u2]
                d = []
                if u1 == u2:
                    for i, n1 in enumerate(i_u1):
                        for n2 in i_u2[i+1:]:
                            d.append(float(atoms.get_distance(n1,n2,self.mic)))
                else:
                    for i, n1 in enumerate(i_u1):
                        for n2 in i_u2:
                            d.append(float(atoms.get_distance(n1,n2,self.mic)))

                d.sort()
                if len(d) == 0:
                    continue
                pair_cor[(u1,u2)] = d

        a.info['features'] = pair_cor
        return pair_cor

    def get_similarity(self,f1,f2):
        """ Method for calculating the similarity between two objects with features
        f1 and f2, respectively.
        """
        # Calculate similarity.
        d1 = np.hstack(zip(*sorted(f1.items()))[1])
        d2 = np.hstack(zip(*sorted(f2.items()))[1])

        d_norm = np.sum(np.mean([d1,d2],axis=0))

        df = np.abs(np.array(d1)-np.array(d2))
        cum_diff = np.sum(df)

        s = cum_diff / d_norm

        return s

    def looks_like(self, a1, a2):
        # Energy criterium
        try:
            dE = abs(a1.get_potential_energy() - a2.get_potential_energy())
        except:
            dE = abs(get_raw_score(a1) - get_raw_score(a2))
        if dE >= self.dE:
            return False

        # Structure criterium
        f1 = self.get_features(a1)
        f2 = self.get_features(a2)

        d1 = sum(f1.values(),[])
        d2 = sum(f2.values(),[])

        max_d = max(np.abs(np.array(d1)-np.array(d2)))
        s = self.get_similarity(f1,f2)

#        print s, max_d
        if s > self.pair_cor_cum_diff or max_d > self.pair_cor_max:
            return False
        else:
            return True

class SeaOfBonds(object):

    def __init__(self, a_opt, n_top=None, rcut=20., binwidth=0.5,
                 excluded_types=[], sigma='auto', nsigma=4,
                 pbc=[True,True,False]):
        self.n_top = n_top or 0
        self.numbers = a_opt[-self.n_top:].numbers
        self.cell = a_opt.get_cell()
        self.rcut = rcut
        self.binwidth = binwidth
        self.excluded_types = excluded_types
        if sigma == 'auto':
            self.sigma = 0.5*binwidth
        else:
            self.sigma = sigma
        self.nsigma = nsigma
        self.pbc = pbc

    def get_features(self, a):
        atoms = a[-self.n_top:]

        unique_types = sorted(list(set(self.numbers)))
        unique_types = [u for u in unique_types if u not in self.excluded_types]

        # Include neighboring unit cells within rcut
        cell_vec_norms = np.apply_along_axis(np.linalg.norm,0,self.cell)
        cell_neighbors = []
        for i in range(3):
            ncell_max = int(np.floor(self.rcut/cell_vec_norms[i]))
            if self.pbc[i]:
                cell_neighbors.append(range(-ncell_max,ncell_max+1))
            else:
                cell_neighbors.append([0])

        # Binning parameters
        m = int(np.ceil(self.nsigma*self.sigma/self.binwidth))
        nbins = int(np.ceil(self.rcut*1./self.binwidth))
        smearing_norm = erf(self.binwidth*(2*m+1)*1./(np.sqrt(2)*self.sigma)) # Correction to nsigma cutoff

        # Get interatomic distances
        pos = atoms.get_positions()

        pair_cor = {}
        for idx, u1 in enumerate(unique_types):
            i_u1 = [i for i,atom in enumerate(atoms) if atom.number == u1]
            for u2 in unique_types[idx:]:
                i_u2 = [i for i,atom in enumerate(atoms) if atom.number == u2]
                rdf = np.zeros(nbins)
                for n1 in i_u1:
                    pi = np.array(pos[n1])
                    for disp_vec in product(*cell_neighbors):
                        displacement = np.dot(self.cell.T,np.array(disp_vec).T)
                        if u1 == u2 and disp_vec == (0,0,0): # Avoid self-counting in unit cell
                            pj = [p for n2,p in enumerate(pos) if n2 in i_u2 and n2 == n1]
                        else:
                            pj = [p for n2,p in enumerate(pos) if n2 in i_u2]
                        if len(pj) == 0:
                            continue

                        displaced_pos = pj + displacement
                        ds = np.apply_along_axis(np.linalg.norm,1,
                                                 displaced_pos-pi) # Interatomic distances

                        for dij in ds:
                            rbin = int(np.floor(dij/self.binwidth)) # Bin for dij
                            for i in range(-m,m+1): # Bins corresponding to (+/-)nsigma*sigma from rbin
                                newbin = rbin + i
                                if newbin < 0 or newbin >= nbins:
                                    continue

                                c = 1./(np.sqrt(2)*self.sigma)
                                value = 0.5*erf(c*((newbin+1)*self.binwidth - dij)) - 0.5*erf(c*(newbin*self.binwidth - dij))

                                value /= smearing_norm
                                rdf[newbin] += value

                pair_cor[(u1,u2)] = rdf

        return pair_cor

    def get_similarity(self,f1,f2):
        """ Method for calculating the similarity between two objects with features
        f1 and f2, respectively.
        """

        # Calculate similarity.
        s = 0.
        for key in sorted(f1.keys()):
            d1 = f1[key]
            d2 = f2[key]

            while len(d1) < len(d2):
                d1.append(0.)
            while len(d2) < len(d1):
                d2.append(0.)

            d_norm = np.sum(np.mean([d1,d2],axis=0))

            df = np.abs(np.array(d1)-np.array(d2))
            cum_diff = np.sum(df)


            s += cum_diff / d_norm

        return s

    def cosine_dist(self,f1,f2,numbers):
        keys = sorted(f1)
        unique_numbers = sorted(list(set(numbers)))
        typedic = {}
        for u in unique_numbers:
            typedic[(u)] = sum([u == i for i in numbers])

        w = {}
        wtot = 0
        for key in keys:
            while len(f1[key]) < len(f2[key]):
                f1[key].append(0.)
            while len(f2[key]) < len(f1[key]):
                f2[key].append(0.)

            weight = typedic[key[0]]*typedic[key[1]]
            wtot += weight
            w[key] = weight
        for key in keys:
            w[key] *= 1./wtot

        norm1 = 0
        norm2 = 0
        for key in keys:
            norm1 += (np.linalg.norm(f1[key])**2)*w[key]
            norm2 += (np.linalg.norm(f2[key])**2)*w[key]
        norm1 = np.sqrt(norm1)
        norm2 = np.sqrt(norm2)

        distance = 0.
        for key in keys:
            distance += np.sum(np.array(f1[key])*np.array(f2[key]))*w[key]/(norm1*norm2)

        distance = 0.5*(1-distance)

        return distance

#        norm1 = 0.
#        norm2 = 0.
#        for key in sorted(f1.keys()):
#            c1 = f1[key]
#            c2 = f2[key]

#            while len(c1) < len(c2):
#                c1.append(0.)
#            while len(c2) < len(c1):
#                c2.append(0.)

#            norm1 += np.linalg.norm(c1)
#            norm2 += np.linalg.norm(c2)

#        for key in sorted(f1.keys()):
#            distance += np.sum(np.array(c1)*np.array(c2))/(norm1*norm2)

#        cos_dist = 0.5*(1-distance)

#        return cos_dist

