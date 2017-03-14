# Definition of cracking patterns
# Aleksander Drenik
# Jozef Stefan Institute
# 2015

# Construction of cracking patterns of hydrogen containing molecules
# (water, methane, ammonia) with arbitrary H/(D+H) ratios.

# Cracking pattern intensities are read out from calibration file
# if reading of calib. file fails, default values are loaded (noted in "version")

# TO DO: CP calibration data quality check

# TO DO: CP weighted with absolute intensities from CP


# coding: utf-8


import scipy as sp
import math
#import sympy
from copy import copy

version = '1.01'

# definition of probability functions for the construction of hydrogen-containing molecules

def b_c(n,k):
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

def b_d(k,n,p):
    """binomial distribution probabiliy, k - number of events, n - number of tries, p - probability of event"""
    return b_c(n,k) * p**k * (1-p)**(n-k)
    
def branch_one(orig):
    "Isotope branching in removal of one atom. Input: dictionary: {atom: number}"
    branch = []
    for atom in orig.keys():
        if orig[atom] == 0:
            continue
        prob = float(orig[atom]) / sum(orig.values())
        br = copy(orig)
        br[atom] -= 1
        branch.append([prob, br])
    return branch

def branch_all(mol):
    "Isotope branching in removal of all atoms. Input: 0-level branch dictionary"
    branches = {0: [[1,mol]]}

    for i in range(1, sum(mol.values())+1):
        branches[i] = []
        for br_pair in branches[i - 1]:
            prob_main, branch = br_pair
            branches_local = branch_one(branch)
            for local_branch in branches_local:
                prob_loc, br_loc = local_branch
                branches[i].append([prob_loc * prob_main, br_loc])
    return branches
    
def branch_to_mass(branches):
    "Calculate masses of isotope branches. Input: branch dictionary. Global: md - atom mass dictionary"
    masses = {}
    for br_n in branches.keys():
        masses[br_n] = {}
        for br_pair in branches[br_n]:
            prob, branch = br_pair
            mass = 0
            for atom in branch.keys():
                mass += atom_mass_d[atom] * branch[atom]
            if mass in masses[br_n].keys():
                masses[br_n][mass] += prob
            else:
                masses[br_n][mass] = prob
    return masses
    

# dictionary of atom masses
atom_mass_d = {'D':2, 'H':1, 'T': 3}

# dictionary of known H-containing molecules:
# {name: name_in_calib, number of H isotopes, total mass of non-H atoms}
H_molecules_d = {}
H_molecules_d['ammonia'] = ['ammonia', 3, 14]
H_molecules_d['ammonia15'] = ['ammonia', 3, 15]
H_molecules_d['water'] = ['water', 2, 16]
H_molecules_d['water18'] = ['water', 2, 18]
H_molecules_d['methane'] = ['methane', 4, 12]
H_molecules_d['methane13'] = ['methane', 4, 13]
H_molecules_d['hydrogen'] = ['hydrogen', 2, 0]

# Default calibration

default_calib = {"version": "default", "device": "generic",
                 "water":[0.23,0.011,0.9],
                 "ammonia":[0.8,0.075,1.3],
                 "methane":[0.858,0.156,1.6],
                 "hydrogen":[0.05,0,1],
                 "N2":[0.072,0.008,1.0],
                 "O2":[0.114,0.008,0.96],
                 "Ar":[0.107,0.003,1.2],
                 "CO2":[0.114,0.085,1.4],
                 "CO":[0.045,0.009,1.05],
                 "Ne":[0.099,0.003,0.23]}

class mass_space:
    
    def __init__(self, max_mass=44):
        max_mass = max(max_mass, 44)
        base = sp.zeros((max_mass+1,max_mass+1))
        for i in range(1,max_mass+1):
            base[i][i]=1
        self.base = base
        self.max_mass = max_mass
        
    def init_CP(self, CP_spec=None):
        
        got_calib_file = False
        # If calibration is provided, the default calibration will be overridden.
        # If the calibration is incomplete, data from default calibration will remain.
        
        self.calib = default_calib
        if CP_spec != None:
            if type(CP_spec) == str:
                filename = CP_spec
                try:
                    calib_file = open(filename, "r").read().splitlines()
                    got_calib_file = True
                except IOError:
                    got_calib_file = False
                # print "calibration file not found"
                # using default values
            elif type(CP_spec) == dict:
                # TO Do - drugacen handling importa kot pri fajlu
                for key in CP_spec.keys():
                    self.calib[key] = CP_spec[key]
                #got_calib_file = True
        
        if got_calib_file:
            #self.calib = {}
            for line in calib_file:
                if len(line) == 0:
                    continue
                if line[0] == "#":
                    continue
                elements = line.split("\t")
                if elements[0] in ["version", "device"]:
                    self.calib[elements[0]] = elements[1]
                else:
                    self.calib[elements[0]] = map(float,elements[1:])    
        
        # Cracking patterns of hydrogen-containing molecules
        # Cracking patterns and relative intensities of all isotope configurations

        CP = {}
        base = self.base
        undefined = []

        # Ammonia
        # Cracking patterns from calibration file
        try:
            am1, am2, amr = self.calib["ammonia"]
            CP["NH3"] = amr*(base[17] + am1*base[16] + am2*base[15])
            CP["NDH2"] = amr*(base[18] + am1*(0.6667*base[17] + 0.3333*base[16])+am2*(0.3333*base[16]+0.6667*base[15]))
            CP["ND2H"] = amr*(base[19] + am1*(0.3333*base[18] + 0.6667*base[17]) + am2*(0.6667*base[16] + 0.3333*base[15]))
            CP["ND3"] = amr*(base[20] + am1*base[18] + am2*base[16])
        except:
            undefined.append("ammonia")

        # Water
        # Cracking patterns from calibration file
        try:
            wa1, wa2, war = self.calib["water"]
            CP["D2O"] = war*(base[20] + wa1*base[18] + wa2*base[16])
            CP["DHO"] = war*(base[19] + .5*wa1*base[18] + .5*wa1*base[17] + wa2*base[16])
            CP["H2O"] = war*(base[18] + wa1*base[17] + wa2*base[16])
        except:
            undefined.append("water")

        # Methane
        # Cracking patterns from calibration file

        try:
            me1, me2, mer = self.calib["methane"]
            CP["CD4"]= mer*(base[20] + me1*base[18] + me2*base[16])
            CP["CD3H"]= mer*(base[19] + me1*(0.25*base[18] + 0.75*base[17]) + me2*(0.5*base[16] + 0.5*base[15]))
            CP["CD2H2"]= mer*(base[18] + me1*(0.5*base[17] + 0.5*base[16]) + me2*(0.1667*base[16] + 0.6667*base[15] + 0.1667*base[14]))
            CP["CDH3"]= mer*(base[17] + me1*(0.75*base[16] + 0.25*base[15]) + me2*(0.5*base[15] + 0.5*base[14]))
            CP["CH4"] = mer*(base[16] + me1*base[15] + me2*base[14])
        except:
            undefined.append("methane")
            
        # Hydrogen
        # Cracking patterns from calibration file
        try:
            hy1, hy2, hyr = self.calib["hydrogen"]
        except:
            undefined.append("hydrogen")

        # Cracking patterns of non-H-containing molecules
        
        # Nitrogen
        # Cracking patterns from calibration file

        try:
            ni1, ni2, nir = self.calib["N2"]
            CP["N2"] = nir*(base[28] + ni1*base[14] + ni2*base[29])
        except:
            undefined.append("N2")
            
        # Nitrogen 15N2
        # Cracking patterns from 14N2 entry

        try:
            ni1, ni2, nir = self.calib["N2"]
            CP["15N2"] = nir*(base[30] + ni1*base[15] + ni2*base[29])
        except:
            undefined.append("N2")
        

        # Oxygen and Oxygen 18O2
        # Cracking patterns from calibration file

        try:
            ox1, ox2, oxr = self.calib["O2"]
            CP["O2"] = oxr*(base[32] + ox1*base[16] + ox2*base[34])
            CP["18O2"] = oxr*(base[36] + ox1*base[18] + ox2*base[34])
        except:
            undefined.append("O2")
           
    

        # Argon
        # Cracking patterns from calibration file

        try:
            ar1, ar2, arr = self.calib["Ar"]
            CP["Ar"] = arr*(base[40] + ar1*base[20] + ar2*base[36])
        except:
            undefined.append("Ar")

        # Carbon dioxide
        # Cracking patterns from calibration file

        try:
            cd1, cd2, cdr = self.calib["CO2"]
            CP["CO2"] = cdr*(base[44] + cd1*base[28] + cd2*base[16])
        except:
            undefined.append("CO2")

        # Carbon monoxide
        # Cracking patterns from calibration file

        try:
            cm1, cm2, cmr = self.calib["CO"]
            CP["CO"] = cmr*(base[28] + cm1*base[12] + cm2*base[16])
        except:
            undefined.append("CO")

        # Neon
        # Cracking patterns from calibration file

        try:
            ne1, ne2, ner = self.calib["Ne"]
            CP["Ne"] = ner*(base[20] + ne1*base[22] + ne2*base[21])
        except:
            undefined.append("Ne")

        self.undefined = undefined
        self.CP = CP

    # Hydrogen containing molecules, cracking pattern construction
    # Construction of the molecule isotopologue composition
    # The molecule is constructed from all isotopologues. The "concentration" of each isotopologue is calculated by the
    # binomial probability distribution:
    # Probability = b_d (number of H atoms, number of H + D atoms, H/(H+D) ratio)
    def construct_full(self, NH_mass, nAt, p, peak_list):
        "construct cracking patterns. input: non-H mass, total number of H isotopes, H/(H+D), CP_list"
        pattern = []
        for nH in range(nAt + 1):
            nH_pattern = []
            prob = b_d(nH, nAt, p)
            mol = {'H': nH, 'D': nAt - nH}
            branches = branch_all(mol)
            mass_branches = branch_to_mass(branches)
            num_of_peaks = min(len(mass_branches), len(peak_list))
            for i in range(num_of_peaks):
                temp = []
                peak = peak_list[i]
                m_br = mass_branches[i]
                for mass in m_br.keys():
                    m_prob = m_br[mass]
                    temp.append(m_prob * self.base[mass + NH_mass])
                #print peak * sum(temp)
                nH_pattern.append(peak * sum(temp))
            pattern.append(prob * sum(nH_pattern))
        return sum(pattern)
    
    def construct(self, molecule_name, ratio):
        calib_name, nAt, NH_mass = H_molecules_d[molecule_name]
        calib_list = self.calib[calib_name]
        intensity = calib_list[-1]
        peak_list = [1] + calib_list[:-1]
        return intensity * self.construct_full(NH_mass, nAt, ratio, peak_list)

known_hydrogen_species = ["ammonia","water","methane"]
known_other_species = ["Ar","O2","CO2","CO","Ne","N2"]


