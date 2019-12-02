# Definition of cracking patterns
# Aleksander Drenik
# Jozef Stefan Institute
# 2015
# IPP Garching
# 2018

# Construction of cracking patterns of hydrogen containing molecules
# (water, methane, ammonia) with arbitrary H/(D+H) ratios.

# Cracking pattern intensities are read out from calibration file
# if reading of calib. file fails, default values are loaded (noted in "version")

# coding: utf-8


import scipy as sp
import math

# import sympy
from copy import deepcopy
from .RGA_calibration_reader import read_calibration_file, read_old_version

version = "2.0b"

# definition of probability functions for the construction of hydrogen-containing molecules


def b_c(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))


def b_d(k, n, p):
    """binomial distribution probabiliy, k - number of events, n - number of tries, p - probability of event"""
    return b_c(n, k) * p ** k * (1 - p) ** (n - k)


def branch_one(orig):
    "Isotope branching in removal of one atom. Input: dictionary: {atom: number}"
    branch = []
    for atom in orig:
        if orig[atom] == 0:
            continue
        prob = float(orig[atom]) / sum(orig.values())
        br = deepcopy(orig)
        br[atom] -= 1
        branch.append([prob, br])
    return branch


def branch_all(mol):
    "Isotope branching in removal of all atoms. Input: 0-level branch dictionary"
    branches = {0: [[1, mol]]}

    for i in range(1, sum(mol.values()) + 1):
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
    for br_n in branches:
        masses[br_n] = {}
        for br_pair in branches[br_n]:
            prob, branch = br_pair
            mass = 0
            for atom in branch:
                mass += atom_mass_d[atom] * branch[atom]
            if mass in masses[br_n]:
                masses[br_n][mass] += prob
            else:
                masses[br_n][mass] = prob
    return masses


# dictionary of atom masses
atom_mass_d = {"D": 2, "H": 1, "T": 3}

# necessary keys for cracking pattern definitions
necessary_keys = {
    "H": {"non-H-mass", "H-atoms", "intensity", "peaks"},
    "non-H": {"intensity", "peaks"},
}


def check_calib(calib):
    # check the definitions of cracking patterns
    # the calib dictionary must include 'H' and 'non-H'
    for key in ["H", "non-H"]:
        if key not in calib:
            calib[key] = {}

    # H-molecules
    # the definition must include the following keys:
    # non-H-mass: total mass of the non-hydrogen atoms
    # H-atoms: number of hydrogen atoms
    # peaks: relative intensities of the peaks
    # intensity: relative intensity of the main peak

    # non-H-molecules
    # the definition must include the following keys:
    # peaks: masses and relative intensities of the peaks
    # intensity: relative intensity of the main peak

    # check definitions of cracking patterns
    for which in ["H", "non-H"]:
        for mol in calib[which]:
            defs = calib[which][mol]
            if len(necessary_keys[which] - set(defs.keys())) > 0:
                # remove the definition from the list
                del calib[which][mol]
    return calib


def check_peaks(CPline):
    """Returns a list of non-zero masses in the cracking pattern, and a list of [mass, intensity] pairs"""
    pks = []
    pk_s = []
    for i in range(len(CPline)):
        if CPline[i] != 0:
            pks.append([i, CPline[i]])
            pk_s.append(i)
    return pk_s, pks


class mass_space:
    def __init__(self, max_mass=44):
        self.make_base(max_mass)

    def make_base(self, max_mass):
        base = sp.zeros((max_mass + 1, max_mass + 1))
        for i in range(1, max_mass + 1):
            base[i][i] = 1
        self.base = base
        self.max_mass = max_mass

    def init_CP(self, CP_spec=None):

        self.CP = {}
        self.calib = {"device": "None", "version": "None", "H": {}, "non-H": {}}
        # read out the cracking pattern definitions from file.
        # ATM no default version
        if type(CP_spec) == str:
            err, calib = read_calibration_file(CP_spec)
            self.calib = check_calib(calib)
            self.calib["source"] = CP_spec
            self.calib["read error"] = err
        elif type(CP_spec) == dict:
            calib = deepcopy(CP_spec)
            self.calib = check_calib(calib)
            self.calib["source"] = "direct entry"
        # if the calibration contains no definitions treat the input as the old version
        if len(self.calib["H"]) == 0 and len(self.calib["non-H"]) == 0:
            calib = read_old_version(CP_spec)
            self.calib = check_calib(calib)

        # check if the base contains sufficiently high masses
        # currently works only for non-H molecules
        used_peaks = []
        for gas in calib["non-H"]:
            peaks = calib["non-H"][gas]["peaks"]
            for peak in peaks:
                used_peaks.append(peak[0])
        # check for H molecules, assuming the heaviest isotope is deuterium
        for gas in calib["H"]:
            CPdef = calib["H"][gas]
            max_mass_H = CPdef["non-H-mass"] + CPdef["H-atoms"] * atom_mass_d["D"]
            used_peaks.append(max_mass_H)
        try:
            max_peak = max(used_peaks)
        except ValueError:
            max_peak = 0
        if max_peak > self.max_mass:
            self.max_mass = max_peak
            self.make_base(max_peak)

        # definition of the non hydrogen molecules
        for gas in self.calib["non-H"]:
            vector = deepcopy(self.base[0])
            CPdef = self.calib["non-H"][gas]
            for peak in CPdef["peaks"]:
                vector += self.base[peak[0]] * peak[1]
            vector *= CPdef["intensity"]
            self.CP[gas] = vector

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
            mol = {"H": nH, "D": nAt - nH}
            branches = branch_all(mol)
            mass_branches = branch_to_mass(branches)
            num_of_peaks = min(len(mass_branches), len(peak_list))
            for i in range(num_of_peaks):
                temp = []
                peak = peak_list[i]
                m_br = mass_branches[i]
                for mass in m_br:
                    m_prob = m_br[mass]
                    temp.append(m_prob * self.base[mass + NH_mass])
                # print peak * sum(temp)
                nH_pattern.append(peak * sum(temp))
            pattern.append(prob * sum(nH_pattern))
        return sum(pattern)

    def construct(self, molecule_name, ratio):
        CPdef = self.calib["H"][molecule_name]
        nAt = CPdef["H-atoms"]
        NH_mass = CPdef["non-H-mass"]
        intensity = CPdef["intensity"]
        # New in this version: the peak list in calibration must inlcude the first peak, too
        # (previously 1 by default)
        peak_list = CPdef["peaks"]
        return intensity * self.construct_full(NH_mass, nAt, ratio, peak_list)
