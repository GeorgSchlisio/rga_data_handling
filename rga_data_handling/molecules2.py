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
import sympy

version = '1.0'

# definition of probability functions for the construction of hydrogen-containing molecules

def b_c(n,k):
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

def b_d(k,n,p):
    """binomial distribution probabiliy, k - number of events, n - number of tries, p - probability of event"""
    return b_c(n,k) * p**k * (1-p)**(n-k)


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
        
        self.calib = {"version": "default", "device": "generic",
                 "water":[0.23,0.011,0.9],
                 "ammonia":[0.8,0.075,1.3],
                 "methane":[0.858,0.156,1.6],
                 "N2":[0.072,0.008,1.0],
                 "O2":[0.114,0.008,0.96],
                 "Ar":[0.107,0.003,1.2],
                 "CO2":[0.114,0.085,1.4],
                 "CO":[0.045,0.009,1.05],
                 "Ne":[0.099,0.003,0.23]}
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
                self.calib = CP_spec
                #got_calib_file = True
        
        if got_calib_file:
            self.calib = {}
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
        

        # Oxygen
        # Cracking patterns from calibration file

        try:
            ox1, ox2, oxr = self.calib["O2"]
            CP["O2"] = oxr*(base[32] + ox1*base[16] + ox2*base[34])
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
    def construct(self, molecule_name, ratio):
        CP = self.CP
        molecule={}
        molecule["water"] = b_d(2,2,ratio)*CP["H2O"] + b_d(1,2,ratio)*CP["DHO"] + b_d(0,2,ratio)*CP["D2O"]
        molecule["ammonia"] = b_d(3,3,ratio)*CP["NH3"] + b_d(2,3,ratio)*CP["NDH2"] + b_d(1,3,ratio)*CP["ND2H"] + b_d(0,3,ratio)*CP["ND3"]
        molecule["methane"] = b_d(4,4,ratio)*CP["CH4"] + b_d(3,4,ratio)*CP["CDH3"] + b_d(2,4,ratio)*CP["CD2H2"] + b_d(1,4,ratio)*CP["CD3H"] + b_d(0,4,ratio)*CP["CD4"]
        molecule["ammonia15"] = sp.delete(sp.insert(molecule["ammonia"],0,0),-1)

        try:
            outmol = molecule[molecule_name]
        except KeyError:
            outmol = []
        return outmol

known_hydrogen_species = ["ammonia","water","methane"]
known_other_species = ["Ar","O2","CO2","CO","Ne","N2"]


