# funkcije za fitanje RGA posnetkov
# Aug 2016
# Aleksander Drenik
# IJS, Ljubljana
# IPP, Garching

# Merged with new candidates (30.11.2018)
# TO DO: Export candidates, consider import candidates

import sympy
import scipy as sp
import scipy.optimize as so
import math
import time
from copy import deepcopy


version = '1.02'

def check_candidates(hydrogen_species,non_H_species):
    #check if format of candidate molecules is OK
    errors = []
    for key in hydrogen_species.keys():
        specimen = hydrogen_species[key]
        pressure = specimen[0]
        mol_name = specimen[1]
        ratio_raw = specimen[2]
        err_val = 0
        if mol_name not in molecules.known_hydrogen_species:
            print "Molecule %s not among hydrogen species" %mol_name
            err_val = 1
        if type(pressure) != str:
            print "Invalid pressure parameter type for %s." %mol_name
            err_val = 10
        if type(mol_name) != str:
            print "Invalid molecule name type for %s." %mol_name
            err_val = 2
        if type(ratio_raw) == list:
            if len(ratio_raw) != 3:
                print "Wrong H/(D+H) ratio definition for %s" %mol_name
                err_val = 50
            else:
                ratio, ratio_min, ratio_max = ratio_raw
                if type(ratio) != str:
                    print "Invalid H/(D+H) ratio parameter type for %s" %mol_name
                    err_val = 51
                if type(ratio_min) not in [int,float] or type(ratio_max) not in [int,float]:
                    print "Invalid H/(D+H) ratio boundary type for %s" %mol_name
                    err_val = 52
                if ratio_min > ratio_max or not(0 <= ratio_min <= 1) or not(0 <= ratio_max <= 1):
                    print "Invalid H/(D+H) ratio boundary definition for %s" %mol_name
                    err_val = 55
        elif type(ratio_raw) not in [str,int,float]:
            print "Invalid H/(D+H) ratio parameter type for %s" %mol_name
            err_val = 40
        errors.append(err_val)
        
    for key in non_H_species.keys():
        specimen = non_H_species[key]
        pressure = specimen[0]
        mol_name = specimen[1]
        err_val = 0
        if mol_name not in molecules.known_other_species:
            print "Molecule %s not among non-hydrogen species" %mol_name
            err_val = 1
        if type(pressure) != str:
            print "Invalid pressure parameter type for %s." %mol_name
            err_val = 10
        if type(mol_name) != str:
            print "Invalid molecule name type for %s." %mol_name
            err_val = 2
        errors.append(err_val)

    return errors


def make_candidates(molecules, hydrogen_species, non_H_species):

    candidates = []
    ratios = []
    temp_parameters = []
    pressures = []
    ratios = []
    boundaries = {'pressure': [], 'ratio': []}
    #init_vals = []
    for key in hydrogen_species.keys():
        specimen = hydrogen_species[key]
        pressure_raw = specimen[0]
        
        if type(pressure_raw) == list:
            if len(pressure_raw) != 3:
                # error, invalid pressure definition
                continue
            pressure = pressure_raw[0]
            pres_bnd = pressure_raw[1:]
        elif type(pressure_raw) == str:
            pressure = pressure_raw
            pres_bnd = [0,None]        
        mol_name = specimen[1]
        ratio_raw = specimen[2]
        if pressure not in temp_parameters:
            local_parameters = [pressure]
            temp_parameters.append(pressure)
            pressures.append(pressure)
            boundaries['pressure'].append(pres_bnd)
            #init_vals.append(1)
        #preveri, da je izotopsko razmerje v redu
        if type(ratio_raw) == list:
            if len(ratio_raw) != 3:
                print "Wrong D/(D+H) ratio definition for %s" %mol_name
                continue
            else:
                ratio = ratio_raw[0]
                ratio_bnd = ratio_raw[1:]
                if ratio not in temp_parameters:
                    local_parameters.append(ratio)
                    temp_parameters.append(ratio)
                    ratios.append(ratio)
                    boundaries['ratio'].append(ratio_bnd)
                    #init_vals.append(0.5)
        elif type(ratio_raw) == str:
            ratio = ratio_raw
            ratio_bnd = [0,1]
            if ratio not in temp_parameters:
                local_parameters.append(ratio)
                temp_parameters.append(ratio)
                ratios.append(ratio)
                boundaries['ratio'].append(ratio_bnd)
                #init_vals_te.append(0.5)
        elif type(ratio_raw) in [int,float]:
            ratio = ratio_raw

        par_string = ", ".join(local_parameters)

        strings_to_exec = []
        strings_to_exec.append("%s = sympy.symbols('%s')" %(par_string, par_string))
        strings_to_exec.append("local_candidate = %s * molecules.construct('%s', %s)" %(pressure, mol_name, ratio))
        for string in strings_to_exec:
            exec(string)

        candidates.append(local_candidate)
        
    for key in non_H_species.keys():
        specimen = non_H_species[key]
        pressure_raw = specimen[0]
        mol_name = specimen[1]
        if type(pressure_raw) == list:
            if len(pressure_raw) != 3:
                # error, invalid pressure definition
                continue
            pressure = pressure_raw[0]
            pres_bnd = pressure_raw[1:]
        elif type(pressure_raw) == str:
            pressure = pressure_raw
            pres_bnd = [0,None]
        if pressure not in temp_parameters:
            temp_parameters.append(pressure)
            pressures.append(pressure)
            boundaries['pressure'].append(pres_bnd)
        #init_vals_temp.append(1)
        exec("%s = sympy.symbols('%s')" %(pressure, pressure))
        exec("local_candidate = %s * molecules.CP['%s']" %(pressure, mol_name))
        
        candidates.append(local_candidate)
        
    # make a list of the fitting parameters and the initial values
    parameters = pressures + ratios
    init_vals = {'pressure': [], 'ratio': []}
    
    for what in ['pressure', 'ratio']:
        for pair in boundaries[what]:
            if pair[1] == None:
                init_val = 0
            else:
                init_val = 0.5 * (pair[1] - pair[0])
            init_vals[what].append(init_val)
    
    candidate_dict = {"candidates": candidates, "parameters": parameters, "pressures": pressures, "ratios": ratios,
    "boundaries": boundaries, "init_vals": init_vals}
        
    return candidate_dict
    
def make_calibration_candidates(molecules, mol_def, ratio_def, peak_defs):
    "Construct the candidates and parameters for fit_line"
    # Za fit potrebujem:
    # candidates
    # parameters
    # init_vals
    # boundaries
    # TO DO - popravi na novo verzijo kandidatov


    # najprej identifikacija molekule, ce je mogoce, in priprava definicije za construct_full
    # construct_full zahteva: NH_mass, nAt, p(AKA ratio), peak_list

    candidates = []
    parameters = []
    boundaries = []
    init_vals = []
    pressures = []
    ratios = []
    peaks = ['1']

    # priprava tlaka:

    parameters.append('pres')
    pressures.append('pres')
    boundaries.append([0,None])
    init_vals.append(1)

    # priprava izotopskega razmerja

    if type(ratio_def) == list:
        if len(ratio_def) != 3:
            print "Wrong D/(D+H) ratio definition"
            #continue
        else:
            ratio = ratio_def[0]
            ratio_bnd = ratio_def[1:]
            if ratio not in parameters:
                parameters.append(ratio)
                ratios.append(ratio)
                boundaries.append(ratio_bnd)
                init_vals.append(0.5 * (ratio_bnd[1] - ratio_bnd[0]))
    elif type(ratio_def) == str:
        ratio = ratio_def
        ratio_bnd = [0,1]
        if ratio not in parameters:
            parameters.append(ratio)
            ratios.append(ratio)
            boundaries.append(ratio_bnd)
            init_vals.append(0.5)
    elif type(ratio_def) in [int,float]:
        ratio = ratio_def


    for peak_def in peak_defs:
        if type(peak_def) == list:
            if len(peak_def) != 3:
                print "Invalid peak definition"
                continue
            if type(peak_def[0]) == str and type(peak_def[1]) in [float, int] and type(peak_def[2]) in [float, int]:
                parameters.append(peak_def[0])
                peaks.append(peak_def[0])
                boundaries.append([peak_def[1], peak_def[2]])
                init_vals.append(0.5 * (peak_def[2] - peak_def[1]))

        if type(peak_def) == str:
            parameters.append(peak_def)
            boundaries.append([0, None])
            init_vals.append(0.5)
            peaks.append(peak_def)
        if type(peak_def) in [float, int]:
            peaks.append(str(peak_def))

    par_string = ", ".join(parameters)
    peak_string = "[%s]" %", ".join(peaks)

    strings_to_exec = []
    strings_to_exec.append("%s = sympy.symbols('%s')" %(par_string, par_string))
    strings_to_exec.append("candidate = %s * %s * molecules.construct_full(%s, %s, %s, %s)" %('pres', mol_def['rel_int'], mol_def['NH_mass'], 
                                                            mol_def['nAt'], ratio, peak_string))
    for string in strings_to_exec:
            exec(string)

    candidates = [candidate]

    cand_dict = {'candidates': candidates, 'parameters': parameters, 'pressures': ['pres'], 'ratios': [ratio], 'boundaries': boundaries,
        'init_vals': init_vals, 'peaks': peaks}
    return cand_dict

def check_disregard(disregard):
    errors = []
    for element in disregard:
        err_val = 0
        if type(element) not in [int,list]:
            err_val = 1
            #wrong definition type
        if type(element) == list:
            if len(element) != 3:
                #Wrong length of list
                err_val = 21
            if type(element[0]) != int:
                #Wrong type of disregarded mass definition
                err_val = 22
            if type(element[1]) != int:
                #Wrong type of trigger mass definition
                err_val = 23
            if type(element[2]) not in [float, int]:
                #Threshold value not a number - rethink the int part
                err_val = 24
        errors.append(err_val)
    return errors

# definicija fit_line funkcije
# input: line, recorded, header_int, Hspecies, nonHspecies, disregard
# output: results(dictionary), calculated_masses

def fit_line(line, recorded_in, header_int, candidates_dict, disregard, n_iter=0):

    candidates = candidates_dict["candidates"]
    par_all = candidates_dict["parameters"]
    par_pres = candidates_dict["pressures"]
    par_rat = candidates_dict["ratios"]
    boundaries = candidates_dict["boundaries"]
    init_vals = candidates_dict["init_vals"]
    masses_of_interest = candidates_dict["masses_of_interest"]
    
    # testno - masses of interest se prebere iz candidates_dict    
    #masses_of_interest = []
    
    # tukaj sem ven vrgel header int - ce dela, ga lahko dam ven iz argumentov.
    #for mass in range(1,len(sum(candidates))):
    #    if sum(candidates)[mass] != 0:
    #        masses_of_interest.append(mass) 
            
    recorded = deepcopy(recorded_in)
    for element in disregard:
        if type(element) == int:
            if len(line) >= element:
                recorded[element] = 0
                #print "disregarding %s allways" %element
        if type(element) == list:
            if len(element) == 3:
                try:
                    if line[element[1]] > element[2]:
                        recorded[element[0]] = 0
                        #print "disregarding %s - on condition" %element[0]
                except:
                    #print "not disregarding shit"
                    pass
            
    try:
        maxvalexp = math.ceil(-math.log10(max(line[1:])))
    except ValueError:
        try:
            maxvalexp = math.ceil(-math.log10(-min(line[1:])))
        except ValueError:
            maxvalexp = 0
    maxval = 10**maxvalexp
    # apply the maxval to the pressure boundaries
    boundaries_final = []
    init_vals_final = []
    for i in range(len(boundaries['pressure'])):
        pair = deepcopy(boundaries['pressure'][i])
        init_val = deepcopy(init_vals['pressure'][i])
        if pair[1] != None:
            for j in (0, 1):
                pair[j] *= maxval
            init_val *= maxval
        boundaries_final.append(pair)
        init_vals_final.append(init_val)
    boundaries_final = boundaries_final + boundaries['ratio']
    #boundaries_final = boundaries['pressure'] + boundaries['ratio']
    init_vals_final = init_vals_final + init_vals['ratio']
    #init_vals_final = init_vals['pressure'] + init_vals['ratio']
    
    ansatz = sum(candidates)
    equations = []
    # ali gre i po masses_of_interest ali header_int?
    # potrebno je sesteti po vseh posnetih masah
    # masses_of_interest je tako ali tako definiran po tem, da je vsota kandidatov (i.e. ansatz) == 0
    for i in header_int:
        recording = line[i] * maxval
        equations.append(recorded[i]*(recording-ansatz[i])*(recording-ansatz[i]))
    residual_string = str(sum(equations))
    par_string = "%s = x" %", ".join(par_all)
    
    def residual(x):
        exec(par_string) in locals()
        return eval(residual_string)
    
    def calculate_masses(x):
        exec(par_string) in locals()
        output = []
        for mass in map(str,ansatz):
            output.append(eval(mass) * maxval**-1)
        return sp.array(output)
    
    start_time = time.clock()
    # if more than 0 iterations are specified, use the basinhopping function
    if n_iter == 0:
        rez = so.minimize(residual,init_vals_final,bounds=boundaries_final,method='L-BFGS-B')
    else:
        rez = so.basinhopping(residual, init_vals_final, niter = n_iter, minimizer_kwargs = dict(method='L-BFGS-B', bounds=boundaries_final))
    residual_value = math.sqrt(rez.fun*maxval**-2)
    loc_duration = time.clock() - start_time

    results={}
    
    for i in range(0,len(par_all)):
        results[par_all[i]]=rez.x[i]
    for prs in par_pres:
        results[prs] = results[prs] * maxval**-1
    results['duration'] = loc_duration
    results['residual_value'] = residual_value
    
    sim_masses = calculate_masses(rez.x)
    
    return results, sim_masses
    
def export_CP_old(CPdict):
    outtab = [["#Cracking patterns used in fit"],[]]
    outheader = ["#Gas","2nd peak","3rd peak","relative intensity"]
    for key in ["device", "version"]:
        try:
            line = [key, CPdict[key]]
        except:
            line = [key, "N/A"]
        outtab.append(line)
    outtab.append([])
    outtab.append(outheader)
    for key in ["water", "ammonia", "methane"]:
        try:
            line = [key] + CPdict[key]
        except:
            line = [key] + ["N/A","N/A","N/A"]
        outtab.append(line)
    for key in CPdict.keys():
        if key not in ["device", "version", "water", "ammonia", "methane"]:
            try:
                line = [key] + CPdict[key]
            except:
                line = [key] + ["N/A","N/A","N/A"]
            outtab.append(line)
    return outtab
    
def export_CP(calib):
    """Prepare the cracking pattern definition for export"""
    # Construction of the output calibration file
    out_new = [['# Definition of cracking patterns for rga_data_handling objects'], 
                ['# File format: Nov 2018']]
    #out_new.append(['# Generated by calibration_convert.py'])
    try:
        out_new.append(['# Converted from %s' %calib['source']])
    except KeyError:
        pass
    
    for what in ['device', 'version', 'source']:
        try:
            val = calib[what]
        except KeyError:
            val = 'unspecified'
            
        #out_new.append('%s\t%s' %(what, val))
        out_new.append([what, val])
        
    for what in ['comment', 'remark']:
        try:
            val = calib['# %s' %what]
        except KeyError:
            continue
            
        out_new.append([val])

    out_new.append([''])
    out_new.append(['# Definitions of cracking patterns start here.'])
    out_new.append([''])

    for gas in calib['H'].keys():
        out_new.append(['H-molecule', gas])
        CP_def = calib['H'][gas]
        for what in ['non-H-mass', 'H-atoms']:
            out_new.append([what, CP_def[what]])
        for peak in CP_def['peaks']:
            out_new.append(['peak', peak])
        for what in ['intensity']:
            out_new.append([what, CP_def[what]])
        out_new.append([''])
        
    for gas in calib['non-H'].keys():
        out_new.append(['non-H-molecule', gas])
        CP_def = calib['non-H'][gas]
        for peak in CP_def['peaks']:
            out_new.append(['peak', peak[0], peak[1]])
        for what in ['intensity']:
            out_new.append([what, CP_def[what]])
        out_new.append([''])
            
    #write_err = write_lines(filename, out_new)
    #return write_err
    return out_new
    
   
def export_candidates(H_species, non_H_species, disregard):

    outtab = []
    outtab.append(["# Hydrogen species"])
    for key in H_species:
        line = ["Hspecies:", key] + H_species[key]
        if type(line[-1]) == list:
            line = line[:-1] + line[-1]
        outtab.append(line)
    outtab.append(["# non-Hydrogen species"])
    for key in non_H_species:
        line = ["nonHspecies:", key] + non_H_species[key]
        outtab.append(line)
    outtab.append(["# Disregard"])
    for element in disregard:
        if type(element) != list:
            element = [element]
        line = ["Disregard:"] + element
        outtab.append(line)
    return outtab
