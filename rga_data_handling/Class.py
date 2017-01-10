# RGA data container class

# Aug 2016
# Aleksander Drenik
# IJS, Ljubljana
# IPP, Garching


version = '2.011'
versions = {'Class': version}

import molecules2
import scipy as sp
import random
try:
    from RGA_fitting import *
    versions['RGA_fitting'] = version
    fitting_loaded = True
except ImportError:
    versions['RGA_fitting'] = version
    fitting_loaded = False
    
from os import path, mkdir

versions['molecules'] = molecules2.version

known_time_labels = ['time', 'hours']

# Import the default values

"""try:
    from __main__ import default_values
    default_loaded = True
except:
    default_loaded = False
    default_values = {}"""

def load_default_values(default_values):
    default = {}
    try:
        default['trace_name'] = default_values['trace_name']
    except KeyError:
        default['trace_name'] = 'Timetrace recording'
    try:
        default['profile_name']  = default_values['profile_name']
    except KeyError:
        default['profile_name'] = "Profile recording"
    try:
        default['export_dir'] = default_values['export_dir']
    except KeyError:
        default['export_dir'] = 'EXPORT'

    try:
        default['pulse_start'] = default_values['pulse_start']
    except KeyError:
        default['pulse_start'] = 0
        
    try:
        default['pulse_stop'] = default_values['pulse_stop']
    except KeyError:
        default['pulse_stop'] = None
        
    return default

def pin_point(input_list, sought_value):
    input_col = sp.array(input_list)
    temp_list = list(map(abs, input_col - sought_value))
    nearest_index = temp_list.index(min(temp_list))
    return nearest_index
    
def write_to_TSV(path, inputlist):
    outfile = open(path,'w')
    for line in inputlist:
        outfile.writelines("%s\r\n" %("\t".join(map(str,line))))
    outfile.close()

def PadRight(inputlist, TargetLength, FillValue):
    if len(inputlist) < TargetLength:
        return inputlist + [FillValue] * (TargetLength - len(inputlist))
    else:
        return inputlist

def perturb_CP(container, perturbation):
    CP_new = container.molecules.calib
    for key in CP_new.keys():
        entry = CP_new[key]
        if type(entry) == list:
            for i in [0, 1]:
                entry[i] = entry[i] * (1 + 2 * perturbation * (random.random() - 0.5))
    CP_new['version'] = "%s perturbed by %s" %(CP_new['version'], perturbation)
    container.replace_CP(CP_new)


#load_default_values({})


class Trace:
    # data handler - napisi vec
    # Internal errors: tag['init_errors']
    # 101: empty or invalid data dictionary
    # 102: empty data columns
    # 103: empty header_int
    
    def __init__(self, columns, tag, default=None):
        
        self.type = "Trace"
        self.columns = {}
        self.filled = False
        self.deconvoluted = False
        self.calibrated = False
        bloc_len = 0
        
        if default == None:
            self.def_val = load_default_values({})
        else:
            self.def_val = load_default_values(default)
        
        if type(tag) != dict:
            self.tag = {'raw_tag':tag}
        else:
            self.tag = tag
        self.tag['init_errors'] = []
        # najprej iz columnsov izlusci stolpce z maso in ostale stolpce
        header_int = []
        lengths = {}
        for key in columns.keys():
            # preverimo tudi, ce so vsi stolpci enake dolzine
            lengths[key] = len(columns[key])
            try:
                value = int(key)
                header_int.append(value)
                self.columns[value] = sp.array(columns[key])
            except ValueError:
                self.columns[key] = sp.array(columns[key])
        if len(header_int) == 0:
            self.tag['init_errors'].append(103)
        try:        
            bloc_len = min(lengths.values())
            for key in self.columns.keys():
                self.columns[key] = self.columns[key][:bloc_len]
            if bloc_len == 0:
                self.tag['init_errors'].append(102)
        except ValueError:
            self.tag['init_errors'].append(101)
        # naredi stolpec z indexi vrstic
        self.columns['index'] = sp.arange(bloc_len)
        if len(header_int) > 0 and bloc_len > 0:
            self.filled = True

            
            time_col_candidates = list(set(self.columns.keys()) & set(known_time_labels))
            try:
                self.set_timecol(time_col_candidates[0])
            except:            
                self.set_timecol('index')

            # Konstrukcija recorded
            self.header_int = sp.array(sorted(header_int))
            self.recorded = sp.zeros(self.header_int[-1] + 1)
            for i in range(self.header_int[-1] + 1):
                if i in self.header_int:
                    self.recorded[i] = 1

            # ime objekta - ce ni v tag-u, se avtoimenuje
            if 'title' not in self.tag.keys():
                self.tag['title'] = self.def_val['trace_name']
                
            # Initialize mass-space and cracking patterns
            self.molecules = molecules2.mass_space(self.header_int[-1])
            #TO DO - by default, initialize cracking patterns with appropriate calibration files
            self.molecules.init_CP()
            
            
    def replace_CP(self, path):
        if self.filled:
            self.molecules.init_CP(path)    
            
    def set_timecol(self, time_key):
        if time_key in self.columns.keys():
            self.time_col = time_key
            self.time_col_name = self.time_col.capitalize()
            if 'time_unit' in self.tag:
                if len(self.tag['time_unit']) > 0:
                    self.time_col_name = "%s [%s]" %(self.time_col_name, self.tag['time_unit'])
            
    def make_line(self, ti):
    # Produce a scipy array of the recorded signals, where index of the element in the array corresponds to the recorded mass
    # e.g. line[28] returns the recording at 28 AMU
    
        line = [self.columns[self.time_col][ti]]
        for mass in range(1,self.header_int[-1]+1):
            if mass in self.header_int:
                line.append(self.columns[mass][[ti]])
            else:
                line.append(0)
        line = sp.array(line)
        return line
    
    def make_MID_col(self, ti):
        # reproduce a line from the original recording - intensities at index ti
        # the corresponding mass column is self.header_int
        line = self.make_line(ti)
        MID_col = [line[mass] for mass in self.header_int]
        return sp.array(MID_col)
    
    def make_results_line(self, ti):
    # Produce a profile-like dictionary with results from the deconvolution.
    # Fro the ti-th line in the results time-trace
        resline = {}
        for gas in self.H_species.keys():
            resline[gas] = {'pressure': self.rescols[gas]['pressure'][ti], 'ratio': self.rescols[gas]['ratio'][ti]}
        for gas in self.non_H_species.keys():
            resline[gas] = {'pressure': self.rescols[gas]['pressure'][ti]}
        resline['residual'] = self.rescols['residual'][ti]
        resline['duration'] = self.rescols['duration'][ti]
        
        sim_MID_col = []
        for mass in self.masses_of_interest:
            sim_MID_col.append(self.simtracecol[mass][ti])
        sim_MID_col = sp.array(sim_MID_col)
        
        return resline, sim_MID_col

    
    def Make_Profile(self, wanted):
        ti = pin_point(self.columns[self.time_col], wanted)
        #line = self.make_line(ti)
        actual = self.columns[self.time_col][ti]
        MID_col = self.make_MID_col(ti)
        datatag = copy.copy(self.tag)
        datatag["mode"] = "export from timetrace"
        datatag["wanted"] = "%s: %s" %(self.time_col, wanted)
        datatag["actual"] = "%s: %s" %(self.time_col, actual)
        datatag["delay"] = "%s" %(actual - wanted)
        datatag["title"] = "%s @ %s" %(self.tag['title'], actual)
        if 'time_unit' in self.tag:
            for key in ['wanted','actual', 'delay', 'title']:
                datatag[key] = datatag[key] + ' %s' %self.tag['time_unit']
        datatag["data point"] = ti
        
        profile = Profile(self.header_int, MID_col, datatag, self.def_val)
        if self.deconvoluted:
            tri = pin_point(self.rescols[self.time_col], wanted)
            profile.deconvoluted = True
            resline, sim_MID_col = self.make_results_line(tri)
            profile.resline = resline
            profile.sim_MID_col = sim_MID_col
            profile.H_species = self.H_species
            profile.non_H_species = self.non_H_species
            profile.disregard = self.disregard
            profile.masses_of_interest = self.masses_of_interest
        return profile
    
    def make_candidates(self):
        
        self.candidates_dict = make_candidates(self.molecules, self.H_species, self.non_H_species)
        
    def deconvolute(self, H_species, non_H_species, disregard, start_time, stop_time, step):
        self.H_species = H_species
        self.non_H_species = non_H_species
        self.disregard = disregard
        
        self.make_candidates()
        self.parameters = self.candidates_dict['parameters']
        self.pressures = self.candidates_dict['pressures']
        self.ratios = self.candidates_dict['ratios']
        self.masses_of_interest = []
        for mass in range(1,self.header_int[-1] + 1):
            if sum(self.candidates_dict['candidates'])[mass] != 0:
                self.masses_of_interest.append(mass)
        
        ti_list = self.columns['index'][(start_time <= self.columns[self.time_col]) * (stop_time >= self.columns[self.time_col])] 
        ti_list = range(ti_list[0],ti_list[-1]+1,step)  
        
                
        simtrace = []
        outcols = {}
        outcols[self.time_col] = []
        outcols['duration'] = []
        outcols['residual'] = []
        for parameter in self.parameters:
            outcols[parameter] = []
            
        glob_start = time.clock()
        for ti in ti_list:
            line = self.make_line(ti)
            
            results, sim_masses = fit_line(line, self.recorded, self.header_int, self.candidates_dict, disregard)
            
            outcols[self.time_col].append(line[0])
            for prs in self.pressures:
                outcols[prs].append(results[prs])
            outcols['residual'].append(results['residual_value'])
            for ratio in self.ratios:
                outcols[ratio].append(results[ratio])
            outcols['duration'].append(results['duration'])

            simtrace.append(sim_masses)
            
        self.glob_duration = time.clock() - glob_start
        
        self.rescols = {}
        self.rescols['index'] = sp.array(ti_list)
        self.rescols[self.time_col] = sp.array(outcols[self.time_col])
        self.rescols['residual'] = sp.array(outcols['residual'])
        self.rescols['duration'] = sp.array(outcols['duration'])
        tablength = len(outcols[self.time_col])
        for key in self.H_species.keys():
            self.rescols[key]={}
            self.rescols[key]['pressure'] = sp.array(outcols[self.H_species[key][0]])
            if type(self.H_species[key][2]) == list:
                self.rescols[key]['ratio'] = sp.array(outcols[self.H_species[key][2][0]])
            elif type(self.H_species[key][2]) == str:
                self.rescols[key]['ratio'] = sp.array(outcols[self.H_species[key][2]])
            elif type(self.H_species[key][2] in [int, float]):
                self.rescols[key]['ratio'] = self.H_species[key][2] * sp.ones(tablength)
        for key in self.non_H_species.keys():
            self.rescols[key]={}
            self.rescols[key]['pressure'] = sp.array(outcols[self.non_H_species[key][0]])
            
        self.simtracecol={}
        for mass in range(1,len(sp.transpose(simtrace))):
            self.simtracecol[mass] = sp.transpose(simtrace)[mass]

        self.outtab = []
        header_line_1 = ['Quantity']
        header_line_2 = [self.time_col_name]
        self.outtab.append(self.rescols[self.time_col])
        for key in self.H_species.keys():
            header_line_2.append(key)
            header_line_1.append('Pressure')
            self.outtab.append(self.rescols[key]['pressure'])
        for key in self.non_H_species.keys():
            header_line_2.append(key)
            header_line_1.append('Pressure')
            self.outtab.append(self.rescols[key]['pressure'])
        for key in self.H_species.keys():
            header_line_2.append(key)
            header_line_1.append('H/(H+D)')
            self.outtab.append(self.rescols[key]['ratio'])
        for key in ['residual','duration']:
            header_line_1.append('Fitting')
            header_line_2.append(key)
            self.outtab.append(outcols[key])

        self.outtab = sp.array(self.outtab)
        self.outtab = sp.transpose(self.outtab)
        self.outtab = list(self.outtab)
        self.outtab = [header_line_2] + [header_line_1] + self.outtab

        massheader = [self.time_col_name]
        self.outsimtracecol = [sp.array(self.rescols[self.time_col])]
        for mass in self.masses_of_interest:
            massheader.append("%s AMU" %mass)
            self.outsimtracecol.append(self.simtracecol[mass])
        self.outsimtrace = list(sp.transpose(sp.array(self.outsimtracecol)))
        self.outsimtrace = [massheader] + self.outsimtrace
        self.deconvoluted = True
        
    def calibrate(self, molecule, ratio_def, peak_defs, disregard, start_time, stop_time, step):
        
        self.calib_tag = {}
        if type(molecule) == str:
            self.calib_tag['title'] = "Calibration of %s for %s" %(self.tag['title'], molecule)
            try:
                mol_def_1 = molecules2.H_molecules_d[molecule]
                rel_int = self.molecules.calib[mol_def_1[0]][-1]
                mol_def = {'NH_mass': mol_def_1[2], 'nAt': mol_def_1[1], 'rel_int': rel_int}
            except KeyError:
                #print "Unknown molecule"
                pass
        if type(molecule) == dict:
            self.calib_tag['title'] = "Calibration of %s" %self.tag['title']
            if not set(['NH_mass', 'nAt', 'rel_int']).issubset(set(mol_def.keys())):
                # print invalid molecule definition
                pass
            
        self.calib_candidates_dict = make_calibration_candidates(self.molecules, mol_def, ratio_def, peak_defs)        
        self.calib_masses_of_interest = []
        for mass in range(1,self.header_int[-1] + 1):
            if sum(self.calib_candidates_dict['candidates'])[mass] != 0:
                self.calib_masses_of_interest.append(mass)
        
        ti_list = self.columns['index'][(start_time <= self.columns[self.time_col]) * (stop_time >= self.columns[self.time_col])] 
        ti_list = range(ti_list[0],ti_list[-1]+1,step)  
        
                
        simtrace = []
        outcols = {}
        outcols[self.time_col] = []
        outcols['duration'] = []
        outcols['residual'] = []
        for parameter in self.calib_candidates_dict['parameters']:
            outcols[parameter] = []
            
        glob_start = time.clock()
        for ti in ti_list:
            line = self.make_line(ti)
            
            results, sim_masses = fit_line(line, self.recorded, self.header_int, self.calib_candidates_dict, disregard)
            
            outcols[self.time_col].append(line[0])
            
            for parameter in self.calib_candidates_dict['parameters']:
                outcols[parameter].append(results[parameter])
            
            outcols['duration'].append(results['duration'])
            outcols['residual'].append(results['residual_value'])

            simtrace.append(sim_masses)
            
        self.calib_glob_duration = time.clock() - glob_start
        
        # do sem sem prisel
        self.calibcols = {}
        self.calibcols['index'] = sp.array(ti_list)
        self.calibcols[self.time_col] = sp.array(outcols[self.time_col])
        self.calibcols['residual'] = sp.array(outcols['residual'])
        self.calibcols['duration'] = sp.array(outcols['duration'])
        self.calibcols['pressure'] = sp.array(outcols['pres'])
        tablength = len(outcols[self.time_col])
        if type(ratio_def) in [float, int]:
            self.calibcols['ratio'] = ratio_def * sp.ones(tablength)
        else:
            self.calibcols['ratio'] = sp.array(outcols[self.calib_candidates_dict['ratios'][0]])
        self.calibcols["peaks"] = {}
        for i in range(len(peak_defs)):
            peak_def = peak_defs[i]
            level = i + 1
            if type(peak_def) in [float, int]:
                self.calibcols['peaks'][level] = peak_def * sp.ones(tablength)
            else:
                if type(peak_def) == list:
                    peak_name = peak_def[0]
                else:
                    peak_name = peak_def
                self.calibcols['peaks'][level] = sp.array(outcols[peak_name])
  
        self.simtracecol={}
        for mass in range(1,len(sp.transpose(simtrace))):
            self.simtracecol[mass] = sp.transpose(simtrace)[mass]

        """self.outtab = []
        header_line_1 = ['Quantity']
        header_line_2 = [self.time_col_name]
        self.outtab.append(self.rescols[self.time_col])
        for key in self.H_species.keys():
            header_line_2.append(key)
            header_line_1.append('Pressure')
            self.outtab.append(self.rescols[key]['pressure'])
        for key in self.non_H_species.keys():
            header_line_2.append(key)
            header_line_1.append('Pressure')
            self.outtab.append(self.rescols[key]['pressure'])
        for key in self.H_species.keys():
            header_line_2.append(key)
            header_line_1.append('H/(H+D)')
            self.outtab.append(self.rescols[key]['ratio'])
        for key in ['residual','duration']:
            header_line_1.append('Fitting')
            header_line_2.append(key)
            self.outtab.append(outcols[key])

        self.outtab = sp.array(self.outtab)
        self.outtab = sp.transpose(self.outtab)
        self.outtab = list(self.outtab)
        self.outtab = [header_line_2] + [header_line_1] + self.outtab

        massheader = [self.time_col_name]
        self.outsimtracecol = [sp.array(self.rescols[self.time_col])]
        for mass in self.masses_of_interest:
            massheader.append("%s AMU" %mass)
            self.outsimtracecol.append(self.simtracecol[mass])
        self.outsimtrace = list(sp.transpose(sp.array(self.outsimtracecol)))
        self.outsimtrace = [massheader] + self.outsimtrace"""
        
        self.calibrated = True
        
    def subtract_background(self, time_window):
        """Subtraction of the background, defined as the average value during the specified time window"""
        
        # This function will override the original columns.     
        if self.filled:
  
            bgr_range = (self.columns[self.time_col] > time_window[0]) * (self.columns[self.time_col] < time_window[1])
            if list(bgr_range).count(True) > 1:
                do_bgr = True
            else:
                do_bgr = False
            if do_bgr:
                for mass in self.header_int:
                    self.columns[mass] = self.columns[mass] - sp.mean(self.columns[mass][bgr_range])

    
    def export(self, name=None, write_path=None, rec=False):
        
        if name == None:
            name = self.tag['title']
        if write_path == None:
            write_path = self.def_val['export_dir']
        if not path.exists(write_path):
            mkdir(write_path)
     
        resultname = path.join(write_path, "%s_TRACE_RESULTS.txt" %name)
        simtracename = path.join(write_path, "%s_TRACE_MASSES.txt" %name)
        cpname = path.join(write_path, "%s_TRACE_CALIBRATION.txt" %name)
        candname = path.join(write_path, "%s_TRACE_PARAMETERS.txt" %name)
        recname = path.join(write_path, "%s_TRACE_RECORDING.txt" %name)
        
        outrecording = [self.columns[self.time_col]]
        outrecheader = [self.time_col_name]
        for mass in self.header_int:
            outrecording.append(self.columns[mass])
            outrecheader.append("%s AMU" %mass)
        outrecording = [outrecheader] + list(sp.transpose(sp.array(outrecording)))       
        
        if self.deconvoluted:
            write_to_TSV(resultname, self.outtab)
            write_to_TSV(simtracename, self.outsimtrace)
            outCPtab = export_CP(self.molecules.calib)        
            write_to_TSV(cpname, outCPtab)
            out_cand_tab = export_candidates(self.H_species, self.non_H_species, self.disregard)
            write_to_TSV(candname, out_cand_tab)
        if rec and self.filled:
            write_to_TSV(recname, outrecording)
            
    def pulse_integrate(self, pulse_start = None, pulse_stop = None):
        if pulse_start == None:
            try:
                pulse_start = self.def_val['pulse_start']
            except KeyError:
                pulse_start = 0
        if pulse_stop == None:
            try:
                pulse_stop = self.columns[self.time_col][-1]
            except AttributeError:
                pulse_stop = pulse_start
        else:
            pulse_stop = self.def_val['pulse_stop']
        
        # Pulse integration of the raw data
        # Datapoints prior to pulse_start will be averaged as background
        if self.filled:
            self.pulse_integral = {}

            bgr_range = (self.columns[self.time_col] < pulse_start)
            if list(bgr_range).count(True) > 1:
                do_bgr = True
                self.background = {}
            else:
                do_bgr = False
            int_range = (self.columns[self.time_col] > pulse_start) * (self.columns[self.time_col] < pulse_stop)
            for mass in self.header_int:
                tempdat = self.columns[mass]
                if do_bgr:
                    bgr = sp.mean(tempdat[bgr_range])
                    self.background[mass] = bgr
                else:
                    bgr = 0
                int_val = sp.trapz(tempdat[int_range] - bgr, self.columns[self.time_col][int_range])
                self.pulse_integral[mass] = int_val
        
        # Pulse integration of the deconvoluted data
        if self.deconvoluted:
            self.pulse_integrated_results = {}
            bgr_range = (self.rescols[self.time_col] < pulse_start)
            if list(bgr_range).count(True) > 1:
                do_bgr = True
                self.background_results = {}
            else:
                do_bgr = False
            int_range = (self.rescols[self.time_col] > pulse_start) * (self.rescols[self.time_col] < pulse_stop)
            
            for gas in self.H_species.keys():                
                temp_pres = self.rescols[gas]['pressure']
                temp_rat = self.rescols[gas]['ratio']
                if do_bgr:
                    bgr_p = sp.mean(temp_pres)
                    bgr_r = sp.trapz(temp_pres[bgr_range] * temp_rat[bgr_range], self.rescols[self.time_col][bgr_range]) /     sp.trapz(temp_rat[bgr_range], self.rescols[self.time_col][bgr_range])
                    self.background_results[gas] = {'pressure': bgr_p, 'ratio': bgr_r}
                else:
                    bgr_p = 0
                int_val_p = sp.trapz(temp_pres[int_range] - bgr_p, self.rescols[self.time_col][int_range])
                int_val_r = sp.trapz(temp_pres[int_range] * temp_rat[int_range], self.rescols[self.time_col][int_range]) / sp.trapz(temp_pres[int_range], self.rescols[self.time_col][int_range])
                self.pulse_integrated_results[gas] = {'pressure': int_val_p, 'ratio': int_val_r}
            
            for gas in self.non_H_species.keys():                
                temp_pres = self.rescols[gas]['pressure']
                if do_bgr:
                    bgr_p = sp.mean(temp_pres)
                    self.background_results[gas] = {'pressure': bgr_p}
                else:
                    bgr_p = 0
                int_val_p = sp.trapz(temp_pres[int_range] - bgr_p, self.rescols[self.time_col][int_range])
                self.pulse_integrated_results[gas] = {'pressure': int_val_p}


class Profile:
    
    def __init__(self, mass_col, intensity_col, tag, default=None):
        
        self.type = "Profile"
        self.mass_col = sp.array(mass_col)
        self.intensity_col = sp.array(intensity_col)
        self.filled = False
        self.deconvoluted = False
        self.calibrated = False
        if type(tag) != dict:
            self.tag = {'raw_tag':tag}
        else:
            self.tag = tag
            
        if default == None:
            self.def_val = load_default_values({})
        else:
            self.def_val = load_default_values(default)
        # Identify recorded masses and extract the maximal intensities
        try:
            first_mass = max(1, int(min(self.mass_col)))
        except:
            first_mass = 0
        try:
            last_mass = int(max(self.mass_col))
        except:
            last_mass = 0
        
        self.header_int = []
        self.recorded = sp.zeros(last_mass + 1)
        self.MID_col = sp.zeros(last_mass + 1)
        for mass in range(first_mass, last_mass + 1):
            temp_area = (self.mass_col > mass - 0.5) * (self.mass_col < mass + 0.5)
            if list(temp_area).count(True) > 0:
                self.recorded[mass] = 1
                self.header_int.append(mass)
                intensity = max(self.intensity_col[temp_area])
                self.MID_col[mass] = intensity
                
        self.header_int = sp.array(sorted(self.header_int))
        
        if len(self.header_int) > 0:
            self.filled = True
        # ime objekta - ce ni v tag-u, se avtoimenuje
        if 'title' not in self.tag.keys():
            self.tag['title'] = self.def_val['profile_name']
            
        if self.filled:
            # Initialize mass-space and cracking patterns
            self.molecules = molecules2.mass_space(self.header_int[-1])
            #TO DO - by default, initialize cracking patterns with appropriate calibration files
            self.molecules.init_CP()
            
    def replace_CP(self, path):
        self.molecules.init_CP(path)


    def make_candidates(self):
        
        self.candidates_dict = make_candidates(self.molecules, self.H_species, self.non_H_species)
        
    def deconvolute(self, H_species, non_H_species, disregard):
        self.H_species = H_species
        self.non_H_species = non_H_species
        self.disregard = disregard
        
        self.make_candidates()
        self.parameters = self.candidates_dict['parameters']
        self.pressures = self.candidates_dict['pressures']
        self.ratios = self.candidates_dict['ratios']
        self.masses_of_interest = []
        for mass in range(1,self.header_int[-1] + 1):
            if sum(self.candidates_dict['candidates'])[mass] != 0:
                self.masses_of_interest.append(mass)
        
        
        results, sim_masses = fit_line(self.MID_col, self.recorded, self.header_int, self.candidates_dict, disregard)
        
        resline = {}
        for prs in self.pressures:
            resline[prs] = results[prs]
        resline['residual'] = results['residual_value']
        for ratio in self.ratios:
            resline[ratio] = results[ratio]
        resline['duration'] = results['duration']

        self.resline = {}
        self.resline['residual'] = resline['residual']
 
        for key in self.H_species.keys():
            self.resline[key]={}
            self.resline[key]['pressure'] = resline[self.H_species[key][0]]
            if type(self.H_species[key][2]) == list:
                self.resline[key]['ratio'] = resline[self.H_species[key][2][0]]
            elif type(self.H_species[key][2]) == str:
                self.resline[key]['ratio'] = resline[self.H_species[key][2]]
            elif type(self.H_species[key][2] in [int, float]):
                self.resline[key]['ratio'] = self.H_species[key][2]
        for key in self.non_H_species.keys():
            self.resline[key]={}
            self.resline[key]['pressure'] = resline[self.non_H_species[key][0]]
            
      
        self.sim_MID_col = [sim_masses[mass] for mass in self.masses_of_interest]

        self.devoncoluted = True
        
    def calibrate(self, molecule, ratio_def, peak_defs, disregard):
        
        self.calib_tag = {}
        if type(molecule) == str:
            try:
                mol_def_1 = molecules2.H_molecules_d[molecule]
                rel_int = self.molecules.calib[mol_def_1[0]][-1]
                mol_def = {'NH_mass': mol_def_1[2], 'nAt': mol_def_1[1], 'rel_int': rel_int}
            except KeyError:
                #print "Unknown molecule"
                pass
        if type(molecule) == dict:
            if not set(['NH_mass', 'nAt', 'rel_int']).issubset(set(mol_def.keys())):
                # print invalid molecule definition
                pass
            
        self.calib_candidates_dict = make_calibration_candidates(self.molecules, mol_def, ratio_def, peak_defs)        
        self.calib_masses_of_interest = []
        for mass in range(1,self.header_int[-1] + 1):
            if sum(self.calib_candidates_dict['candidates'])[mass] != 0:
                self.calib_masses_of_interest.append(mass)
        
        results, sim_masses = fit_line(self.MID_col, self.recorded, self.header_int, self.calib_candidates_dict, disregard)
        
        self.calib_line = {}
        self.calib_line['pressure'] = results['pres']
        self.calib_line['duration'] = results['duration']
        self.calib_line['residual'] = results['residual_value']
        if type(ratio_def) in [float, int]:
            self.calib_line['ratio'] = ratio_def
        else:
            self.calib_line['ratio'] = results[self.calib_candidates_dict['ratios'][0]]
        self.calib_line["peaks"] = {}
        for i in range(len(peak_defs)):
            peak_def = peak_defs[i]
            level = i + 1
            if type(peak_def) in [float, int]:
                self.calib_line['peaks'][level] = peak_def
            else:
                if type(peak_def) == list:
                    peak_name = peak_def[0]
                else:
                    peak_name = peak_def
                self.calib_line['peaks'][level] = results[peak_name]
        
           
      
        self.sim_MID_col = [sim_masses[mass] for mass in self.calib_masses_of_interest]
        calib.tag['title'] = "Calibration for %s" %molecule
        self.calibrated = True

      

    def export(self, name=None, write_path=None, rec=False):
    
        
        # Contruction of data for export to a TSV file
        
        outlist = []
        outheader = []
        
        if rec and self.filled:
            outheader.append("Mass [AMU]")
            outheader.append("Intensity [arb.u.]")        
            outlist.append(self.mass_col)
            outlist.append(self.intensity_col)
            
        if self.deconvoluted:
            outheader = outheader + ["Mass [AMU]","Simulated intensity [arb.u.]","Gas","Pressure [arb.u.]","Gas","H/(D+H)"]
            outlist.append(self.masses_of_interest)
            outlist.append(self.sim_MID_col)

            prescol_legend = []
            prescol_value = []
            for gas in self.H_species.keys() + self.non_H_species.keys():
                prescol_legend.append(gas)
                prescol_value.append(self.resline[gas]['pressure'])
            ratiocol_legend = []
            ratiocol_value = []
            for gas in self.H_species.keys():
                ratiocol_legend.append(gas)
                ratiocol_value.append(self.resline[gas]['ratio'])
            prescol_legend.append('residual')
            prescol_value.append(self.resline['residual'])

            outlist.append(prescol_legend)
            outlist.append(prescol_value)
            outlist.append(ratiocol_legend)
            outlist.append(ratiocol_value)
       
        blocklength = max(map(len, outlist))
        outlist = map(list, outlist)
        for i in range(0,len(outlist)):
            outlist[i] = PadRight(outlist[i], blocklength, "")
        self.outtab = sp.transpose(sp.array(outlist))
        self.outtab = [outheader] + list(self.outtab)
        
        if write_path == None:
            write_path = self.def_val['export_dir']
        if name == None:
            name = self.tag['title']
        if not path.exists(write_path):
            mkdir(write_path)

        out_path = path.join(export_dir, "%s_PROFILE.txt" % name)
        write_to_TSV(out_path, self.outtab)    
        if self.deconvoluted:
            cpname = path.join(export_dir, "%s_PROFILE_CALIBRATION.txt" %name)
            candname = path.join(export_dir, "%s_PROFILE_PARAMETERS.txt" %name)
            outCPtab = export_CP(self.molecules.calib)
            write_to_TSV(cpname, outCPtab)
            out_cand_tab = export_candidates(self.H_species, self.non_H_species, self.disregard)
            write_to_TSV(candname, out_cand_tab)


# zdruzevanje rezultatov fita v enem trejsu

def join_traces(in_trace_list):

    tracelist = []
    for trace in in_trace_list:
        # najprej sortiranje tracov po zacetnem casu
        tracelist.append([trace.rescols[trace.time_col][0], trace, trace.time_col])
    tracelist = sorted(tracelist)    
    tc = tracelist[0][2]
    for entry in tracelist[1:]:
        if entry[2] != tc:
            print "Time col not same in all traces."

    rescols_j = {}
    HM_j = {}
    NHM_j = {}
    moi_j = []
    stc_j = {}

    rescols_j[tc] = sp.array([])
    rescols_j['residual'] = sp.array([])
    rescols_j['duration'] = sp.array([])


    for entry in tracelist:
        trace = entry[1]
        rc_t = trace.rescols
        stc_t = trace.simtracecol

        HM_t = trace.H_species
        NHM_t = trace.non_H_species
        moi_t = trace.masses_of_interest

        zerocol_j = sp.zeros(len(rescols_j[tc]))
        zerocol_t = sp.zeros(len(rc_t[tc]))

        print "trace starts at: %s, H molecules: %s, non-H molecules: %s" %(entry[0], ", ".join(HM_t.keys()), ", ".join(NHM_t.keys()))

        for gas in HM_j.keys():
            if gas not in HM_t.keys():
                HM_t[gas] = HM_j[gas]
                rc_t[gas] = {'pressure': zerocol_t, 'ratio': zerocol_t}
        for gas in NHM_j.keys():
            if gas not in NHM_t.keys():
                NHM_t[gas] = NHM_j[gas]
                rc_d[gas] = {'pressure': zerocol_t}

        for gas in HM_t.keys():
            if gas not in HM_j.keys():
                HM_j[gas] = HM_t[gas]
                rescols_j[gas] = {'pressure': zerocol_j, 'ratio': zerocol_j}
            for what in ['pressure', 'ratio']:
                rescols_j[gas][what] = sp.concatenate((rescols_j[gas][what], rc_t[gas][what]))
        for gas in NHM_t.keys():
            if gas not in NHM_j.keys():
                NHM_j[gas] = NHM_t[gas]
                rescols_j[gas] = {'pressure': zerocol_j}
            rescols_j[gas]['pressure'] = sp.concatenate((rescols_j[gas]['pressure'], rc_t[gas]['pressure']))
        for what in ['duration', 'residual', tc]:
            rescols_j[what] = sp.concatenate((rescols_j[what], rc_t[what]))

        for mass in moi_t:
            if mass not in moi_j:
                moi_j.append(mass)
        for mass in stc_j.keys():
            if mass not in stc_t.keys():
                stc_t[mass] = zerocol_t
        for mass in stc_t.keys():
            if mass not in stc_j.keys():
                stc_j[mass] = zerocol_j
            stc_j[mass] = sp.concatenate((stc_j[mass], stc_t[mass]))

    jtrace = copy(tracelist[0][1])
    jtrace.rescols = rescols_j
    jtrace.H_species = HM_j
    jtrace.non_H_species = NHM_j
    jtrace.masses_of_interest = moi_j
    jtrace.simtracecol = stc_j
    
    return jtrace
    
    
# zdruzevanje kalibracijskih rezultatov

def join_traces_calib(in_trace_list):

    tracelist = []
    for trace in in_trace_list:
        # najprej sortiranje tracov po zacetnem casu
        tracelist.append([trace.calibcols[trace.time_col][0], trace, trace.time_col])
    tracelist = sorted(tracelist)    
    tc = tracelist[0][2]
    for entry in tracelist[1:]:
        if entry[2] != tc:
            print "Time col not same in all traces."

    calibcols_j = {}
    peaks_j = []
    moi_j = []
    stc_j = {}

    calibcols_j[tc] = sp.array([])
    calibcols_j['residual'] = sp.array([])
    calibcols_j['duration'] = sp.array([])
    calibcols_j['pressure'] = sp.array([])
    calibcols_j['ratio'] = sp.array([])
    calibcols_j['peaks'] = {}


    for entry in tracelist:
        trace = entry[1]
        cc_t = trace.calibcols
        stc_t = trace.simtracecol

        peaks_t = trace.calibcols['peaks'].keys()
        moi_t = trace.calib_masses_of_interest

        zerocol_j = sp.zeros(len(calibcols_j[tc]))
        zerocol_t = sp.zeros(len(cc_t[tc]))

        print "trace starts at: %s, peaks: %s" %(entry[0], ", ".join(map(str, peaks_t)))

        for peak in peaks_j:
            if peak not in peaks_t:
                peaks_t.append(peak)
                cc_t['peaks'][peak] = zerocol_t


        for peak in peaks_t:
            if peak not in peaks_j:
                peaks_j.append(peak)
                calibcols_j['peaks'][peak] = zerocol_j
                
        for peak in peaks_t:
            if peak not in peaks_j:
                simcols_j['peaks'][peak] = zerocol_j
            calibcols_j['peaks'][peak] = sp.concatenate((calibcols_j['peaks'][peak], cc_t['peaks'][peak]))

        for what in ['duration', 'residual', 'pressure', 'ratio', tc]:
            calibcols_j[what] = sp.concatenate((calibcols_j[what], cc_t[what]))

        for mass in moi_t:
            if mass not in moi_j:
                moi_j.append(mass)
        for mass in stc_j.keys():
            if mass not in stc_t.keys():
                stc_t[mass] = zerocol_t
        for mass in stc_t.keys():
            if mass not in stc_j.keys():
                stc_j[mass] = zerocol_j
            stc_j[mass] = sp.concatenate((stc_j[mass], stc_t[mass]))

    jtrace = copy(tracelist[0][1])
    jtrace.calibcols = calibcols_j
    jtrace.calib_masses_of_interest = moi_j
    jtrace.simtracecol = stc_j
    
    return jtrace
