# RGA data container class
# Two RGA data containter classes: profile (intensity vs mass) and time trace (intensity vs time at specific masses)
# Basic data-processing functions are included in this module
# Mass-space and cracking pattern functions are imported from molecules2
# Fitting imported from RGA_fitting

# Aug 2016
# Aleksander Drenik
# IJS, Ljubljana
# IPP, Garching


version = '2.1'
versions = {'Class': version}

from . import molecules3 as molecules2
import scipy as sp
import random
from copy import deepcopy
import time
from .RGA_calibration_reader import write_to_TSV

try:
    #from RGA_fitting import version, check_candidates, make_candidates, check_disregard, export_CP, export_candidates, make_calibration_candidates, fit_line
    from .RGA_fitting import version, RGA_fitting, export_CP, export_candidates
    versions['RGA_fitting'] = version
    fitting_loaded = True
except ImportError:
    versions['RGA_fitting'] = version
    fitting_loaded = False
    
from os import path, mkdir

versions['molecules'] = molecules2.version

known_time_labels = ['time', 'hours'] # labels of columns which are automatically recognized as time columns in the Trace class

# Import the default values

def load_default_values(new_defaults=None):
    """
    Complete default dict
    :param new_defaults:    (optional) dict with default values for 'trace_name', 'profile_name', 'export_dir', 'pulse_start' and 'pulse_stop'
    :return:                input dict completed with standard values
    """
    default = {'trace_name': 'Timetrace recording',
               'profile_name': "Profile recording",
               'export_dir': 'EXPORT',
               'pulse_start': 0,
               'pulse_stop': None,
               }
    if isinstance(new_defaults, dict):
        default.update(new_defaults)

    return default

def pin_point(input_list, sought_value):
    """Return the index of the element in input_list closest to sought_value"""
    input_col = sp.array(input_list)
    temp_list = list(map(abs, input_col - sought_value))
    nearest_index = temp_list.index(min(temp_list))
    return nearest_index


def PadRight(inputlist, TargetLength, FillValue):
    """Append FillValue to inputlist until length(inputlist) == TargetLength"""
    if len(inputlist) < TargetLength:
        return inputlist + [FillValue] * (TargetLength - len(inputlist))
    else:
        return inputlist

# TO DO: error handling - raise warning

def perturb_CP(container, perturbation):
    CP_new = deepcopy(container.molecules.calib)
    perturb_d = {}
    if type(perturbation) in [int, float]:
        # construct the dictionary of perturbation for each molecule in the CP
        for what in ['H', 'non-H']:
            for mol in CP_new[what].keys():
                perturb_d[mol] = perturbation
        CP_new['version'] = "%s perturbed by %s" %(CP_new['version'], perturbation)
    elif type(perturbation) == dict:
        # take only the keys that appear in the CP
        for what in ['H', 'non-H']:
            for mol in CP_new[what].keys():
                try:
                    perturb_d[mol] = perturbation[mol]
                except KeyError:
                    perturb_d[mol] = 0
        CP_new['version'] = "%s perturbed by various degrees" %(CP_new['version'])
    else:
        # invalid perturbation definition, perturb all by 0
        for what in ['H', 'non-H']:
            for mol in CP_new[what].keys():
                perturb_d[mol] = 0
        CP_new['version'] = "%s not perturbed - invalid definition of perturbation" %(CP_new['version'])
                
    for mol in CP_new['H'].keys():
        mol_d = CP_new['H'][mol]
        pks = mol_d['peaks']
        for i in range(1, len(pks)):
            pks[i] = pks[i] * (1 + 2 * perturb_d[mol] * (random.random() - 0.5))
    for mol in CP_new['non-H'].keys():
        mol_d = CP_new['non-H'][mol]
        pks = mol_d['peaks']
        for i in range(1, len(pks)):
            pks[i][1] = pks[i][1] * (1 + 2 * perturb_d[mol] * (random.random() - 0.5))
            
    container.replace_CP(CP_new)
    

# old version of perturb_CP
# kept for legacy reasons and to be removed
def perturb_CP_old(container, perturbation):
    """Perturb the cracking pattern loaded in the container with the magnitude of perturbation"""
    CP_new = container.molecules.calib
    for key in CP_new.keys():
        entry = CP_new[key]
        if type(entry) == list:
            for i in [0, 1]:
                entry[i] = entry[i] * (1 + 2 * perturbation * (random.random() - 0.5))
    CP_new['version'] = "%s perturbed by %s" %(CP_new['version'], perturbation)
    container.replace_CP(CP_new)

def read_molecule_definition(molecules, molecule):
    # reads the molecule definition and returns number of H atoms, total mass of non H atoms and the relative intensity
    # the molecule name, if known
    # and the error
    # error codes
    # 5001 - uknown specified molecule
    # 5002 - invalid specified molecule definition
    #self.calib_tag = {}
    err = 0
    mol_name = 'unknown'
    if type(molecule) == str:
        # TO DO move calib_tag to container
        #self.calib_tag['title'] = "Calibration of %s for %s" %(self.tag['title'], molecule)
        # try to read the properties of the molecules from the loaded calibration
        try:
            mol_def = molecules.calib['H'][molecule]
            mol_name = molecule
        except KeyError:
            err = 5001
            mol_def = {'non-H-mass':0, 'H-atoms':0, 'intensity':0}
    if type(molecule) == dict:
        mol_def = molecule
        #self.calib_tag['title'] = "Calibration of %s" %self.tag['title']
        if not set(['non-H-mass', 'H-atoms', 'intensity']).issubset(set(mol_def.keys())):
            # print invalid molecule definition
            err = 5002
            mol_def = {'non-H-mass':0, 'H-atoms':0, 'intensity':0}
            
    return err, {'mol_def': mol_def, 'molecule': mol_name}


class RGA_base_container:
    # Base class for RGA data
    # provides initializing cracking patterns, preparation of fitting parameters
    # etc

    def __init__(self):
        # Initialize mass-space and cracking patterns
        self.molecules = molecules2.mass_space(self.header_int[-1])
        #TO DO - by default, initialize cracking patterns with appropriate calibration files
        self.molecules.init_CP()
        self.make_recorded()
        
    def make_recorded(self): # not Trace specific
        # make the 'recorded' array
        # an array the lengh of the mass_space dimension
        # at each mass index, element = 1 if there are recordings at this mass, 0 otherwise
        try:
            max_mass = self.molecules.max_mass
        except AttributeError:
            max_mass = max(self.header_int)
        self.recorded = sp.zeros(max_mass + 1)
        for i in range(max(self.header_int) + 1):
            if i in self.header_int:
                self.recorded[i] = 1
            
    def replace_CP(self, path): # Not Trace specific
        """Replace the cracking patterns. path can be string with the calibration file location or calibration dictionary"""
        if self.filled:
            self.molecules.init_CP(path)
            self.make_recorded()


class Trace(RGA_base_container, RGA_fitting):
    # RGA data container for time trace shaped data
    # Data is provided as a dictionary in 'columns'
    # column names: if in known_time_labels the column will be assigned as a time column
    # column names: integers - will be assigned as intensity recordings and added to header_int
    # Other information about the data is provided in the tag (as dictionary)
    # Internal errors: tag['init_errors']
    # 101: empty or invalid data dictionary
    # 102: empty data columns
    # 103: empty header_int
    
    def __init__(self, columns, tag, default=None):
        """Create a Trace object and fill it with data
        Parameters:
        :param columns: dict of 1D data arrays. column names either:
                            "time" (or any other value in known_time_labels)    to be automatically recognized as time vector
                            mass (int)                                          to be automatically recognized as data corresponding to given mass
                            arbitrary string                                    to be manually assigned to timetrace via self.set_timecol
        :param tag: dict of metadata
                            time_unit   unit of time vector
                            title       set data title
        :param default: dict of default values (see set_default_values for understood options)
        :return None
        """
        
        self.type = "Trace"
        self.columns = {}
        self.filled = False # Trace contains data (Boolean)
        self.deconvoluted = False # deconvolute function has been used on the data
        self.calibrated = False # calibrate function has been used on the data
        bloc_len = 0 # length of shortest data column
        self.echo = True
        
        # Load the default values for the object
        if default == None:
            self.def_val = load_default_values({})
        else:
            self.def_val = load_default_values(default)
        
        if type(tag) != dict: # if tag is not provided in the correct format
            self.tag = {'raw_tag':tag} # store it anyway
        else:
            self.tag = tag
        self.tag['init_errors'] = [] # Errors encountered in Trace.__init__
        header_int = [] # List of recorded masses, represented by integers
        lengths = {} # Lengths of all of the columns
        for key in columns.keys():
            # Check if all columns have the same length
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
            bloc_len = min(lengths.values()) # shortest column length
            for key in self.columns.keys():
                # If a column is longer than the shortest column (bloc_len)
                # drop any datapoints beyond bloc_len
                self.columns[key] = self.columns[key][:bloc_len]
            if bloc_len == 0: # at least one zero-lenght column was included in the data
                self.tag['init_errors'].append(102)
        except ValueError:
            self.tag['init_errors'].append(101)
        # Create the column of indices of datapoints
        self.columns['index'] = sp.arange(bloc_len)
        if len(header_int) > 0 and bloc_len > 0:
            self.filled = True # Trace contains technically OK data
            # set the time column                
            time_col_candidates = list(set(self.columns.keys()) & set(known_time_labels))
            try:
                self.set_timecol(time_col_candidates[0])
            except:
                # if no label is recognized, set the index column as the time column
                self.set_timecol('index')

            self.header_int = sp.array(sorted(header_int))
            # Title of the Trace object
            # If 'title' is not provided in the tag, use default value
            if 'title' not in self.tag.keys():
                self.tag['title'] = self.def_val['trace_name']
            # Initialize mass-space and cracking patterns
            # This has now been moved to RGA_base_container

            RGA_base_container.__init__(self)
            
    def set_timecol(self, time_key):
        """Designate one of the columns as the time column"""
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
        # testno - mase morajo iti do zadnje mase v masses_of_interest
        # last_mass = max(self.masses_of_interest)  
        # testno 2 - max maso vzami iz mass_space
        last_mass = self.molecules.max_mass      
        #for mass in range(1,self.header_int[-1]+1):
        for mass in range(1, last_mass + 1):
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
    # From the ti-th line in the results time-trace
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
        """Create a Profile object from a timeslice of the data, at time closest to wanted"""
        ti = pin_point(self.columns[self.time_col], wanted)
        #line = self.make_line(ti)
        actual = self.columns[self.time_col][ti]
        MID_col = self.make_MID_col(ti)
        datatag = deepcopy(self.tag)
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
    
    # now in RGA_fitting
    #def make_candidates(self): # Not Trace specific
    #    """Make data fitting candidates"""
    #    self.candidates_dict = make_candidates(self.molecules, self.H_species, self.non_H_species)
    
    # TO DO: keep track of make_recorded    
    def deconvolute(self, H_species, non_H_species, disregard, start_time, stop_time, step, n_iter=0):
        """Deconvolute the data with H_species and non_H_species"""
        self.H_species = H_species # Not Trace specific
        self.non_H_species = non_H_species # Not Trace specific
        self.disregard = disregard # Not Trace specific
        
        self.make_candidates() # Now from RGA_fitting

        # construction of the recorded array
        # now based on the masses_of_interest        
        
        self.recoded = self.make_recorded() 
        # Trace specific part starts here
        ti_list = self.columns['index'][(start_time <= self.columns[self.time_col]) * (stop_time >= self.columns[self.time_col])] 
        if len(ti_list) == 0:
            self.glob_duration = 0
            return
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

            results, sim_masses = self.fit_line(line, n_iter=n_iter) # now from RGA_fitting
            
            outcols[self.time_col].append(line[0])
            for prs in self.pressures:
                outcols[prs].append(results[prs])
            outcols['residual'].append(results['residual_value'])
            for ratio in self.ratios:
                outcols[ratio].append(results[ratio])
            outcols['duration'].append(results['duration'])

            simtrace.append(sim_masses)
        #return outcols    
        self.glob_duration = time.clock() - glob_start
        
        self.rescols = {}
        self.rescols['index'] = sp.array(ti_list)
        self.rescols[self.time_col] = sp.array(outcols[self.time_col])
        self.rescols['residual'] = sp.array(outcols['residual'])
        self.rescols['duration'] = sp.array(outcols['duration'])
        tablength = len(outcols[self.time_col])
        for key in self.H_species.keys():
            self.rescols[key]={}
            if type(self.H_species[key][0]) == list:
                pres_name = self.H_species[key][0][0]
            elif type(self.H_species[key][0]) == str:
                pres_name = self.H_species[key][0]
            self.rescols[key]['pressure'] = sp.array(outcols[pres_name])
            if type(self.H_species[key][2]) == list:
                self.rescols[key]['ratio'] = sp.array(outcols[self.H_species[key][2][0]])
            elif type(self.H_species[key][2]) == str:
                self.rescols[key]['ratio'] = sp.array(outcols[self.H_species[key][2]])
            elif type(self.H_species[key][2] in [int, float]):
                self.rescols[key]['ratio'] = self.H_species[key][2] * sp.ones(tablength)
        for key in self.non_H_species.keys():
            self.rescols[key]={}
            pres_def = self.non_H_species[key][0]
            if type(pres_def) == list:
                pres_name = pres_def[0]
            elif type(pres_def) == str:
                pres_name = pres_def
            self.rescols[key]['pressure'] = sp.array(outcols[pres_name])
            
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
        if self.echo:
            # TO DO: print -> report: Fitting {title} for {candidates}, duration, errors
            print("Deconvolutiuon for %s done, duration %s s." %(self.tag['title'], self.glob_duration))
        
    def calibrate(self, molecule, ratio_def, peak_defs, disregard, start_time, stop_time, step, n_iter=0):
        
        self.calib_tag = {}
        err, mol_read = read_molecule_definition(self.molecules, molecule)
        self.calib_tag['molecule_read_error'] = err
        mol_def = mol_read['mol_def']
        mol_name = mol_read['molecule']
        if mol_name != 'unknown':
            self.calib_tag['title'] = "Calibration of %s for %s" %(self.tag['title'], mol_name)
        else:
            self.calib_tag['title'] = "Calibration of %s" %self.tag['title']
        self.calib_tag = {}
        # pass the molecule definition to the make_calibration_candidates function
        self.make_calibration_candidates(mol_def, ratio_def, peak_defs)
        # make calibration candidates up to here
        
        ti_list = self.columns['index'][(start_time <= self.columns[self.time_col]) * (stop_time >= self.columns[self.time_col])] 
        if len(ti_list) == 0:
            self.glob_duration = 0
            return
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
            
            #results, sim_masses = fit_line(line, self.recorded, self.header_int, self.calib_candidates_dict, disregard, n_iter=n_iter)
            results, sim_masses = self.fit_line(line, n_iter=n_iter)
            outcols[self.time_col].append(line[0])
            
            for parameter in self.calib_candidates_dict['parameters']:
                outcols[parameter].append(results[parameter])
            
            outcols['duration'].append(results['duration'])
            outcols['residual'].append(results['residual_value'])

            simtrace.append(sim_masses)
            
        self.calib_glob_duration = time.clock() - glob_start
        
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
                    
    def correct_offset(self, offset_mass):
        """Correction of the offset, defined by the offset signal
        Args: offset_mass - integer or list of integers
        If list of integers is given, the offset signal is the average signal value from specified masses.
        This function overrides the original columns"""
        
        if type(offset_mass) == int:
            offset_mass = [offset_mass]
        self.offset_mass = list(set(offset_mass) & set(self.header_int))
        if len(self.offset_mass) > 0:
            self.columns['offset'] = sp.zeros(len(self.columns['index']))
            for mass in offset_mass:
                self.columns['offset'] += self.columns[mass]
            self.columns['offset'] *= 1.0 / len(self.offset_mass)
        
            for mass in self.header_int:
                self.columns[mass] -= self.columns['offset']
            
    def clip_negative(self, limit=0):
        """Wrapper for the .clip function for all mass columns"""
        for mass in self.header_int:
            self.columns[mass] = self.columns[mass].clip(limit)
  
    
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
            
    def pulse_integrate(self, pulse_start = None, pulse_stop = None, remove_bgr = True):
        if pulse_start == None:
            try:
                pulse_start = self.def_val['pulse_start']
            except KeyError:
                pulse_start = 0
        if pulse_stop == None:
            pulse_stop = self.def_val['pulse_stop']
        if pulse_stop == None:   
            try:
                pulse_stop = self.columns[self.time_col][-1]
            except AttributeError:
                pulse_stop = pulse_start
            
            
        
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
                if not remove_bgr:
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
                
    def join_isotopes(self, isotope_list, common):
        #if common not in self.rescols.keys():
        p_col = sp.zeros(len(self.rescols['index']))
        r_col = sp.zeros(len(self.rescols['index']))
        in_rescols = list(set(self.rescols.keys()) & set(isotope_list))
        if len(in_rescols) > 0:
            self.rescols[common] = {}
        for isotope in in_rescols:
            p_col = p_col + self.rescols[isotope]['pressure']
        self.rescols[common]['pressure'] = p_col
        try:
            for isotope in in_rescols:
                r_col = r_col + self.rescols[isotope]['ratio'] * self.rescols[isotope]['pressure'] / p_col
            self.rescols[common]['ratio'] = r_col
        except KeyError:
            pass


class Profile(RGA_base_container, RGA_fitting):
    
    def __init__(self, mass_col, intensity_col, tag, default=None):
        """Create a Profile object and fill it with data
        Parameters:
        :param massc_col: mass column (list or array)
        :param intensity_col: intensity column (list or array), must have same dimension as mass column
        :param tag: dict of metadata
                            time_unit   unit of time vector
                            title       set data title
        :param default: dict of default values (see set_default_values for understood options)
        :return None
        """
        
        
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
        #self.recorded = sp.zeros(last_mass + 1)
        self.MID_col = sp.zeros(last_mass + 1)
        for mass in range(first_mass, last_mass + 1):
            temp_area = (self.mass_col > mass - 0.5) * (self.mass_col < mass + 0.5)
            if list(temp_area).count(True) > 0:
                #self.recorded[mass] = 1
                self.header_int.append(mass)
                intensity = max(self.intensity_col[temp_area])
                self.MID_col[mass] = intensity
                
        self.header_int = sp.array(sorted(self.header_int))
        
        if len(self.header_int) > 0:
            self.filled = True
        # Title of the object. If not in tag, then name automatically with the default title
        if 'title' not in self.tag.keys():
            self.tag['title'] = self.def_val['profile_name']
            
        if self.filled:
            RGA_base_container.__init__(self)
                   
    def deconvolute(self, H_species, non_H_species, disregard, n_iter=0):
        self.H_species = H_species
        self.non_H_species = non_H_species
        self.disregard = disregard
        
        self.make_candidates() # Now from RGA_fitting
        results, sim_masses = self.fit_line(self.MID_col, n_iter=n_iter)
        
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
        
    def calibrate(self, molecule, ratio_def, peak_defs, disregard, n_iter=0):
        
        self.calib_tag = {}
        err, mol_read = read_molecule_definition(self.molecules, molecule)
        self.calib_tag['molecule_read_error'] = err
        mol_def = mol_read['mol_def']
        mol_name = mol_read['molecule']
        if mol_name != 'unknown':
            self.calib_tag['title'] = "Calibration of %s for %s" %(self.tag['title'], mol_name)
        else:
            self.calib_tag['title'] = "Calibration of %s" %self.tag['title']
            
        # pass the molecule definition to the make_calibration_candidates function
        self.make_calibration_candidates(mol_def, ratio_def, peak_defs)
        
        #results, sim_masses = fit_line(self.MID_col, self.recorded, self.header_int, self.calib_candidates_dict, disregard, n_iter=n_iter)
        results, sim_masses = self.fit_line(line, n_iter=n_iter)
        
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

def join_traces(in_trace_list, echo=False):
    # error messages:
    # 101 - Time column not the same in all traces
    err = 0
    tracelist = []
    parameters = []
    for trace in in_trace_list:
        # najprej sortiranje tracov po zacetnem casu
        tracelist.append([trace.rescols[trace.time_col][0], trace, trace.time_col])
    tracelist = sorted(tracelist)    
    tc = tracelist[0][2]
    for entry in tracelist[1:]:
        if entry[2] != tc:
            err = 101
            print("Time col not same in all traces.")

    rescols_j = {}
    HM_j = {}
    NHM_j = {}
    moi_j = []
    stc_j = {}

    rescols_j[tc] = sp.array([])
    rescols_j['residual'] = sp.array([])
    rescols_j['duration'] = sp.array([])


    for entry in tracelist:
        sub_param = {}
        start_t = entry[0]
        sub_param['start time'] = start_t
        trace = entry[1]
        rc_t = trace.rescols
        stc_t = trace.simtracecol

        HM_t = trace.H_species
        NHM_t = trace.non_H_species
        moi_t = trace.masses_of_interest
        sub_param['H_molecules'] = HM_t
        sub_param['non_H_molecules'] = NHM_t
        
        parameters.append(sub_param)

        zerocol_j = sp.zeros(len(rescols_j[tc]))
        zerocol_t = sp.zeros(len(rc_t[tc]))
        
        if echo:
            print( "trace starts at: %s, H molecules: %s, non-H molecules: %s" %(start_t, ", ".join(HM_t.keys()), ", ".join(NHM_t.keys())))

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

    jtrace = deepcopy(tracelist[0][1])
    jtrace.rescols = rescols_j
    jtrace.H_species = HM_j
    jtrace.non_H_species = NHM_j
    jtrace.masses_of_interest = moi_j
    jtrace.simtracecol = stc_j
    jtrace.joined_parameters = parameters
    
    return err, jtrace
    
    
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
            print( "Time col not same in all traces.")

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

        print("trace starts at: %s, peaks: %s" %(entry[0], ", ".join(map(str, peaks_t))))

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

    jtrace = deepcopy(tracelist[0][1])
    jtrace.calibcols = calibcols_j
    jtrace.calib_masses_of_interest = moi_j
    jtrace.simtracecol = stc_j
    
    return jtrace
    
    
# Construction of a common molecule signals from different non-H isotopes
# such as a commomn ammonia signal from 14N- and 15N ammonia signals
# input:
# isotope list: [['ammonia14', 'ammonia15']]
# common: 'ammonia'

# Moved to Trace as method
# TO DO: consider a universal version to work both on rescols and resline
"""
def total_isotopes_trace(trace, isotope_list, common):
    #if common not in trace.rescols.keys():
    p_col = sp.zeros(len(trace.rescols['index']))
    r_col = sp.zeros(len(trace.rescols['index']))
    in_rescols = list(set(trace.rescols.keys()) & set(isotope_list))
    if len(in_rescols) > 0:
        trace.rescols[common] = {}
    for isotope in in_rescols:
        p_col = p_col + trace.rescols[isotope]['pressure']
    trace.rescols[common]['pressure'] = p_col
    try:
        for isotope in in_rescols:
            r_col = r_col + trace.rescols[isotope]['ratio'] * trace.rescols[isotope]['pressure'] / p_col
        trace.rescols[common]['ratio'] = r_col
    except KeyError:
        pass
        """
