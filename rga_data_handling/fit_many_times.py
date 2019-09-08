# Collection of functions to fit one RGA timetrace many times, with randomized cracking patterns
# Thereby obtaining average values and errobars.
# Also added other functions...

# Aleksander Drenik
# IPP Garching
# October 2016
# aleksander.drenik@ipp.mpg.de

from .Class import perturb_CP
from copy import deepcopy
import matplotlib.pyplot as plt
import scipy as sp
from os import path
import glob
import pickle

# default size definition for plotting of results

default_sz = {'legend': 10, 'label_x': 14, 'label_y': 14, 'tick_x': 14, 'tick_y': 14}

def fit_many_times(container,times,perturb,to_join,*args,**kwargs):
    if args == None:
        return []
    traces = []
    tot_dur = 0
    num_of_goes = times
    # store the original containers cracking patterns
    CP_orig = deepcopy(container.molecules.calib)
    container_orig = deepcopy(container)
    try:
        common_list = to_join[1]
        isotopes_list = to_join[0]
        if len(common_list) == len(isotopes_list):
            make_total = True
        else:
            make_total = False
    except:
        make_total = False
    print("Doing %s for %s :" %(num_of_goes, container.tag['title']))
    for i in range(num_of_goes):
        print(" | %s" %(i + 1)),
        # make a new copy of the original CP,
        # otherwise perturb_CP will override the original copy as well
        CP_new = deepcopy(CP_orig)
        container_new = deepcopy(container_orig)
        container_new.echo = False
        container_new.replace_CP(CP_new)
        perturb_CP(container_new, perturb)
        container_new.deconvolute(*args, **kwargs)
        if make_total:
            for j in range(len(common_list)):
                container_new.join_isotopes(isotopes_list[j], common_list[j])
        traces.append(container_new)
        tot_dur += container_new.glob_duration
    print("deconvolution done, total duration %s seconds" %tot_dur)
    return traces
    
def povpreci(sims):
    
    #HM = sims[0].H_species
    #NHM = sims[0].non_H_species
    
    HM = {}
    NHM = {}
    for key in sims[0].rescols.keys():
        rc_i = sims[0].rescols[key]
        try:
            if 'pressure' in rc_i.keys():
                if 'ratio' in rc_i.keys():
                    HM[key] = ''
                else:
                    NHM[key] = ''
        except:
            pass
    title = sims[0].tag['title']
    tc_name = sims[0].time_col
    tc = sims[0].rescols[tc_name]
    povp_res = {'time': tc}
    gas_list = HM.keys() + NHM.keys()
    for gas in gas_list:
        qt_list = ['pressure']
        if gas in HM.keys():
            qt_list.append('ratio')
        
        povp_res[gas] = {}
        for qt in qt_list:
            povp_res[gas][qt] ={'val': [], 'std': []}
        for ti in range(len(tc)):
            temp = {}
            for qt in qt_list:
                temp[qt] = []
            for simtrace in sims:
                for what in qt_list:
                    temp[what].append(simtrace.rescols[gas][what][ti])
            for what in qt_list:
                povp_res[gas][what]['val'].append(sp.mean(temp[what]))
                povp_res[gas][what]['std'].append(sp.std(temp[what]))
    for gas in gas_list:
        for what in qt_list:
            for how in ['val', 'std']:
                povp_res[gas][what][how] = sp.array(povp_res[gas][what][how])
                
    # same for residual and duration
    for what in ['residual', 'duration']:
        povp_res[what] = {'val': [], 'std': []}
        for ti in range(len(tc)):
            temp = []
            for simtrace in sims:
                temp.append(simtrace.rescols[what][ti])
            povp_res[what]['val'].append(sp.mean(temp))
            povp_res[what]['std'].append(sp.std(temp))
        for how in ['val', 'std']:
            povp_res[what][how] = sp.array(povp_res[what][how])
    povp_res['defs'] = {'HM': HM, 'NHM': NHM, 'title': title, 'time_col': tc_name}
    return povp_res
    
def errorbar_results(povprecni, sz=default_sz, lp=2, residual=True, gl=None, ratio_limit=None):
    
    # TO DO - kaj s casovno skalo...
    
    slika = plt.figure()
    pres = slika.add_subplot(211)
    rat = slika.add_subplot(212)
    tc = povprecni['time']
    HM = povprecni['defs']['HM']
    NHM = povprecni['defs']['NHM']
    title = povprecni['defs']['title']
    
    if gl != None:
        to_plot_H = list(set(HM.keys()) & set(gl))
        to_plot_NH = list(set(NHM.keys()) & set(gl))
    else:
        to_plot_H = HM.keys()
        to_plot_NH = NHM.keys()

            
    for gas in to_plot_H:
        pres.errorbar(tc, povprecni[gas]['pressure']['val'], yerr = povprecni[gas]['pressure']['std'],
                       marker = 'o', label = gas)
        rat.errorbar(tc, povprecni[gas]['ratio']['val'], yerr = povprecni[gas]['ratio']['std'],
                       marker = 'o', label = gas)
    for gas in to_plot_NH:
        pres.errorbar(tc, povprecni[gas]['pressure']['val'], yerr = povprecni[gas]['pressure']['std'],
                       marker = 'o', label = gas)
        
    if residual:
        pres.errorbar(tc, povprecni['residual']['val'], yerr=povprecni['residual']['std'],
                        marker = 'o', label = 'residual')
    pres.set_ylabel("Pressure", fontsize = sz['label_y'])
    rat.set_ylabel("H/(H+D)", fontsize = sz['label_y'])
    
    for plot in [pres, rat]:
        #plot.set_xlim(plot_range)
        plot.legend(loc=lp, fontsize = sz['legend'])
        #plot.yaxis.get_major_formatter().set_powerlimits((0, 1))
        plot.tick_params(axis = 'x', labelsize = 0)
        plot.tick_params(axis = 'y', labelsize = sz['tick_y'])
    rat.set_xlabel('Time [s]', fontsize = sz['label_x'])
    if ratio_limit != None:
        rat.set_ylim(ratio_limit)
    #rat.set_ylim([0, 1.05])
    rat.tick_params(axis = 'x', labelsize = sz['tick_x'])
    slika.suptitle(title, fontsize = sz['label_y'])
    slika.tight_layout(h_pad = -0.25)
    
def one_trace(traces):
    povp = povpreci(traces)
    trc = traces[0]
    rescols = trc.rescols
    rescols['time'] = povp['time']
    
    for gas in trc.H_species.keys():
        rescols[gas] = {'pressure': povp[gas]['pressure']['val'], 'ratio': povp[gas]['ratio']['val']}
    for gas in trc.non_H_species.keys():
        rescols[gas] = {'pressure': povp[gas]['pressure']['val']}
        
    trc.rescols = rescols
    return trc

        
def save_traces(data_folder, traces):
    title = traces[0].tag['title']
    existing = glob.glob(path.join(data_folder,"%s*.pkl" %title))
    if len(existing) == 0:
        t_num = '01'
    else:
        num_list = []
        for file_name in existing:
            basename = path.basename(file_name)
            to_num = basename.replace('%s_' %title, '')
            to_num = to_num.replace('.pkl', '')
            try:
                num = int(to_num)
            except:
                num = 1
            num_list.append(num)
        num = max(num_list)
        t_num = str(num + 1).zfill(2)
                    
    file_name = path.join(data_folder, "%s_%s.pkl" %(title, t_num))
    try:
        p_file = open(file_name, 'wb')
    except:
        print("Error opening %s" %file_name)
        return [10, t_num]
    try:
        pickle.dump(traces, p_file)
    except:
        print("Error saving to %s" %file_name)
        return [11, t_num]
    p_file.close()
    print("Saved traces as %s" %file_name)
    return [0, t_num]
    
def retrieve_traces(data_folder, title, device='RGA3', num=None):
    existing = glob.glob(path.join(data_folder,"%s*.pkl" %title))
    if len(existing) == 0:
        print("No files found for %s" %title)
        return []

    num_list = []
    for file_name in existing:
        basename = path.basename(file_name)
        to_num = basename.replace('%s_' %title, '')
        to_num = to_num.replace('.pkl', '')
        try:
            to_num = int(to_num)
        except:
            to_num = 1
        num_list.append(to_num)
    if num == None:
        t_num = max(num_list)
    else:
        try:
            s_num = int(num)
            if s_num in num_list:
                t_num = s_num
            else:
                print('Number %s not found.' %num),
                t_num = max(num_list)
        except ValueError:
            t_num = max(num_list)
            print('Invalid number specification: %s' %num),
    file_name = path.join(data_folder, "%s_%s.pkl" %(title, str(t_num).zfill(2)))
    if file_name not in existing:
        print("Caramba! %s" %file_name)
        return []
    else:
        print("Loading %s" %file_name)
        p_file = open(file_name, 'rb')
        trs = pickle.load(p_file)
        
        p_file.close()
        return trs
            
        
