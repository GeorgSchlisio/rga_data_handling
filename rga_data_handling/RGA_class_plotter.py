# coding: utf-8
# Prikazovalnik meritev in naracunanih rezultatov
# Aleksander Drenik, IJS
# Oktober 2015

import scipy as sp
import matplotlib.pyplot as plt
from copy import copy
version = '1.2'

sz = {'legend': 12, 'legend_mass': 10, 'napis': 12,'label_x': 16, 'label_y': 16, 'tick_x': 14, 'tick_y': 14}

# Risanje posnetih meritev

def plot_data_profile(profile):
    plot_mode = "plot"
    if len(profile.mass_col) == len(profile.header_int):
        plot_mode = "bar"
    mass_col, intensity_col = profile.mass_col, profile.intensity_col
    try:
        title = profile.tag['title']
    except KeyError:
        title = "Profile plot"
        
  
    slika = plt.figure()
    specplot = slika.add_subplot(1,1,1)
    if plot_mode == "bar":
        plt.bar(mass_col, intensity_col)
    elif plot_mode == "plot":
        plt.plot(mass_col, intensity_col)
    specplot.set_xlabel("Mass [AMU]")
    specplot.set_ylabel("Intensity [arb.u.]")
    specplot.set_title(title)
    specplot.ticklabel_format(useOffset=False)
    #plt.show()
    
def plot_data_trace(trace, masslist=None):
    
    try:
        title = trace.tag['title']
    except KeyError:
        title = "Timetrace plot"
        
    # x axis title
    x_label = trace.time_col_name
    #if 'time_unit' in trace.tag.keys():
    #    if len(trace.tag['time_unit']) > 0:
    #        x_label = "%s [%s]" %(x_label, trace.tag['time_unit'])
    
    slika = plt.figure()
    traceplot = slika.add_subplot(1,1,1)
    if masslist == None:
        list_to_plot = trace.header_int
    else:
        list_to_plot = list (set(masslist) & set(trace.header_int))
    for mass in list_to_plot:
        traceplot.plot(trace.columns[trace.time_col], trace.columns[mass], label = "%s AMU" %mass)
    traceplot.legend(loc=0)
    traceplot.yaxis.get_major_formatter().set_powerlimits((0, 1))
    traceplot.set_xlabel(x_label)
    traceplot.set_ylabel("Intensity [arb.u.]")
    traceplot.set_title(title)
    #plt.show()

# Risanje rezultatov

def show_results_profile(profile):

    masses_of_interest = profile.masses_of_interest
    plot_mode = "plot"
    if len(profile.mass_col) == len(profile.header_int):
        plot_mode = "bar"
    mass_col, intensity_col = profile.mass_col, profile.intensity_col
    title = profile.tag['title']
    resline = profile.resline
    sim_MID_col = profile.sim_MID_col
    H_species = profile.H_species
    non_H_species = profile.non_H_species
    mass_frame = [masses_of_interest[0] - 1, masses_of_interest[-1] + 1]
    mass_area = (mass_col > mass_frame[0]) * (mass_col < mass_frame[-1])
    max_int = max(max(intensity_col[mass_area]), max(sim_MID_col))
    int_frame = [0, 1.1 * max_int]
    
    fig = plt.figure()
    # Plot of recorded spectrum and simulated intensities
    mass_plot = fig.add_subplot(3,1,1)
    #plt.bar(sp.array(masses_of_interest) - 0.4, [sim_MID_col[x] for x in masses_of_interest], color = "red", label = "simulated")
    mass_plot.bar(sp.array(masses_of_interest) - 0.4, sim_MID_col, color = "red", label = "simulated")
    if plot_mode == "bar":
        mass_plot.bar(mass_col, intensity_col, label = "recording", color = "blue")
    elif plot_mode == "plot":    
        mass_plot.plot(mass_col, intensity_col, label = "recording")
    mass_plot.legend(loc=2)
    mass_plot.set_xlim(mass_frame)
    mass_plot.set_ylim(int_frame)
    #plt.plot(masses_of_interest, [sim_MID_col[x] for x in masses_of_interest], label = "simulated", linewidth = 0, marker = 'o')

    mass_plot.set_ylabel("Intensity [arb.u.]")
    mass_plot.set_xlabel("Mass [AMU]")

    # Plot of fitted partial pressures
    
    pressures = {}
    ratios = {}
    for key in resline.keys():
        try:
            pressures[key] = resline[key]['pressure']
        except:
            pass
        try:
            ratios[key] = resline[key]['ratio']
        except:
            pass
    
    pressure_scale = sp.array(range(len(pressures.keys())))
    ratio_scale = sp.array(range(len(ratios.keys())))
    
    pres_plot = fig.add_subplot(3,1,2)
    pres_plot.bar(pressure_scale, pressures.values(), color ='r')
    pres_plot.set_xticks(pressure_scale + 0.5)
    pres_plot.set_xticklabels(pressures.keys())
    pres_plot.set_ylabel('Partial pressure [arb.u.]')
    
    rat_plot = fig.add_subplot(3,1,3)
    rat_plot.bar(ratio_scale, ratios.values(), color = 'r')
    rat_plot.set_xticks(ratio_scale + 0.5)
    rat_plot.set_xticklabels(ratios.keys())
    
    #min_val = 2.2E-10
    #for specimen in H_species:
    #    pressure_bar.append(resline[specimen]["pressure"] + min_val)
    #    ratio_bar.append(resline[specimen]["ratio"] + min_val)
    #    pressure_scale.append(specimen)
    #    ratio_scale.append(specimen)
    #for specimen in non_H_species:
    #    pressure_bar.append(resline[specimen]["pressure"] + min_val)
    #    pressure_scale.append(specimen)
    #pressure_scale.append('residual')
    #pressure_bar.append(resline['residual'])
    #pres_x = sp.arange(len(pressure_bar))
    #scale2 = sp.transpose(sp.array([sp.array(pres_x), sp.array(pressure_scale)]))

    #width = 0.35       # the width of the bars
    #pres_plot = fig.add_subplot(3,1,2)
    #pres_plot.bar(pres_x, pressure_bar, color='r')
    #pres_plot.set_xticklabels(pressure_scale)
    #pres_plot.set_ylabel('Partial pressure [arb.u.]')

    # Plot of H/(H+D) ratios
    #rat_plot = fig.add_subplot(3,1,3)
    #ratio_x = sp.arange(len(ratio_bar))
    #rat_plot.bar(ratio_x, ratio_bar, color='r')
    #rat_plot.set_xticklabels(ratio_scale)
    rat_plot.set_ylabel("H/(H+D)")
    fig.suptitle(title, fontsize = 16)
    for plot in [mass_plot, pres_plot]:
        plot.yaxis.get_major_formatter().set_powerlimits((0, 1))
    
    #fig.tight_layout()
    #plt.show()
    return fig


def show_results_trace(in_trace, gas_list=None):
    trace = copy(in_trace)
    masses_of_interest = trace.masses_of_interest
    columns = trace.columns
    rescols = trace.rescols
    simtracecol = trace.simtracecol
    H_species = trace.H_species
    non_H_species = trace.non_H_species
    time_col = trace.time_col
    if gas_list != None:
        to_plot_H = list(set(H_species.keys()) & set(gas_list))
        to_plot_NH = list(set(non_H_species.keys()) & set(gas_list))
        new_H = {}
        new_NH = {}
        for key in to_plot_H:
            new_H[key] = trace.H_species[key]
        for key in to_plot_NH:
            new_NH[key] = trace.non_H_species[key]
        trace.H_species = new_H
        trace.non_H_species = new_NH
        trace.make_candidates()
        candidates = trace.candidates_dict['candidates']
        masses_of_interest = []
        for mass in range(1,trace.header_int[-1] + 1):
            if sum(candidates)[mass] != 0:
                masses_of_interest.append(mass) 
        
        
    else:
        to_plot_H = H_species.keys()
        to_plot_NH = non_H_species.keys()
    
    try:
        title = trace.tag['title']
    except KeyError:
        title = "Timetrace plot"
    
    # x axis title
    x_label = trace.time_col_name
    #if 'time_unit' in trace.tag.keys():
    #    if len(trace.tag['time_unit']) > 0:
    #        x_label = "%s [%s]" %(x_label, trace.tag['time_unit'])
    
    fig = plt.figure()
    time_step = sp.mean(sp.diff(columns[time_col]))
    time_frame = [rescols[time_col][0] - time_step/2, rescols[time_col][-1] + time_step/2]
    
    # Plot of recorded and simulated masses
    mass_plot = fig.add_subplot(3,1,1)
    for mass in masses_of_interest:
        try:
            mass_plot.plot(columns[time_col],columns[mass],label = "%s AMU" %mass)
        except KeyError:
            mass_plot.plot(columns[time_col],sp.zeros(len(columns[time_col])),label = "%s AMU" %mass)
    plt.gca().set_color_cycle(None)
    for mass in masses_of_interest:
        mass_plot.plot(rescols[time_col],simtracecol[mass],marker = "x", ls = '--')
    mass_plot.plot(rescols[time_col],rescols['residual'],label = 'residual', marker = 'o')
    mass_plot.set_ylabel("Intensity [arb.u.]")
    mass_plot.set_title("RGA intensities")
    
    # Plot of fitted partial pressures
    pres_plot = fig.add_subplot(3,1,2)
    for key in to_plot_H:
        pres_plot.plot(rescols[time_col],rescols[key]['pressure'],marker = 'o', label = key)
    for key in to_plot_NH:
        pres_plot.plot(rescols[time_col],rescols[key]['pressure'],marker = 'o', label = key)
    pres_plot.set_ylabel('p.p. [arb.u.]')
    pres_plot.set_title("Partial pressures")

    # Plot of H/(H+D) ratios
    rat_plot = fig.add_subplot(3,1,3)
    for key in to_plot_H:
        rat_plot.plot(rescols[time_col],rescols[key]['ratio'],marker = 'o', label = key)
    rat_plot.set_ylabel('H/(H+D)')
    rat_plot.set_title("H/(H+D) ratios")
    
    mass_plot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 10, labelspacing = 0)
    pres_plot.legend(bbox_to_anchor=(1.05, 0.7), loc=2, borderaxespad=0., fontsize = 12, labelspacing = 0)
    rat_plot.legend(bbox_to_anchor=(1.05, 0.9), loc=2, borderaxespad=0., fontsize = 12, labelspacing = 0)
    for plot in [mass_plot, pres_plot, rat_plot]:
        plot.set_xlim(time_frame)
        plot.set_xlabel(x_label)
        plot.ticklabel_format(useOffset=False)
    for plot in [mass_plot, pres_plot]:
        plot.yaxis.get_major_formatter().set_powerlimits((0, 1))

    fig.suptitle(title, fontsize = 16)
    fig.tight_layout()
    fig.subplots_adjust(right = 0.8, top = 0.88)
    #plt.show()
    return fig
    
    
def show_calibration_trace(in_trace):
    trace = copy(in_trace)
    masses_of_interest = trace.calib_masses_of_interest
    columns = trace.columns
    calibcols = trace.calibcols
    simtracecol = trace.simtracecol
    
    time_col = trace.time_col
    
    
    try:
        title = trace.calib_tag['title']
    except KeyError:
        title = "Timetrace plot"
    
    # x axis title
    x_label = trace.time_col_name

    
    fig = plt.figure()
    time_step = sp.mean(sp.diff(columns[time_col]))
    time_frame = [calibcols[time_col][0] - time_step/2, calibcols[time_col][-1] + time_step/2]
    
    # Plot of recorded and simulated masses
    mass_plot = fig.add_subplot(3,1,1)
    for mass in masses_of_interest:
        try:
            mass_plot.plot(columns[time_col],columns[mass],label = "%s AMU" %mass)
        except KeyError:
            mass_plot.plot(columns[time_col],sp.zeros(len(columns['time'])),label = "%s AMU" %mass)
    plt.gca().set_color_cycle(None)
    for mass in masses_of_interest:
        mass_plot.plot(calibcols[time_col],simtracecol[mass],marker = "x", ls = '--')
    mass_plot.plot(calibcols[time_col],calibcols['residual'],label = 'residual', marker = 'o')
    mass_plot.set_ylabel("Intensity [arb.u.]", fontsize = sz['label_y'])
    #mass_plot.set_title("RGA intensities")
    
    # Plot of fitted partial pressure and HD ratio
    pres_plot = fig.add_subplot(3,1,2)
    HD_plot = pres_plot.twinx()
    pres_plot.plot(calibcols[time_col], calibcols['pressure'], marker = 'o', label = 'pressure', color = 'blue')
    HD_plot.plot(calibcols[time_col], calibcols['ratio'], marker = 'o', label = 'H/(H+D)', color = 'red')
    pres_plot.set_ylabel('pressure [arb.u.]', fontsize = sz['label_y'], color = 'blue')
    #pres_plot.set_title("Partial pressures")
    HD_plot.set_ylabel("H/(H+D)", fontsize = sz['label_y'], color = 'red')
    
    h1, l1 = pres_plot.get_legend_handles_labels()
    h2, l2 = HD_plot.get_legend_handles_labels()
    pres_plot.legend(h1+h2, l1+l2, bbox_to_anchor=(1.05, 0.7), loc=2, borderaxespad=0., fontsize = sz['legend'], labelspacing = 0)    
    HD_plot.tick_params(axis='y', colors='red')
    HD_plot.spines['right'].set_color('red')
    pres_plot.tick_params(axis='y', colors='blue')
    HD_plot.spines['left'].set_color('blue')

    # Plot of cracking pattern peaks
    peak_plot = fig.add_subplot(3,1,3)
    for peak in sorted(calibcols['peaks'].keys()):
        peak_plot.plot(calibcols[time_col],calibcols['peaks'][peak],marker = 'o', label = "Peak #%s" %peak)
    peak_plot.set_ylabel('Relative peak', fontsize = sz['label_y'])
    #rat_plot.set_title("H/(H+D) ratios")
    
    mass_plot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = sz['legend_mass'], labelspacing = 0)
    #pres_plot.legend(bbox_to_anchor=(1.05, 0.7), loc=2, borderaxespad=0., fontsize = sz['legend'], labelspacing = 0)
    peak_plot.legend(bbox_to_anchor=(1.05, 0.9), loc=2, borderaxespad=0., fontsize = sz['legend'], labelspacing = 0)
    for plot in [mass_plot, pres_plot, peak_plot]:
        plot.set_xlim(time_frame)
        plot.tick_params(axis = 'y', labelsize = sz['tick_y'])
        plot.tick_params(axis = 'x', labelsize = 0)
        plot.ticklabel_format(useOffset=False)
    for plot in [mass_plot, pres_plot]:
        plot.yaxis.get_major_formatter().set_powerlimits((0, 1))
    peak_plot.tick_params(axis = 'x', labelsize = sz['tick_x'])
    peak_plot.set_xlabel(x_label)
    fig.suptitle(title, fontsize = 16)
    fig.tight_layout()
    fig.subplots_adjust(right = 0.8, top = 0.88)
    #plt.show()
    return fig
    
def show_calibration_profile(profile):
    masses_of_interest = profile.calib_masses_of_interest
    plot_mode = "plot"
    if len(profile.mass_col) == len(profile.header_int):
        plot_mode = "bar"
    mass_col, intensity_col = profile.mass_col, profile.intensity_col
    title = profile.calib_tag['title']
    sim_MID_col = profile.sim_MID_col
    mass_frame = [masses_of_interest[0] - 1, masses_of_interest[-1] + 1]
    mass_area = (mass_col > mass_frame[0]) * (mass_col < mass_frame[-1])
    max_int = max(max(intensity_col[mass_area]), max(sim_MID_col))
    int_frame = [0, 1.1 * max_int]
    
    fig = plt.figure()
    # Plot of recorded spectrum and simulated intensities
    mass_plot = fig.add_subplot(2,1,1)
    #plt.bar(sp.array(masses_of_interest) - 0.4, [sim_MID_col[x] for x in masses_of_interest], color = "red", label = "simulated")
    mass_plot.bar(sp.array(masses_of_interest) - 0.4, sim_MID_col, color = "red", label = "simulated")
    if plot_mode == "bar":
        mass_plot.bar(mass_col, intensity_col, label = "recording", color = "blue")
    elif plot_mode == "plot":    
        mass_plot.plot(mass_col, intensity_col, label = "recording")
    mass_plot.legend(loc=2)
    mass_plot.set_xlim(mass_frame)
    mass_plot.set_ylim(int_frame)
    #plt.plot(masses_of_interest, [sim_MID_col[x] for x in masses_of_interest], label = "simulated", linewidth = 0, marker = 'o')

    mass_plot.set_ylabel("Intensity [arb.u.]")
    mass_plot.set_xlabel("Mass [AMU]")

    # Plot of fitted Cracking Pattern peaks
   
    peaks = profile.calib_line['peaks'].keys()
    peak_vals = [1]
    for peak in peaks:
        peak_vals.append(profile.calib_line['peaks'][peak])
    peaks = sp.array([0] + peaks) + 1

    peak_plot = fig.add_subplot(2,1,2)
    peak_plot.bar(peaks, peak_vals, color ='r')
    peak_plot.set_ylabel('Relative peak intensity')
    
    fig.suptitle(title, fontsize = 16)
    for plot in [mass_plot]:
        plot.yaxis.get_major_formatter().set_powerlimits((0, 1))
    
    #fig.tight_layout()
    #plt.show()
    return fig

# skupna funkcija za profile in trace
def plot_data(container, masslist=None, save_name=None):
    if container.type == "Trace":
        fig = plot_data_trace(container, masslist)
    if container.type == "Profile":
        fig = plot_data_profile(container)
    if save_name != None:
        plt.savefig(save_name, dpi=300)
        plt.close()
    else:
        plt.show()
    
def show_results(container, additional=None, save_name=None):
    if container.type == "Trace":
        fig = show_results_trace(container, additional)
    if container.type == "Profile":
        fig = show_results_profile(container)
    if save_name != None:
        plt.savefig(save_name, dpi=300)
        plt.close()
    else:
        plt.show()

def show_calibration(container, save_name=None):
    if container.type == "Trace":
        fig = show_calibration_trace(container)
    if container.type == "Profile":
        fig = show_calibration_profile(container)
    if save_name != None:
        plt.savefig(save_name, dpi=300)
        plt.close()
    else:
        plt.show()

