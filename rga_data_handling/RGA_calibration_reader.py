# Reader of RGA calibration files / definitions of cracking patterns
# Handles definitions of both new type (Nov 2018) and old type
# Nov 2018
# Aleksander Drenik
# IPP Garching

from copy import deepcopy

# from .RGA_fitting import export_CP

version = "1.0b"


# expected keys in cracking pattern definitions read out from the new version (Nov 2018) of the file
definition_keys = {
    "H": ["non-H-mass", "H-atoms", "intensity", "peak"],
    "non-H": ["intensity", "peak"],
    "clone": ["source", "non-H-mass", "H-atoms"],
    "general": ["device", "version", "source", "comment", "remark"],
}

# data for the old version file readout
H_mol = ["hydrogen", "ammonia", "water", "methane"]
non_H_mol = ["N2", "O2", "CO2", "CO", "Ar", "Ne", "15N2"]

# Default calibration for old type cracking patterns

default_calib = {
    "version": "default",
    "device": "generic",
    "water": [0.23, 0.011, 0.9],
    "ammonia": [0.8, 0.075, 1.3],
    "methane": [0.858, 0.156, 1.6],
    "hydrogen": [0.05, 0, 1],
    "N2": [0.072, 0.008, 1.0],
    "O2": [0.114, 0.008, 0.96],
    "Ar": [0.107, 0.003, 1.2],
    "CO2": [0.114, 0.085, 1.4],
    "CO": [0.045, 0.009, 1.05],
    "Ne": [0.099, 0.003, 0.23],
}

# dictionary of known non-H containing molecules
# {name: [peak mass 1, peak mass 2, peak mass 3]}

non_H_molecules_d = {
    "N2": [28, 14, 29],
    "15N2": [30, 15, 29],
    "O2": [32, 16, 34],
    "Ar": [40, 20, 36],
    "CO2": [44, 28, 16],
    "CO": [28, 12, 26],
    "Ne": [20, 22, 21],
}

# dictionary of known H-containing molecules:
# {name: name_in_calib, number of H isotopes, total mass of non-H atoms}
H_molecules_d = {}
H_molecules_d["ammonia"] = ["ammonia", 3, 14]
H_molecules_d["ammonia15"] = ["ammonia", 3, 15]
H_molecules_d["water"] = ["water", 2, 16]
H_molecules_d["water18"] = ["water", 2, 18]
H_molecules_d["methane"] = ["methane", 4, 12]
H_molecules_d["methane13"] = ["methane", 4, 13]
H_molecules_d["hydrogen"] = ["hydrogen", 2, 0]


# this function is not in use by any of the modules
# can be used by the end user
def convert_calibration(CP_spec, new_filename=None):
    """Converts an old type definition of cracking patterns to the new type (Nov 2018)"""
    new_filename_used = "default_calibration_new.txt"
    if CP_spec != None:
        if type(CP_spec) == str:
            filename = CP_spec
            if new_filename == None:
                fn, ext = filename.split(".")
                new_filename_used = "%s_new.%s" % (fn, ext)
            else:
                new_filename_used = new_filename

        elif type(CP_spec) == dict:

            if new_filename == None:
                new_filename_used = "unnamed_calibration_new.txt"
            else:
                new_filename_used = new_filename

    calib = read_old_version(CP_spec)
    # write_CP_calibration(CP_spec, new_filename_used)
    calib_as_tab = export_CP(calib)
    write_to_TSV(new_filename_used, calib_as_tab)


def read_old_version(CP_spec):
    """Read an old version of the cracking pattern defition, ouput as new version. If none is provided, default values are returned"""
    calib_new = {"H": {}, "non-H": {}}
    got_calib_file = False
    # If calibration is provided, the default calibration will be overridden.
    # If the calibration is incomplete, data from default calibration will remain.
    calib = deepcopy(default_calib)
    calib["source"] = "default values"
    if CP_spec != None:
        if type(CP_spec) == str:
            filename = CP_spec
            try:
                calib_file = open(filename, "r").read().splitlines()
                got_calib_file = True
                # calib['source'] = filename
            except IOError:
                got_calib_file = False
            # using default values
        elif type(CP_spec) == dict:
            for key in CP_spec:
                calib[key] = CP_spec[key]
            calib["source"] = "direct entry"

    if got_calib_file:
        calib = {}
        calib["source"] = filename
        for line in calib_file:
            if len(line) == 0:
                continue
            if line[0] == "#":
                continue
            elements = line.split("\t")
            if elements[0] in ["version", "device"]:
                calib[elements[0]] = elements[1]
            else:
                try:
                    calib[elements[0]] = map(float, elements[1:])
                except ValueError:
                    continue

    for what in ["device", "version", "source"]:
        try:
            val = calib[what]
        except KeyError:
            val = "unspecified"
        calib_new[what] = val

    calib_new["comment"] = "Converted from %s" % calib["source"]

    for gas in calib:
        if gas in H_mol:
            calib_new["H"][gas] = {"peaks": [1]}
            calib_new["H"][gas]["non-H-mass"] = H_molecules_d[gas][2]
            calib_new["H"][gas]["H-atoms"] = H_molecules_d[gas][1]
            for idx in [0, 1]:
                calib_new["H"][gas]["peaks"].append(calib[gas][idx])
            calib_new["H"][gas]["intensity"] = calib[gas][2]

        if gas in non_H_mol:
            calib_new["non-H"][gas] = {"peaks": [[non_H_molecules_d[gas][0], 1]]}
            for idx in [0, 1]:
                calib_new["non-H"][gas]["peaks"].append(
                    [non_H_molecules_d[gas][idx + 1], calib[gas][idx]]
                )
            calib_new["non-H"][gas]["intensity"] = calib[gas][2]

    return calib_new


def read_calibration_file(filename):
    """Read a new (Nov 2018) version of the calibration file"""
    # errors
    # 1001 - error opening file
    calib = {"H": {}, "non-H": {}, "clone": {}}
    err = 0
    active_name = "dump"
    active_type = "general"
    try:
        with open(filename, "r") as f:
            raw_feed = f.read().splitlines()
    except:
        err = 1001
        return err, calib

    for line in raw_feed:
        if len(line) == 0:
            continue
        while line[0] == " ":
            line = line[1:]
        if line[0] == "#":
            continue

        elements = line.split("\t")
        if len(elements) not in [2, 3]:
            active_type = "general"
            active_name = "dump"
            continue
        if elements[0] in ["H-molecule", "non-H-molecule"]:
            # start reading out the cracking pattern of a H-containing molecule
            active_name = elements[1]
            active_type = elements[0].strip("-molecule")
            calib[active_type][active_name] = {"peaks": [], "intensity": 1}
            continue
        if elements[0] == "H-clone":
            # start reading out a cloned definition
            active_name = elements[1]
            active_type = "clone"
            calib[active_type][active_name] = {}
            continue

        if elements[0] not in definition_keys[active_type]:
            # if elements[0] in definition_keys['general']:
            #    calib[elements[0]] = elements[1]
            active_type = "general"
            active_name = "dump"
            # continue

        if active_type == "general":
            calib[elements[0]] = elements[1]
            continue

        if (elements[0] == "peak") and (active_type in ["H", "non-H"]):
            peak = "unread"
            if active_type == "H":
                peak = float(elements[1])
            elif active_type == "non-H":
                peak = [int(elements[1]), float(elements[2])]
            calib[active_type][active_name]["peaks"].append(peak)
            continue

        try:
            to_val = float(elements[1])
        except:
            # not a number
            calib[active_type][active_name][elements[0]] = elements[1]
            continue

        to_int = int(to_val)
        if to_int == to_val:
            to_val = to_int
        calib[active_type][active_name][elements[0]] = to_val

    # process the H-molecule clone entries
    for clone in calib["clone"]:
        try:
            source = calib["clone"][clone]["source"]
            clone_def = deepcopy(calib["H"][source])

            # replacement entries
            new_keys = set(calib["clone"][clone].keys()) - {"source"}

            for new_key in new_keys:
                clone_def[new_key] = deepcopy(calib["clone"][clone][new_key])
            calib["H"][clone] = clone_def
        except KeyError:
            pass

    return err, calib


# the following functions have been replaced with export_CP from RGA_fitting (already in use by Class)


def write_to_TSV(path, inputlist, line_light="False"):
    """Write a 2D matrix as a TSV file, saved to path"""
    # if line_light:
    #    linebreak = '\n'
    # else:
    #    linebreak = 'test\r\n'
    linebreak = "\r\n"
    outfile = open(path, "wb")
    for line in inputlist:
        outfile.writelines("%s%s" % ("\t".join(map(str, line)), linebreak))
    outfile.close()


'''def write_CP_calibration(calib, filename):
    """Write the cracking pattern definition to a specified filename"""
    # Construction of the output calibration file
    
    out_new = ['# Definition of cracking patterns for rga_data_handling objects']
    out_new.append('# Generated by calibration_convert.py')
    try:
        out_new.append('# Converted from %s' %calib['source'])
    except KeyError:
        pass
    
    for what in ['device', 'version', 'source']:
        try:
            val = calib[what]
        except KeyError:
            val = 'unspecified'
            
        out_new.append('%s\t%s' %(what, val))
        
    for what in ['comment', 'remark']:
        try:
            val = calib[what]
        except KeyError:
            continue
            
        out_new.append('# %s' %(val))

    out_new.append('')
    out_new.append('# Definitions of cracking patterns start here.')
    out_new.append('')

    for gas in calib['H']:
        out_new.append('H-molecule\t%s' %gas)
        CP_def = calib['H'][gas]
        for what in ['non-H-mass', 'H-atoms']:
            out_new.append('%s\t%s' %(what, CP_def[what]))
        for peak in CP_def['peaks']:
            out_new.append('peak\t%s' %peak)
        for what in ['intensity']:
            out_new.append('%s\t%s' %(what, CP_def[what]))
        out_new.append('')
        
    for gas in calib['non-H']:
        out_new.append('non-H-molecule\t%s' %gas)
        CP_def = calib['non-H'][gas]
        for peak in CP_def['peaks']:
            out_new.append('peak\t%s\t%s' %(peak[0], peak[1]))
        for what in ['intensity']:
            out_new.append('%s\t%s' %(what, CP_def[what]))
        out_new.append('')
            
    write_err = write_lines(filename, out_new)
    return write_err
    

def write_lines(path, inputlist, line_light='False'):
    """write a list of strings as lines in a file"""
    # errors:
    # 101 - error opening the file
    err = 0
    try:
        outfile = open(path,'w')
    except IOError:
        err = 101
        return err
    if line_light:
        linebreak = '\n'
    else:
        linebreak = '\r\n'
    for line in inputlist:
        outfile.writelines("%s%s" %(line, linebreak))
    outfile.close()
    return err'''
