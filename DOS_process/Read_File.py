#!/usr/bin/env python
# encoding: utf-8
# Author: Computational Material Physics (Feng & Zhang)
# Website: https://2dmaterials.nus.edu.sg/feng/
# Email: e0001020@u.nus.edu


#Application: This python script can read partial DOS files that derive from SplitDos.sh and format the DOS which
#             can be used the input DOS data of class 'DOS' in Process_Partial_Dos.py.
#             On the other hand, it can extract certain params bellow from a setting, which defines the DOS manipulations.
#             You can look at Process_Partial_Dos.py or Run.py to figure out the function of this script.
#How to use:
#
#          1. Use SplitDos.sh to split a DOSCAR into a series of sub-DOS files for each atom with
#             filenames 'DOS + the atom index in the POSCAR', i.e. DOS001, DOS101, ...
#          2. In order to implement this code, you need to supply a input filename, say setting.
#             In this file, each line is a individual job. The partial DOS combination or the total DOS
#             calculation of the given atoms will be implemented and written into a file as per parameters bellow.
#             Parameters required:
#                                - ISPIN=1 | 2: (optional, default: 1)
#                                             whether it is a spin-polarized calculation: set "ISPIN=2" if it is; otherwise "ISPIN=1".
#                                - ATOMS= indice (i.e. 1,4,102,205): (Must be set)
#                                             List the atom indice as in the POSCAR.
#                                             The indice are ONLY separated by a comma. Whitespaces and tabs are not allowed.
#                                - MAX_INDEX=len(atom_index) : (Must be set)
#                                             The maximum atom index
#                                - OPERATION=COM | SUB: (optional, default: COM)
#                                                              If COM, DOS combination calculation is made ;
#                                                              If SUB, DOS subtraction calculation is made;
#                                - SUM_DOS= TRUE | FALSE: (optional, default: FALSE)
#                                                              If TRUE and ISPIN=2, the total spin polarized DOSs are calculated:
#                                                                      1st col, 2nd col, 3rd col: energy, spin up DOS, spin down DOS
#                                                              If TRUE, the total DOS is calculated: 1st col, 2nd col: energy, total DOS
#                                                              If FALSE, no action
#                                - FERMI_ENERGY: (optional, default: 0)
#                                                              Reset energy with respect to FERMI_ENERGY
#                                - FILENAME= e.g setting_file: (optional, default: output_dos)
#                                             Give a filename which you want to store the calculated DOS
#                                             AFTER ALL DOS MANIPULATIONS. --->This filename will be formated into lower cases.
#                 ***The three paramters should be separted by whitespaces or tabs
#                 ***The parameters could be in either the lower case or the upper case
#
#          >>>Note that the first and the last lines of the split DOS file for each are the separator which separates the partial DOS of each atom
#              in the overall DOSCAR. Between the two lines are the real partial DOS for each atom.<<<



import re, sys, pprint

def read_DOS_file(filename):
    """
    read and format DOS obtained from SplitDos.sh.
    remove the first line and last line which is the operator used by SplitDos.sh
    """

    f = open(filename, 'r')
    data = [ ]
    for line in f:
        line = line.strip("\n ")
        m = re.findall("[\.\+\-E0-9]+", line)
        if m:
            m = [float(num) for num in m]
            data.append(tuple(m))
    f.close()
    data.pop(0); data.pop(-1)
    return tuple(data)

def parse_setting_file(filename):
    """parse the setting file. Refer to the leading comment for more details"""

    setting = [ ]
    f = open(filename)
    for ind, line in enumerate(f):
        ind += 1

        params = {
            "ispin": {
                "reg_exp":re.compile("ispin\s*=\s*([1-2])"),
                "necessity": False,
                "default_value": "1",
                "valid_values": ["1", "2"]
            },
            "atoms": {
                "reg_exp": re.compile(r"atoms\s*=\s*([0-9\,]+)"),
                "necessity": True,
                "valid_values": [ ]
            },
            "max_index": {
                "reg_exp": re.compile("max_index\s*=\s*(\d+)"),
                "necessity": True,
                "valid_values": [ ]
            },
            "file_name": {
                "reg_exp": re.compile("filename\s*=\s*([a-zA-Z0-9\_]+)"),
                "necessity": False,
                "default_value": "output_dos",
                "valid_values": [ ]
            },
            "operation": {
                "reg_exp": re.compile("operation\s*=\s*([a-z]+)"),
                "necessity": False,
                "default_value": "COM",
                "valid_values": ["c", "s", "com", "sub"]
            },
            "sum_dos": {
                "reg_exp": re.compile("sum_dos\s*=\s*([a-z]+)"),
                "necessity": False,
                "default_value": "false",
                "valid_values": ["t", "f", "true", "false"]
            },
            "fermi_energy": {
                "reg_exp": re.compile("fermi_energy\s*=\s*([0-9e\.\+-]+)"),
                "necessity": False,
                "default_value": "0",
                "valid_values": [ ]
            }
        }

        if not line.strip():
            continue
        line = line.lower()
        correct_set = True
        for key, value in params.items():
            m = value["reg_exp"].search(line)
            if m:
                m = m.group(1)
                m = check_param_validation(m, params[key]["valid_values"])
                if m:
                    params[key]["value"] = m
                else:
                    params[key]["value"] = value["default_value"]
            elif value["necessity"] == False:
                params[key]["value"] = value["default_value"]
            else:
                correct_set = False
                print("\n[Error]: <line %d> tag %s must be set! Skip the job\n\tline %d: %r" % (ind, key.upper(), ind, line))
                break
        if not correct_set:
            continue

        m = generate_dos_filename([ind, line], params["atoms"]["value"], params["max_index"]["value"])
        if not m:
            continue
        else:
            params["dos_filenames"] = {"value": m}


        setting = dict([(key, params[key]["value"]) for key in params])
        print("\n[INFO]:<line %d> has been parsed successfully\nline %d: %r" % (ind, ind, line))
        print("Parsed tags:")
        pprint.pprint(setting)
        yield setting
    #return setting

def check_param_validation(parsed_param, valid_values):
    """Check whether parsed_param contains any valid param value lised in valid_values.
    if yes, return the valid param value appearing in parsed_param;
    if no, return False"""

    if not valid_values:
        return parsed_param
    for value in valid_values:
        if value in parsed_param:
            return value
    return False


def generate_dos_filename(line_info, atoms, max_index, *argvs, **kargvs):

    """
    atom_indice is a string containing the atom indice that you want to process.
    This function generate dos filenames in the form of 'DOS + index', i.e. DOS013 if
    max_index = 3
    """

    line_index = line_info[0]
    line = line_info[1]

    atoms = re.findall("\d+", atoms)
    if atoms:
        if len(atoms) <= 1:
            print("\n[Error]: <line %d> The number of atoms specified using ATOMS should be more than one!\n\tline %d: %r" % (line_index, line_index, line))
            return [ ]
        atoms = [int(index) for index in atoms]
        for index in atoms:
            if index > int(max_index):
                print("\n[Error]: <line %d> The atom indice in ATOMS should be less than or equal to MAX_INDEX! Skip the job.\n\tline %d: %r" % (line_index, line_index, line))
                continue
    else:
        print("\n[Error]: <line %d> AMTOMS is not set correctly! Skip this job.\n\tline %d: %r" % (line_index, line_index, line))
        return [ ]
    string_format = "DOS%0" + str(len(max_index)) + "d"
    return [string_format % index for index in atoms]



if __name__ == '__main__':
    for setting in parse_setting_file(sys.argv[1]):
        pass
