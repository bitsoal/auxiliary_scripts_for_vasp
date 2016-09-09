#!/usr/bin/env python
# encoding: utf-8
# Author: Computational Material Physics (Feng & Zhang)
# Website: https://2dmaterials.nus.edu.sg/feng/
# Email: e0001020@u.nus.edu


#Application: This python script with another bash script named SplitDos.sh is coded
#             for the post-process of a partial DOS calculation using VASP when the
#             the total (spin-polarized) DOSs of specific atoms in a crystal system
#             is required.
#How to use:
#             Suppose this script is named 'Process_Partial_Dos.py' and the setting file is named 'setting'
#             command: python Process_Partial_Dos.py setting
#
#          1. Use SplitDos.sh to split a DOSCAR into a series of sub-DOS files for each atom with
#             filenames 'DOS + the atom index in the POSCAR', i.e. DOS001, DOS101, ...
#          2. In order to implement this code, you need to supply a input filename, say setting.
#             In this file, each line is a individual job. The partial DOS combination or the total DOS
#             calculation of the given atoms will be implemented and written into a file as per parameters bellow.
#             Parameters required:
#                                - ISPIN=1 | 2: whether it is a spin-polarized calculation: set "ISPIN=2" if it is; otherwise "ISPIN=1".
#                                - ATOMS= indice (i.e. 1,4,102,205): List the atom indice as in the POSCAR.
#                                         The indice are ONLY separated by a comma. Whitespaces and tabs are not allowed.
#                                - MAX_INDEX=len(atom_index) : The maximum atom index
#                                - FILENAME= e.g setting_file: Give a filename which you want to store the calculated total DOS or
#                                                              the combined partial DOS. --->This filename will be formated into lower cases.
#                                - SUMMATION=TRUE | FALSE: If TRUE, sum partial DOS over orbitals and write the total DOS into FILENAME;
#                                                          If FALSE, no summation
#                 ***The three paramters should be separted by whitespaces or tabs
#                 ***The parameters could be in either the lower case or the upper case
#
#          >>>Note that the first and the last lines of the split DOS file for each are the separator which separates the partial DOS of each atom
#              in the overall DOSCAR. Between the two lines are the real partial DOS for each atom.<<<



import re, sys

class DOS(object):

    """
    ---The data structure of the input object DOS_data should be a list or tuple,
            each entry of which is also a list or tuple with the same length.

    ---For a non-spin polarized calculation, the first element of each entry is energy
        and is followed by the other elements (partial DOS): s, p, d, f DOS, which are
        in the order of the angular momentum quantum numbers. Elements are separted by a tab.

    ---For a spin polarized calculation, each entry has a form of energy, s_spin_up DOS,
        s_spin_down DOS, p_spin_up DOS, p_spin_down DOS, d_spin_up DOS, d_spin_down DOS,
        f_spin_up DOS, f_spin_down DOS... Elements are separted by a tab.

    ---The data sturcutre of the input object DOS_data is intialized as a tuple of tuples
        such that it cannot be changed.
    """

    def __init__(self, DOS_data):
        """
        convert the data structure of DOS_data to a tuple of tuples
        and check if the input data is empty.
        """

        assert isinstance(DOS_data, (tuple, list)), "The input argument should a list or tuple."
        assert isinstance(DOS_data, (tuple, list)), "Each entry of the input argument should be a list or tuple"
        if not DOS_data:
            raise Exception, "the input data is empty."
        if not DOS_data[0]:
            raise Exception, "the input data is empty."

        self.DOS = list(DOS_data)
        self.DOS_len = len(self.DOS)
        self.entry_len = len(self.DOS[0])
        for ind, entry in enumerate(self.DOS):
            self.DOS[ind] = tuple(entry)
        self.DOS = tuple(self.DOS)

    def combine_DOS(self, other):

        """
        Combine this DOS with another DOS that must have the same data structure and format.
        """

        assert isinstance(other, self.__class__), "Cannot combine the two DOS files of different formats"
        if len(self.DOS) != len(other.DOS) or len(self.DOS[0]) != len(other.DOS[0]):
            print("Cannot combine <--- The two DOS files have different data structures.")
        combined_DOS = [ ]
        for i, j in zip(self.DOS, other.DOS):
            combined_DOS.append([i[0]] + list([m + n for m, n in zip(i[1:], j[1:])]))
        return DOS(combined_DOS)

    def subtract_DOS(self, other):

        """
        subtract the DOS from another DOS that must have the same data structure and format.
        """

        assert isinstance(other, self.__class__), "Cannot combine the two DOS files of different formats"
        if len(self.DOS) != len(other.DOS) or len(self.DOS[0]) != len(other.DOS[0]):
            print("Cannot combine <--- The two DOS files have different data structures.")
        subtracted_DOS = [ ]
        for i, j in zip(self.DOS, other.DOS):
            subtracted_DOS.append([i[0]] + list([m - n for m, n in zip(i[1:], j[1:])]))
        return DOS(subtracted_DOS)

    def cal_total_DOS(self):

        """
        Calculate the total DOS and return as a DOS instance.
        """
        total_DOS = [ ]
        for entry in self.DOS:
            total_DOS.append([entry[0], sum(entry[1:])])
        return DOS(total_DOS)

    def cal_spin_polarized_DOS(self):

        """
        Calculate the total spin up and spin down DOS separately
        and return them as a tuple of two DOS instances.
        """

        spin_polarized_DOS = [ ]
        for entry in self.DOS:
            spin_polarized_DOS.append([entry[0], sum(entry[1::2]), sum(entry[2::2])])
        return DOS(spin_polarized_DOS)

    def write(self, filename):

        """
        Write the DOS to a given filename. Each line is composed of energy followed by:
        case 1: partial DOS in s, p, d, f, ...
        case 2: total DOS only if it is an object returned by self.cal_total_DOS()
        case 3: total spin up DOS and total spin down DOS if it is an object returned by self.cal_spin_polarized_DOS()

        Numbers in each line are separated by a tab.
        """

        f = open(filename, 'w')
        for entry in self.DOS:
            line_format = "%.8f\t" * self.entry_len + '\n'
            f.write(line_format % tuple(entry))
        f.close()

    def __call__(self):
        return self.DOS

    def __add__(self, other):
        return self.combine_DOS(other)


def read(filename):
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

        params = {
            "ispin": {
                "reg_exp":re.compile("ispin\s*=\s*([1-2])"),
            },
            "atoms": {
                "reg_exp": re.compile(r"atoms\s*=\s*([0-9\,]+)"),
            },
            "max_index": {
                "reg_exp": re.compile("max_index\s*=\s*(\d+)"),
            },
            "file_name": {
                "reg_exp": re.compile("filename\s*=\s*([a-zA-Z0-9\_]+)"),
            },
            "summation": {
                "reg_exp": re.compile("summation\s*=\s*([a-z]+)"),
            }
        }

        if not line.strip():
            continue
        line = line.lower()
        correct_set = True
        for key, value in params.items():
            m = value["reg_exp"].search(line)
            if m:
                params[key]["value"] = m.group(1)
            else:
                correct_set = False
                print("[Error]: The job setting in line %d is not correct! Skip the job\n\tline %d: %r" % (ind, ind, line))
                break
        if not correct_set:
            continue

        m = generate_dos_filename([ind, line], params["atoms"]["value"], params["max_index"]["value"])
        if not m:
            continue
        else:
            params["dos_filenames"] = {"value": m}


        setting.append(dict([(key, params[key]["value"]) for key in params]))
    return setting

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
            print("[Error]: <line %d> The number of atoms specified using ATOMS should be more than one!\n\tline %d: %r" % (line_index, line_index, line))
            return [ ]
        atoms = [int(index) for index in atoms]
        for index in atoms:
            if index > int(max_index):
                print("[Error]: <line %d> The atom indice in ATOMS should be less than or equal to MAX_INDEX! Skip the job.\n\tline %d: %r" % (line_index, line_index, line))
                continue
    else:
        print("[Error]: <line %d> AMTOMS is not set correctly! Skip this job.\n\tline %d: %r" % (line_index, line_index, line))
        return [ ]
    string_format = "DOS%0" + str(len(max_index)) + "d"
    return [string_format % index for index in atoms]

def no_action(line_info, value, *argvs, **kargvs):
    return value



def test():
    """This is a testing"""

    dos1 = DOS(read("DOS123"))
    dos2 = DOS(read("DOS052"))

    print("total_dos_cal_test ---> total_dos_cal_test")
    dos1.cal_total_DOS().write("total_dos_cal_test")

    print("spoin_polarized_DOS_cal_test --> cal_spin_polarized_dos_test")
    dos1.cal_spin_polarized_DOS().write("cal_spin_polarized_dos_test")

    print("dos_combination_test --> dos_combination_test")
    (dos1 + dos2).write("dos_combination_test_test")


if __name__ == '__main__':
    for setting in parse_setting_file(sys.argv[1]):
        dos = DOS(read(setting["dos_filenames"][0]))
        for file_name in setting["dos_filenames"][1:]:
            dos = dos + DOS(read(file_name))
        if setting["summation"] ==  "true":
            if setting["ispin"] == '2':
                dos = dos.cal_spin_polarized_DOS()
            if setting["ispin"] == '1':
                dos = dos.cal_total_DOS()
        dos.write(setting["file_name"])


########################################################################################################################################
#                                                   >>>SplitDos.sh<<<
##!/bin/bash
#
##Application:
##           This script splits the overall DOSCAR into a set of sub-DOS files for each atom in the crystal
##How to use:
##          . SplitDos.sh
#
#
#sep=$(cat test_DOSCAR | sed -n '6p')   # the 6th line is the separator line which separates the partial DOS of each atom from each other
#sep_lines=$(nl -n "ln" test_DOSCAR | grep "$sep" | cut -d " " -f 1) #get the line indice of the separator
#i=0
#
#for ind in $sep_lines #This loop assign the line indice obtained abobe to an array named sep_indice
#do
#    sep_indice[$i]=$ind
#    ((i++))
#done
#sep_indice[$i]=$(wc -l test_DOSCAR | cut -d " " -f 1) #append the line index of the last line of the DOSCAR to sep_indice
#i_max=$(($i-1))
#
#
#for ((i=0; i<i_max; i++))
#do
#    filename=DOS$(printf "%0${#i_max}d" $i) #format filenames
#    cat test_DOSCAR | sed -n "${sep_indice[$i]},${sep_indice[$(($i+1))]}p" >> $filename
#done
#filename=DOS$(printf "%0${#i_max}d" $i)
#cat >$filename <<!
#$(cat test_DOSCAR | sed -n "${sep_indice[$i]},${sep_indice[$(($i+1))]}p")
#$sep
#!
#
