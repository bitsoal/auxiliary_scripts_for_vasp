#!/bin/bash
# Author: Computational Material Physics (Feng & Zhang)
# Website: https://2dmaterials.nus.edu.sg/feng/
# Email: e0001020@u.nus.edu


#Application:
#           This script splits the overall DOSCAR into a set of sub-DOS files for each atom in the crystal
#How to use:
#          . SplitDos.sh


sep=$(cat test_DOSCAR | sed -n '6p')   # the 6th line is the separator line which separates the partial DOS of each atom from each other
sep_lines=$(nl -n "ln" test_DOSCAR | grep "$sep" | cut -d " " -f 1) #get the line indice of the separator
i=0

for ind in $sep_lines #This loop assign the line indice obtained abobe to an array named sep_indice
do
    sep_indice[$i]=$ind
    ((i++))
done
sep_indice[$i]=$(wc -l test_DOSCAR | cut -d " " -f 1) #append the line index of the last line of the DOSCAR to sep_indice
i_max=$(($i-1)) 


for ((i=0; i<i_max; i++))
do
    filename=DOS$(printf "%0${#i_max}d" $i) #format filenames
    cat test_DOSCAR | sed -n "${sep_indice[$i]},${sep_indice[$(($i+1))]}p" >> $filename
done
filename=DOS$(printf "%0${#i_max}d" $i)
cat >$filename <<!
$(cat test_DOSCAR | sed -n "${sep_indice[$i]},${sep_indice[$(($i+1))]}p")
$sep
!


#############################################################################################################################################
##                                                  >>>Process_Partial_Dos.py<<<
##!/usr/bin/env python
## encoding: utf-8
## Author: yang tong
## Email: e0001020@u.nus.edu
#
#
##Application: This python script with another bash script named SplitDos.sh is coded
##             for the post-process of a partial DOS calculation using VASP when the
##             the total (spin-polarized) DOSs of specific atoms in a crystal system
##             is required.
##How to use:
##          1. Use SplitDos.sh to split a DOSCAR into a series of sub-DOS files for each atom with
##             filenames DOS + the atom index in the POSCAR
##          2. In order to implement this code, you need to supply a input filename, say setting_file.
##             In this file, each line is a summation action and the total DOS of the given atoms will be written
##             into a five filename.
##             Parameters required:
##                                - ISPIN=1 | 2: whether it is a spin-polarized calculation: set "ISPIN=2" if it is; otherwise "ISPIN=1".
##                                - ATOMS= indice: List the atom indice as in the POSCAR. The indice are separated by a comma.
##                                - MAX_INDEX=len(atom_index) : The maximum atom index
##                                - FILENAME= e.g setting_file: Give a filename which you want to store the calculated total DOS
##                 ***The three paramters should be separted by whitespaces or tabs
##                 ***There shold not be any whitespace or tabs between parameter names and =, as well as between = and your given values
##
##             command: python + the name of this script + the setting file (e.g. setting_file)
##
##
##          >>>Note that the first and the last lines of the split DOS file for each are the separator which separates the partial DOS of each atom
##              in the overall DOSCAR. Between the two lines are the real partial DOS for each atom.<<<
#
#
#
#import re, sys
#
#class DOS(object):
#
#    """
#    ---The data structure of the input object DOS_data should be a list or tuple,
#            each entry of which is also a list or tuple with the same length.
#
#    ---For a non-spin polarized calculation, the first element of each entry is energy
#        and is followed by the other elements (partial DOS): s, p, d, f DOS, which are
#        in the order of the angular momentum quantum numbers. Elements are separted by a tab.
#
#    ---For a spin polarized calculation, each entry has a form of energy, s_spin_up DOS,
#        s_spin_down DOS, p_spin_up DOS, p_spin_down DOS, d_spin_up DOS, d_spin_down DOS,
#        f_spin_up DOS, f_spin_down DOS... Elements are separted by a tab.
#
#    ---The data sturcutre of the input object DOS_data is intialized as a tuple of tuples
#        such that it cannot be changed.
#    """
#
#    def __init__(self, DOS_data):
#        """
#        convert the data structure of DOS_data to a tuple of tuples
#        and check if the input data is empty.
#        """
#
#        assert isinstance(DOS_data, (tuple, list)), "The input argument should a list or tuple."
#        assert isinstance(DOS_data, (tuple, list)), "Each entry of the input argument should be a list or tuple"
#        if not DOS_data:
#            raise Exception, "the input data is empty."
#        if not DOS_data[0]:
#            raise Exception, "the input data is empty."
#
#        self.DOS = list(DOS_data)
#        self.DOS_len = len(self.DOS)
#        self.entry_len = len(self.DOS[0])
#        for ind, entry in enumerate(self.DOS):
#            self.DOS[ind] = tuple(entry)
#        self.DOS = tuple(self.DOS)
#
#    def combine_DOS(self, other):
#
#        """
#        Combine this DOS with another DOS that must have the same data structure and format.
#        """
#
#        assert isinstance(other, self.__class__), "Cannot combine the two DOS files of different formats"
#        if len(self.DOS) != len(other.DOS) or len(self.DOS[0]) != len(other.DOS[0]):
#            print("Cannot combine <--- The two DOS files have different data structures.")
#        combined_DOS = [ ]
#        for i, j in zip(self.DOS, other.DOS):
#            combined_DOS.append([i[0]] + list([m + n for m, n in zip(i[1:], j[1:])]))
#        return DOS(combined_DOS)
#
#    def cal_total_DOS(self):
#
#        """
#        Calculate the total DOS and return as a DOS instance.
#        """
#        total_DOS = [ ]
#        for entry in self.DOS:
#            total_DOS.append([entry[0], sum(entry[1:])])
#        return DOS(total_DOS)
#
#    def cal_spin_polarized_DOS(self):
#
#        """
#        Calculate the total spin up and spin down DOS separately
#        and return them as a tuple of two DOS instances.
#        """
#
#        spin_polarized_DOS = [ ]
#        for entry in self.DOS:
#            spin_polarized_DOS.append([entry[0], sum(entry[1::2]), sum(entry[2::2])])
#        return DOS(spin_polarized_DOS)
#
#    def write(self, filename):
#
#        """
#        Write the DOS to a given filename. Each line is composed of energy followed by:
#        case 1: partial DOS in s, p, d, f, ...
#        case 2: total DOS only if it is an object returned by self.cal_total_DOS()
#        case 3: total spin up DOS and total spin down DOS if it is an object returned by self.cal_spin_polarized_DOS()
#
#        Numbers in each line are separated by a tab.
#        """
#
#        f = open(filename, 'w')
#        for entry in self.DOS:
#            line_format = "%.8f\t" * self.entry_len + '\n'
#            f.write(line_format % tuple(entry))
#        f.close()
#
#    def __call__(self):
#        return self.DOS
#
#    def __add__(self, other):
#        return self.combine_DOS(other)
#
#
#def read(filename):
#    f = open(filename, 'r')
#    data = [ ]
#    for line in f:
#        line = line.strip("\n ")
#        m = re.findall("[\.\+\-E0-9]+", line)
#        if m:
#            m = [float(num) for num in m]
#            data.append(tuple(m))
#    f.close()
#    data.pop(0); data.pop(-1)
#    return tuple(data)
#
#def parse_setting_file(filename):
#    """parse the setting file. Refer to the leading comment for more details"""
#
#    setting = [ ]
#    f = open(filename)
#    for line in f:
#        ispin = re.search("ISPIN=([1-2])", line)
#        atom_indice = re.search("ATOMS=([0-9\,]+)", line)
#        max_index = re.search("MAX_INDEX=(\d+)", line)
#        file_name = re.search("FILENAME=([a-zA-Z0-9\_]+)", line)
#        if not (ispin and atom_indice and file_name and max_index):
#            raise Exception, "Yours setting is not proper. Please check!"
#        else:
#            ispin = ispin.group(1)
#            atom_indice = atom_indice.group(1)
#            max_index= int(max_index.group(1))
#            file_name = file_name.group(1)
#
#        atom_indice = re.findall("\d+", atom_indice)
#        if atom_indice:
#            atom_indice = [int(index) for index in atom_indice]
#            for index in atom_indice:
#                assert index < max_index, "The given atom indice should be less than or equal to the maximum atom index"
#        else:
#            raise Exception, "Yours setting is not proper. Please check!"
#        string_format = "DOS%0" + str(len(str(max_index))) + "d"
#        dos_filenames = [string_format % index for index in atom_indice]
#        setting.append([ispin, dos_filenames, file_name])
#    return setting
#
#def test():
#    """This is a testing"""
#
#    dos1 = DOS(read("DOS123"))
#    dos2 = DOS(read("DOS052"))
#
#    print("total_dos_cal_test ---> total_dos_cal_test")
#    dos1.cal_total_DOS().write("total_dos_cal_test")
#
#    print("spoin_polarized_DOS_cal_test --> cal_spin_polarized_dos_test")
#    dos1.cal_spin_polarized_DOS().write("cal_spin_polarized_dos_test")
#
#    print("dos_combination_test --> dos_combination_test")
#    (dos1 + dos2).write("dos_combination_test_test")
#
#
#if __name__ == '__main__':
#    for setting in parse_setting_file(sys.argv[1]):
#        dos = DOS(read(setting[1][0]))
#        for file_name in setting[1][1:]:
#            dos = dos + DOS(read(file_name))
#        if setting[0] == '2':
#            dos = dos.cal_spin_polarized_DOS()
#        if setting[0] == '1':
#            dos = dos.cal_total_DOS()
#        dos.write(setting[-1])
#

