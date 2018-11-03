#!/usr/bin/env python
# encoding: utf-8
# Author: Computational Material Physics (Feng & Zhang)
# Website: https://2dmaterials.nus.edu.sg/feng/
# Email: e0001020@u.nus.edu


#Application: This python script with another bash script named SplitDos.sh is coded
#             for the post-process of a partial DOS calculation using VASP when the
#             manipulation (i.e. DOS combination or subtraction, total (spin-polarized)
#             DOS calculation, fermi energy correction) of the (spin-polarized) DOSs
#             of specific atoms in a crystal system is required.


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
        self.fermi_energy_correction = False

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
        subtract other's DOS from the DOS that must have the same data structure and format.
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

    def reset_fermi_energy(self, fermi_energy):
        """reset energy with respect to fermi_energy"""

        if not self.fermi_energy_correction: #check whether the fermi energy has been set to be zero
            self.fermi_energy_correction = True
            fermi_energy = float(fermi_energy)
            self.DOS = list(self.DOS)
            for ind, entry in enumerate(self.DOS):
                entry = list(entry)
                entry[0] -= fermi_energy
                self.DOS[ind] = tuple(entry)
            self.DOS = tuple(self.DOS)

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
        print("***Write DOS into %r" % filename)

    def __call__(self):
        return self.DOS

    def __add__(self, other):
        return self.combine_DOS(other)

    def __sub__(self, other):
        return self.subtract_DOS(other)


def test():
    """This is a testing"""

    import sys

    from Process_Partial_Dos import DOS
    from Read_File import read_DOS_file, parse_setting_file


    dos1 = DOS(read_DOS_file("DOS123"))
    dos2 = DOS(read_DOS_file("DOS052"))

    print("total_dos_cal_test ---> write into total_dos_cal_test")
    dos1.cal_total_DOS().write("total_dos_cal_test")

    print("spoin_polarized_DOS_cal_test --> write into cal_spin_polarized_dos_test")
    dos1.cal_spin_polarized_DOS().write("cal_spin_polarized_dos_test")

    print("dos_combination_test --> write into dos_combination_test")
    (dos1 + dos2).write("dos_combination_test_test")

    print("dos_subtraction_test --> write into dos_subtraction_test")
    (dos1 - dos2).write("dos_subtraction_test")

    print("fermi_energy_correction_test (suppose Ef=100) --> write into fermi_energy_correction_test")
    dos1.reset_fermi_energy("100")
    dos1.write("fermi_energy_correction_test")


if __name__ == '__main__':
    test()
