#!/usr/bin/env python
# encoding: utf-8
# Author: Computational Material Physics (Feng & Zhang)
# Website: https://2dmaterials.nus.edu.sg/feng/
# Email: e0001020@u.nus.edu

import sys

from Process_Partial_Dos import DOS
from Read_File import read_DOS_file, parse_setting_file

def test():
    """This is a testing"""

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
#    test()
    for setting in parse_setting_file(sys.argv[1]):
        dos = DOS(read_DOS_file(setting["dos_filenames"][0]))
        if "c" in setting["operation"]:
            for file_name in setting["dos_filenames"][1:]:
                dos = dos + DOS(read_DOS_file(file_name))
        elif setting["operation"] == "sub":
            for file_name in setting["dos_filenames"][1:]:
                dos = dos - DOS(read_DOS_file(file_name))

        if "t" in setting["sum_dos"]:
            if setting["ispin"] == '2':
                dos = dos.cal_spin_polarized_DOS()
            if setting["ispin"] == '1':
                dos = dos.cal_total_DOS()

        dos.reset_fermi_energy(setting["fermi_energy"])

        print("***Dos manipulations specified by the above setting are done.")
        dos.write(setting["file_name"])

