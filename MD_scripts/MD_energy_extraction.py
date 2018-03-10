#!/usr/bin/env python

#how to use: python MD_energy_extraction.py OSZICAR1 [OSZICAR2] [OSZICAR3] ...

import re
import sys

def extract_energy_from_OSZICAR(step_index, average_energy, output_file, input_filename):
    f = open(input_filename, "r")
    for line_no, line in enumerate(f):
        if "T" in line and "E" in line and "F" in line and "E0" in line:
            m = re.search("E0=\s([-+E0-9\.]+)", line)
            assert m, line
            average_energy = (average_energy * (step_index-1) + float(m.group(1)))/step_index
            output_file.write("%5d\t%.8f\t%.8f\t%r\t%d\n" % (step_index, float(m.group(1)), average_energy, input_filename, line_no))
            step_index += 1
    f.close()

    return output_file, step_index, average_energy
    
    output.close()

if __name__ == "__main__":

    OSZICARs = sys.argv[1:]
    print(OSZICARs)
    #print(sys.argv)
    output_filename = raw_input("Enter the filename into which the extracted energy will be stored: ")
    output_file = open(output_filename, "w")
    output_file.write("step_index\tenrgy(E0)\t<E0>\tinputfile\tline no in inputfile\n")

    step_index = 1
    average_energy = 0
    for OSZICAR in OSZICARs:
#        print("output_file:", type(output_file))
        #print("atom1, atom2", atom1, atom2)
        #print("bond index", bond_index)
        print("input_filename:", OSZICAR)
        output_file, step_index, average_energy = extract_energy_from_OSZICAR(step_index, average_energy, output_file, input_filename=OSZICAR)

    output_file.close()
