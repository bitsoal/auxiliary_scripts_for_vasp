#!/usr/bin/env python

#how to use: python MD_bond_length_extraction.py XDATCAR1 [XDATCAR2] [XDATCAR3] ...

import re
import sys

def extract_bond_length(atom1, atom2, bond_index, output_file, input_filename, average_bond_length):
    f = open(input_filename, "r")
    f.next()
    scale = float(f.next().strip())
    basis_x = [float(item) for item in re.findall("[-0-9\.]+", f.next())][0]
    basis_y = [float(item) for item in re.findall("[-0-9\.]+", f.next())][1]
    basis_z = [float(item) for item in re.findall("[-0-9\.]+", f.next())][2]
    print(scale)
    print(basis_x)
    print(basis_y)
    print(basis_z) 

    increase_cursor = False
    for line_no, line in enumerate(f):
        if "Direct configuration" in line:
            cursor = 0
            increase_cursor = True
            continue
        if increase_cursor:
            cursor += 1
            if cursor == atom1:
                x1, y1, z1 = [float(item) for item in re.findall("[-0-9\.]+", line)][:3]
            if cursor == atom2:
                x2, y2, z2 = [float(item) for item in re.findall("[-0-9\.]+", line)][:3]
    
                bond_length = scale * pow(((x1-x2)*basis_x)**2 + ((y1-y2)*basis_y)**2 + ((z1-z2)*basis_z)**2, 0.5)
                average_bond_length = (average_bond_length * (bond_index-1) + bond_length)/bond_index
                output_format = "%5d\t" + "%.5f\t"*8 + "%r\t%d\n"
                output_file.write(output_format % (bond_index, x1, y1, z1, x2, y2, z2, bond_length, average_bond_length, input_filename, line_no))        
                bond_index += 1
    
                increase_cursor = False
    f.close()

    return output_file, bond_index, average_bond_length
    
    output.close()

if __name__ == "__main__":
    atoms = raw_input("Enter the indexes of the two atoms whose bond length will be calculated (separate the indexes by a space): ")
    m = [int(item) for item in re.findall("[-0-9\.]+", atoms) if item]
    if len(m)>=2:
        atom1, atom2 = sorted(m[:2])
    

    XDATCARs = sys.argv[1:]
    print(XDATCARs)
    print(sys.argv)
    output_filename = raw_input("Enter the filename into which the extracted bond length of %d with %d will be stored: " % (atom1, atom2))
    output_file = open(output_filename, "w")
    output_file.write("bond_index\tatom1_x\tatom1_y\tatom1_z\tatom2_x\tatom2_y\tatom2_z\tbond_length\taverage_bond_length\tinputfile\tline no in inputfile\n")

    bond_index = 1
    average_bond_length = 0
    for XDATCAR in XDATCARs:
#        print("output_file:", type(output_file))
        print("atom1, atom2", atom1, atom2)
        print("bond index", bond_index)
        print("input_filename:", XDATCAR)
        output_file, bond_index, average_bond_length = extract_bond_length(atom1, atom2, bond_index, output_file, input_filename=XDATCAR, average_bond_length=average_bond_length)

    output_file.close()
