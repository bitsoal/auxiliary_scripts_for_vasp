#!/usr/bin/env python

#how to use: python MD_bond_length_extraction.py XDATCAR1 [XDATCAR2] [XDATCAR3] ...

import re
import sys
import math

def number_times_list(number, a_list):
    return [number*ele for ele in a_list]

def elementwise_list_plus_list(list_1, list_2, list_3):
    return [ele_1+ele_2+ele_3 for ele_1, ele_2, ele_3 in zip(list_1, list_2, list_3)]

def map_a_number_into_minus_half_to_half(number):
    """map a number into [-0.5, 0.5]"""
    while number > 0.5:
        number -= 1
    while number < -0.5:
        number += 1
    return number

def cal_angle(vec_1, vec_2):
    inner_product = sum([i*j for i, j in zip(vec_1[:3], vec_2[:3])])
    norm_1 = pow(sum([i*i for i in vec_1[:3]]), 0.5)
    norm_2 = pow(sum([i*i for i in vec_2[:3]]), 0.5)
    return math.acos(inner_product / norm_1 / norm_2)/(2*math.pi)*360

def extract_bond_length(atom1, atom2, bond_index, output_file, input_filename, average_bond_length):
    f = open(input_filename, "r")
    next(f)
    scale = float(next(f).strip())
    basis_a = [float(item) for item in re.findall("[-0-9\.]+", next(f))]
    basis_b = [float(item) for item in re.findall("[-0-9\.]+", next(f))]
    basis_c = [float(item) for item in re.findall("[-0-9\.]+", next(f))]
    print(scale)
    print(basis_a)
    print(basis_b)
    print(basis_c) 

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

                frac_x_diff = map_a_number_into_minus_half_to_half(x1-x2)
                frac_y_diff = map_a_number_into_minus_half_to_half(y1-y2)
                frac_z_diff = map_a_number_into_minus_half_to_half(z1-z2)

                cart_diff_vec = elementwise_list_plus_list(number_times_list(frac_x_diff, basis_a), number_times_list(frac_y_diff, basis_b), number_times_list(frac_z_diff, basis_c))
                bond_length = scale * pow(sum([vec_ele**2 for vec_ele in cart_diff_vec]), 0.5)
                average_bond_length = (average_bond_length * (bond_index-1) + bond_length)/bond_index
                output_format = "%5d\t" + "%.5f\t"*8 + "%r\t%d\n"
                output_file.write(output_format % (bond_index, x1, y1, z1, x2, y2, z2, bond_length, average_bond_length, input_filename, line_no))        
                bond_index += 1
    
                increase_cursor = False
    f.close()

    return output_file, bond_index, average_bond_length


def extract_bond_angle(atom1, atom2, atom3,  bond_index, output_file, input_filename, average_bond_length_1, average_bond_length_2, average_bond_angle):
    f = open(input_filename, "r")
    next(f)
    scale = float(next(f).strip())
    basis_a = [float(item) for item in re.findall("[-0-9\.]+", next(f))]
    basis_b = [float(item) for item in re.findall("[-0-9\.]+", next(f))]
    basis_c = [float(item) for item in re.findall("[-0-9\.]+", next(f))]
    print(scale)
    print(basis_a)
    print(basis_b)
    print(basis_c) 

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
            if cursor == atom3:
                x3, y3, z3 = [float(item) for item in re.findall("[-0-9\.]+", line)][:3]

                frac_x_diff_1 = (x1-x2)
                frac_y_diff_1 = (y1-y2)
                frac_z_diff_1 = (z1-z2)

                frac_x_diff_2 = (x3-x2)
                frac_y_diff_2 = (y3-y2)
                frac_z_diff_2 = (z3-z2)

                cart_diff_vec_1 = elementwise_list_plus_list(number_times_list(frac_x_diff_1, basis_a), number_times_list(frac_y_diff_1, basis_b), number_times_list(frac_z_diff_1, basis_c))
                cart_diff_vec_2 = elementwise_list_plus_list(number_times_list(frac_x_diff_2, basis_a), number_times_list(frac_y_diff_2, basis_b), number_times_list(frac_z_diff_2, basis_c))
                bond_length_1 = scale * pow(sum([vec_ele**2 for vec_ele in cart_diff_vec_1]), 0.5)
                bond_length_2 = scale * pow(sum([vec_ele**2 for vec_ele in cart_diff_vec_2]), 0.5)
                bond_angle = cal_angle(cart_diff_vec_1, cart_diff_vec_2) 
                average_bond_length_1 = (average_bond_length_1 * (bond_index-1) + bond_length_1)/bond_index
                average_bond_length_2 = (average_bond_length_2 * (bond_index-1) + bond_length_2)/bond_index
                average_bond_angle = (average_bond_angle * (bond_index-1) + bond_angle)/bond_index
                output_format = "%5d\t" + "%.5f\t"*15 + "%r\t%d\n"
                output_file.write(output_format % (bond_index, x1, y1, z1, x2, y2, z2, x3, y3, z3, bond_length_1, bond_length_2, bond_angle, average_bond_length_1, average_bond_length_2, average_bond_angle, input_filename, line_no))        
                bond_index += 1
    
                increase_cursor = False
    f.close()

    return output_file, bond_index, average_bond_length_1, average_bond_length_2, average_bond_angle

    

if __name__ == "__main__":
    atoms = input("Enter the 1-based indexes of the two (three) atoms whose bond length (angle: atom_1-atom_2-atom_3) will be calculated (separate the indexes by a space): ")
    m = [int(item) for item in re.findall("[-0-9\.]+", atoms) if item]
    if len(m) == 2:
        atom1, atom2 = m[:2]
    if len(m) >= 3:
        atom1, atom2, atom3 = m[:3]
    

    XDATCARs = sys.argv[1:]
    print(XDATCARs)
    print(sys.argv)
    if len(m) == 2:
        output_filename = input("Enter the filename into which the extracted bond length of %d with %d will be stored: " % (atom1, atom2))
    else:
        output_filename = input("Enter the filename into which the extracted bond angle %d-%d-%d will be stored: " % (atom1, atom2, atom3))
    output_file = open(output_filename, "w")
    if len(m) == 2:
        output_file.write("bond_index\tatom1_x\tatom1_y\tatom1_z\tatom2_x\tatom2_y\tatom2_z\tbond_length\taverage_bond_length\tinputfile\tline no in inputfile\n")
    else:
        output_file.write("bond_index\tatom1_x\tatom1_y\tatom1_z\tatom2_x\tatom2_y\tatom2_z\tatom3_x\tatom3_y\tatom3_z\tbond_length_1\tbond_length_2\tbond_angle\taverage_bond_length_1\taverage_bond_length_2\taverage_bond_angle\tinputfile\tline no in inputfile\n")


    bond_index = 1
    if len(m) == 2:
        average_bond_length = 0
        for XDATCAR in XDATCARs:
    #        print("output_file:", type(output_file))
            print("atom1, atom2", atom1, atom2)
            print("bond index", bond_index)
            print("input_filename:", XDATCAR)
            output_file, bond_index, average_bond_length = extract_bond_length(atom1, atom2, bond_index, output_file, input_filename=XDATCAR, average_bond_length=average_bond_length)
    else:
        for XDATCAR in XDATCARs:
            average_bond_length_1, average_bond_length_2, average_bond_angle = 0, 0, 0
    #        print("output_file:", type(output_file))
            print("atom1, atom2, atom3", atom1, atom2, atom3)
            print("bond index", bond_index)
            print("input_filename:", XDATCAR)
            output_file, bond_index, average_bond_length_1, average_bond_length_2, average_bond_angle = extract_bond_angle(atom1, atom2, atom3, bond_index, output_file, input_filename=XDATCAR, average_bond_length_1=average_bond_length_1, average_bond_length_2=average_bond_length_2, average_bond_angle=average_bond_angle)

    output_file.close()
