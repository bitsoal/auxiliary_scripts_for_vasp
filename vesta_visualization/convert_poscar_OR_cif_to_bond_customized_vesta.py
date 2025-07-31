#!/usr/bin/env python
# coding: utf-8

# - author: YANG Tong 
# - email: tong-ap.yang@polyu.edu.hk
# - date of creation: 31 July 2025  
# - dependent package: pymatgen (to read the structure file in the cif or poscar format.)  
# - purpose: This python script read either cif-formated or poscar-formated structure files and convert them into the vesta format with customized interatomic bond lengths. Please refer to the documentation of function "read_customization" for how the to-be-converted cif-formated or poscar-formated structure files and the associated bond length tolerance   
# - how to run this script:   
# -- step 1. set the following parameters in the __name__ block:    
#         ** customization_filename (string): the path of the customization filename that function "read_customization" read the to-be-converted structure files and the associated bond length tolerance. For details, please refer to the documentation of function "read_customization"   
#         * structure_folder (string): the path to the folder under which the to-be-converted structures are stored.  
#         * customized_vesta_file_loc (string): the path to the folder under which the generated vesta files are saved.    
#         * bond_tol_shift (float): a tiny shift to the bond length tolerances prescribed in customization_filename. Just set it to zero if not shift is needed.   
# -- step 2. run "python convert_poscar_OR_cif_to_bond_customized_vesta.py"   

# In[1]:


import os
import json

from pymatgen.core import Structure


# In[2]:


#from Cordero, B. et al. Covalent radii revisited. Dalton Trans. 0, 2832–2838 (2008).
atomic_radius = {'Er': 1.89, 'Os': 1.44, 'Ne': 0.58, 'Pb': 1.46, 'Br': 1.2,  'Rb': 2.2,  'U': 1.96,  'Al': 1.21, 'Eu': 1.98, 
                 'Sr': 1.95, 'Hf': 1.75, 'S': 1.05,  'Pu': 1.87, 'Xe': 1.4,  'Lu': 1.87, 'As': 1.19, 'Ge': 1.2,  'F': 0.57, 
                 'Ga': 1.22, 'Cd': 1.44, 'Sn': 1.39, 'Rn': 1.5,  'Pd': 1.39, 'Si': 1.11, 'Pt': 1.36, 'Na': 1.66, 'Yb': 1.87, 
                 'Ca': 1.76, 'Tl': 1.45, 'Gd': 1.96, 'Zr': 1.75, 'Ar': 1.06, 'Fe': 1.52, 'Kr': 1.16, 'Mn': 1.61, 'In': 1.42, 
                 'Nd': 2.01, 'Th': 2.06, 'I': 1.39,  'Np': 1.9,  'N': 0.71,  'He': 0.28, 'Pa': 2.0,  'Mg': 1.41, 'B': 0.84, 
                 'Dy': 1.92, 'Pm': 1.99, 'Po': 1.4,  'Ho': 1.92, 'Zn': 1.22, 'Cr': 1.39, 'Te': 1.38, 'Tb': 1.94, 'Rh': 1.42, 
                 'C':  0.76, 'Fr': 2.6,  'Y': 1.9,   'H': 0.31,  'Bi': 1.48, 'P': 1.07,  'Ti': 1.6,  'Sb': 1.39, 'Tm': 1.9, 
                 'Hg': 1.32, 'Cm': 1.69, 'Li': 1.28, 'Ag': 1.45, 'Be': 0.96, 'W': 1.62,  'V': 1.53,  'Co': 1.5,  'Ir': 1.41, 
                 'Cu': 1.32, 'Ta': 1.7,  'Am': 1.8,  'K': 2.03,  'Cl': 1.02, 'Ru': 1.46, 'Pr': 2.03, 'Ra': 2.21, 'Sc': 1.7, 
                 'Ac': 2.15, 'La': 2.07, 'Mo': 1.54, 'At': 1.5,  'Se': 1.2,  'Re': 1.51, 'Au': 1.36, 'Ni': 1.24, 'Sm': 1.98, 
                 'Tc': 1.47, 'Ce': 2.04, 'Nb': 1.64, 'Cs': 2.44, 'Ba': 2.15, 'O': 0.66}


# In[3]:


def read_customization(customization_filename):
    """
    read customization_filename which has the following format:
    1. Whatever follows sign "#" will be ignored.
    2. Each line comprises no less than 2 entries which are separated by whitespaces.
        2.1. All entries but the last one correspond to a series of to-be-converted structure filenames;
        2.2. The last entry specifies the customized bond tolerance. If there are more than bond tolerances, 
            join them using sign "," only. Note that no whitespaces are allowed in between bond tolerances.
            For atom i and atom j, the bond length threshold between them is (1 + bond_tol)  * (R_i + R_j), 
            where bond_tol is the specified bond length tolerance. 
            R_i and R_j are the covalent radii of atom i and j, respectively.
    output: A list of tuples:
            1. Each tuple corresponds to one valid line in the customization_filename
            2. For each tuple, all entries except the last one are the parsed to-be-convereted structure filenames
            3. The last entry is a tuple of parsed customized bond tolerances.
    """
    with open(customization_filename, "r") as f:
        lines = list(f)
        
    customizations = []
    for line in lines:
        line = line.split("#")[0].strip()
        if line:
            #print(line)
            temp_list = line.split()
            temp_list[-1] = tuple((tol for tol in temp_list[-1].split(",")))
            for filename in temp_list[:-1]:
                full_filename = os.path.join(structure_folder, filename)
                assert os.path.isfile(full_filename), "{} does not exist".format(full_filename)
            customizations.append(tuple(temp_list))
    return customizations


# In[4]:


def get_vesta_format(full_filename, bond_tol, bond_tol_shift, boundary=[0, 1, 0, 1, 0, 1], is_reversed_sbond_included=True):
    struct = Structure.from_file(full_filename)
    
    vesta_formated_list = ["""#VESTA_FORMAT_VERSION 3.5.4


CRYSTAL

TITLE
{}                         

GROUP
1 1 P 1
SYMOP
 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1
 -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0
TRANM 0
 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1
LTRANSL
 -1
 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
LORIENT
 -1   0   0   0   0
 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000
 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000
LMATRIX
 1.000000  0.000000  0.000000  0.000000
 0.000000  1.000000  0.000000  0.000000
 0.000000  0.000000  1.000000  0.000000
 0.000000  0.000000  0.000000  1.000000
 0.000000  0.000000  0.000000
CELLP
 {:.6f}   {:.6f}  {:.6f}  {:.6f}  {:.6f}  {:.6f}
  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000""".format(struct.composition.formula, 
                                                                            struct.lattice.a,     struct.lattice.b,    struct.lattice.c, 
                                                                            struct.lattice.alpha, struct.lattice.beta, struct.lattice.gamma)]
    vesta_formated_list.append("STRUC")
    species_list = [str(spec) for spec in struct.species]
    for atom_ind, spec in enumerate(species_list):
        if atom_ind == 0:
            spec_ind = 1
        else:
            if species_list[atom_ind-1] != spec:
                spec_ind = 1
        frac_x, frac_y, frac_z = struct.frac_coords[atom_ind]        
        vesta_formated_list.append("  {} {}        {}{}{}  1.0000   {:.6f}   {:.6f}   {:.6f}    1a       1".format(atom_ind+1, spec, 
                                                                                                                   atom_ind+1, spec, spec_ind, 
                                                                                                                   frac_x, frac_y, frac_z))
        vesta_formated_list.append("                            0.000000   0.000000   0.000000  0.00")
        spec_ind += 1
    vesta_formated_list.append("  0 0 0 0 0 0 0")
    
    vesta_formated_list.append("BOUND")
    vesta_formated_list.append("       {}        {}         {}        {}         {}        {}".format(*boundary))
    vesta_formated_list.append("  0   0   0   0  0")
    
    vesta_formated_list.append("SBOND")
    bond_count = 1
    sorted_spec_list = sorted(set(species_list))
    for spec_1_ind, spec_1 in enumerate(sorted_spec_list):
        for spec_2 in sorted_spec_list[spec_1_ind:]:
            customized_bond_tol = (1 + bond_tol + bond_tol_shift) * (atomic_radius[spec_1] + atomic_radius[spec_2])
            vesta_formated_list.append("  {}    {}    {}    0.00000    {:.6f}  0  1  1  0  1  0.250  2.000 127 127 127".format(bond_count, spec_1, spec_2, 
                                                                                                                               customized_bond_tol))
            bond_count += 1
            if is_reversed_sbond_included and spec_1 != spec_2:
                vesta_formated_list.append("  {}    {}    {}    0.00000    {:.6f}  0  1  1  0  1  0.250  2.000 127 127 127".format(bond_count, spec_2, spec_1, 
                                                                                                                                   customized_bond_tol))
                bond_count += 1
    vesta_formated_list.append("  0 0 0 0")
    
    return "\n".join(vesta_formated_list)
    


# In[5]:


if __name__ == "__main__":
    #Refer to the documentation of functional read_customization for the format of customization_filename
    customization_filename="input_data/test_customization.dat"
    structure_folder = "input_data/"
    customized_vesta_file_loc="customized_vesta_files"
    bond_tol_shift = 1.0e-5
    
    if not os.path.isdir(customized_vesta_file_loc):
        os.mkdir(customized_vesta_file_loc)
    else:
        for file in os.listdir(customized_vesta_file_loc):
            os.remove(os.path.join(customized_vesta_file_loc, file))
    
    customization_details = read_customization(customization_filename=customization_filename)
    for cust_detail in customization_details:
        for str_filename in cust_detail[:-1]:
            if str_filename.endswith(".cif"):
                vesta_filename_prefix = str_filename.replace(".cif", "")
            elif str_filename.endswith(".vasp"):
                vesta_filename_prefix = str_filename.replace(".vasp", "")
            else:
                raise Exception("The structure filename must have a file extension of either '.cif' or '.vasp'. But {} does not fulfill this requirement".format(str_filename))
            for bond_tol in cust_detail[-1]:
                full_vesta_filename = os.path.join(customized_vesta_file_loc, "{}_tol_{}.vesta".format(vesta_filename_prefix, bond_tol))
                if not os.path.isfile(full_vesta_filename):
                    formated_vesta_output = get_vesta_format(full_filename=os.path.join(structure_folder, str_filename), 
                                                             bond_tol=float(bond_tol), bond_tol_shift=1e-3, 
                                                             boundary=[0, 1, 0, 1, 0, 1], 
                                                             is_reversed_sbond_included=False)
                    with open(full_vesta_filename, "w") as f:
                        f.write(formated_vesta_output)
                    print("Finished preparing {} from {}".format(full_vesta_filename, str_filename))
                else:
                    print("{} already exists. Skipped it".format(full_vesta_filename))

