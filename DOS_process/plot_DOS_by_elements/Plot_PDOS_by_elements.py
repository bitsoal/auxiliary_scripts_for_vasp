
# coding: utf-8

# In[1]:


from pymatgen.io.vasp.outputs import Vasprun
#from pymatgen.electronic_structure.plotter
from pymatgen.io.vasp.outputs import Orbital, Spin
from pymatgen import Structure

import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

import os, glob


# In[2]:


def sum_pdos(pdos_obj, atom_ind_list, orbital_list, spin_list):
    total_pdos = np.zeros(pdos_obj[0][Orbital.s][Spin.up].shape)
    for atom_ind in atom_ind_list:
        for orbital in orbital_list:
            for spin in spin_list:
                total_pdos = total_pdos + pdos_obj[atom_ind][orbital][spin]
    return total_pdos


# In[3]:


def get_atom_ind_of_an_ele(ele_list, target_ele):
    atom_ind_list = []
    for atom_ind, ele in enumerate(ele_list):
        if target_ele == ele:
            atom_ind_list.append(atom_ind)
    return atom_ind_list


# In[27]:


def save_pdos_data(tdos, spin_up_pdos_list, spin_down_pdos_list, projected_ele_list, ISPIN=2, wrt_fermi=True, *args, **argvs):
    if wrt_fermi:
        energies = tdos.energies - tdos.efermi
    else:
        energies = tdos.energies
        
    spin_up_tdos = tdos.densities[Spin.up]
        
    no_of_columns = 1 + len(projected_ele_list)
    no_of_data_points = len(list(spin_up_tdos))
    head_line = "#Energy\tTotal\t" + "\t".join(projected_ele_list) + "\n"
    with open("pdos.dat", "w") as pdos_f:
        pdos_f.write(head_line)
        for data_point_ind in range(no_of_data_points):
            line = "%f\t%f\t" % (energies[data_point_ind], spin_up_tdos[data_point_ind])
            for spin_up_pdos in spin_up_pdos_list:
                line += "%f\t" % spin_up_pdos[data_point_ind]
            line += "\n"
            pdos_f.write(line)
            
        if ISPIN == 2:
            pdos_f.write("\n")
            spin_down_tdos = tdos.densities[Spin.down]
            for data_point_ind in range(no_of_data_points):
                line = "%f\t%f\t" % (energies[data_point_ind], spin_up_tdos[data_point_ind])
                for spin_down_pdos in spin_down_pdos_list:
                    line += "%f\t" % spin_down_pdos[data_point_ind]
                line += "\n"
                pdos_f.write(line)


# In[30]:


def plot_pdos_by_ele(tdos, spin_up_pdos_list, spin_down_pdos_list, projected_ele_list, ISPIN=2, plot_tdos=False, wrt_fermi=True, energy_range_wrt_fermi=3, title=None):
    if wrt_fermi:
        energies = tdos.energies - tdos.efermi
        valid_energy_ind = np.logical_and(energy_range_wrt_fermi[0]+0.1 <= energies, energies<=energy_range_wrt_fermi[1]-0.1)
        energy_range = energy_range_wrt_fermi
    else:
        energies = tdos.energies
        valid_energy_ind = np.logical_and(energy_range_wrt_fermi[0]+tdos.efermi+0.1 <= energies, 
                                          energies<=energy_range_wrt_fermi[1]+tdos.efermi-0.1)
        energy_range = [energy_range_wrt_fermi[0] + tdos.efermi, energy_range_wrt_fermi[1] + tdos.efermi]
     
    spin_up_max_y = max([np.max(pdos[valid_energy_ind]) for pdos in spin_up_pdos_list])
    spin_down_max_y = 0
    if ISPIN == 2:
        spin_down_max_y = max([np.max(np.abs(pdos[valid_energy_ind])) for pdos in spin_down_pdos_list])
    max_y = max([spin_up_max_y, spin_down_max_y])
    
    spin_up_tdos = tdos.densities[Spin.up]
    spin_up_tdos_max = np.max(spin_up_tdos[valid_energy_ind])
    spin_down_tdos = np.zeros(spin_up_tdos.shape)
    spin_down_tdos_max = 0
    if ISPIN==2:
        spin_down_tdos = -1*tdos.densities[Spin.down]
        spin_down_tdos_max = np.max(np.abs(spin_down_tdos[valid_energy_ind]))
        
    if plot_tdos:
        max_y = max([max_y, spin_up_tdos_max, spin_down_tdos_max])
    max_y *= 1.1
    min_y = -max_y if ISPIN == 2 else 0
    
        
    color_list = ["blue", "green", "red", "cyan", "magenta", "yellow", "purple", "orange"]
    
    line_width = 0.5
    plt.cla()
    if plot_tdos:
        plt.plot(energies, spin_up_tdos, label='Total', color="black", linewidth=line_width)
        plt.plot(energies, spin_down_tdos, color="black", linewidth=line_width)
        
    for ele_ind, ele in enumerate(projected_ele_list):
        plt.plot(energies, spin_up_pdos_list[ele_ind], label=ele, color=color_list[ele_ind], linewidth=line_width)
        if ISPIN == 2:
            plt.plot(energies, spin_down_pdos_list[ele_ind], color=color_list[ele_ind], linewidth=line_width)
            
    
    if wrt_fermi:
        plt.plot([0, 0], [min_y, max_y], "--", color="brown", linewidth=line_width)
        plt.xlim(energy_range_wrt_fermi)
    else:
        plt.plot([tdos.efermi, tdos.efermi], [min_y, max_y], "--", color="brown", linewidth=line_width)
        plt.xlim([energy_range_wrt_fermi[0]+tdos.efermi, tdos.efermi+energy_range_wrt_fermi[1]])
    plt.ylim([min_y, max_y])
    plt.legend()
    plt.xlabel("Energy (eV)")
    plt.ylabel("Density of States (states/eV)")
    if title != None:
        plt.title(title)
    plt.tick_params(axis='both', which='both', direction='in')
    plt.savefig("pdos.png", format='png', dpi=500)


# In[31]:


if __name__ == "__main__":
    projected_ele_list = ["Co", "N", "C", "H"]
    plot_total = False
    
    vasprun = Vasprun("vasprun.xml")
    pdos = vasprun.pdos
    tdos = vasprun.tdos
    
    orbital_list = [Orbital.s, Orbital.px, Orbital.py, Orbital.pz, Orbital.dx2, Orbital.dxy, Orbital.dxz, Orbital.dyz, Orbital.dz2]
    spin_up_list, spin_down_list = [], []
    ISPIN = vasprun.incar["ISPIN"]
    for ele in projected_ele_list:
        atom_ind_list = get_atom_ind_of_an_ele(vasprun.atomic_symbols, ele)
        spin_up_list.append(sum_pdos(pdos_obj=pdos, atom_ind_list=atom_ind_list, orbital_list=orbital_list, spin_list=[Spin.up]))
        if ISPIN == 2:
            spin_down_list.append(-1*sum_pdos(pdos_obj=pdos, atom_ind_list=atom_ind_list, orbital_list=orbital_list, spin_list=[Spin.down]))
            
    input_arguments = {"tdos": tdos, 
                       "spin_up_pdos_list": spin_up_list, 
                       "spin_down_pdos_list": spin_down_list, 
                       "projected_ele_list": projected_ele_list, 
                       "ISPIN": 2, 
                       "plot_tdos": True, 
                       "wrt_fermi": True, 
                       "energy_range_wrt_fermi": [-2, 2],
                       "title": "test"}
    save_pdos_data(**input_arguments)
    plot_pdos_by_ele(**input_arguments)
    


