
# coding: utf-8

# In[1]:


#Application: This script reads the calculated band structure from vasprun.xml and estimates effective masses of electrons and holes.
#             You need to specify the number of the k-points next to a local extremum that are used for mass estimation. We call
#             it k_window. Note k_window is just the number of the k-points included in either the right or the left half parabola. 
#             The total number of k-points that are used for effective mass estimation is actually 2*k_window + 1
#How to use it:
#             python cal_effective_mass_VASP.py
#             You will be prompted to input k_window
#Output:
#            - for a non-spin-polarized calculation, a band structure figure will be generated, in which the fitted parabolas
#                    and estimated effective masses are shown.
#            - for a spin-polarized calculation, besides the figure of the total band structure, spin-up and spin-down band 
#              structures are plotted separately, where the fitted parabolas and estimated effective masses are also shown.
#Note that the script can only process the band structures from the spin-polarized or non-spin-polarized VASP calculations. 
#     Not sure if it is valid for the band structures based on VASP non-collinear calculations.
#
#Python Package requirements: numpy, matplotlib

#Author: YANG Tong
#Email: e0001020@u.nus.edu, bitsoal@gmail.com


# In[2]:


import os
import numpy as np

import matplotlib.pyplot as plt
#%matplotlib inline


# def read_recp_latt_from_vasprun_xml(cal_loc='.'):
#     rec_basis = []
#     with open(os.path.join(cal_loc, "vasprun.xml"), "r") as f:
#         for line in f:
#             if "<varray name=\"rec_basis\" >" in line:
#                 break
#         rec_basis_lines = [next(f) for i in range(3)]
#     [print(line) for line in rec_basis_lines]
#     for line in rec_basis_lines:
#         line = line.split(">")[1].split("<")[0].strip()
#         rec_basis.append([float(num) for num in line.split()])
#     return np.array(rec_basis)
# read_recp_latt_from_vasprun_xml("data/")

# In[3]:


class VasprunXml():
    
    def __init__(self, cal_loc=".", set_fermi_to_zero=True, k_window=10, skip_no_of_kpoints=0):
        """
        optional arguments:
            cal_loc (string): the absolute path to vasprun.xml from which the band structure is extracted and thereby the effective masses are estimated.
                            default: "."
            set_fermi_to_zero (Boolean): If the Fermi level is set to zero.
                            default: True
            k_window (integer): The number of the k-points next to a local extrema that are used to estimate the effective mass.
                                This is just the number of the k-points in the either left or right half parabola. So the total
                                number of k-points around a local extrema that are used for the effective mass estimation is 2*k_window+1.
                                If a local extrema is located nearby the left (right) boundary of the k-path such that its left
                                (right) side has less than k_window k-points, all k-points on the left (right) side of the local
                                extrema are used for the estimation.
                            default: 10
            skip_no_of_kpoints (integer): For the band structure calculation at the HSE06 level of theory, the k-path is
                                the k-pionts of a scf calculation + the k-points of a routine DFT band structure calculation using
                                a segment-like k-point setting. In this case, the first part of k-points originating from the scf
                                calculation should be abandoned. So we can use this tag to skip those useless k-points.
                            default: 0
        """
        assert os.path.isfile(os.path.join(cal_loc, "vasprun.xml")), "cannot find vasprun.xml under directory %s" % cal_loc
        self.set_fermi_to_zero = set_fermi_to_zero
        self.cal_loc = cal_loc
        self.skip_no_of_kpoints = skip_no_of_kpoints
        self.eigenvalues, self.occupations = self.get_eigenvalues(set_fermi_to_zero=set_fermi_to_zero)
        self.kpoint_list = self.get_kpoint_list()
        self.kpath = self.get_kpath()
        self.k_window = k_window
                                

        
    def get_recp_latt(self):
        """
        parse the reciprocal lattice vectors from vasprun.xml. 
        It should be noted that VASP implicitly factorizes 2*pi out of the reciprocal lattice vectors, that is, (a_i, b_i)=1.
        """
        rec_basis = []
        with open(os.path.join(self.cal_loc, "vasprun.xml"), "r") as f:
            for line in f:
                if "<varray name=\"rec_basis\" >" in line:
                    break
            rec_basis_lines = [next(f) for i in range(3)]
        #[print(line) for line in rec_basis_lines]
        for line in rec_basis_lines:
            line = line.split(">")[1].split("<")[0].strip()
            rec_basis.append([float(num) for num in line.split()])
        return rec_basis
    
    def get_kpoint_list(self):
        kpoint_list = []
        with open(os.path.join(self.cal_loc, "vasprun.xml"), "r") as f:
            for line in f:
                if "<varray name=\"kpointlist\" >" in line:
                    break
            for line in f:
                if "</varray>" not in line:
                    line = line.split(">")[1].split("<")[0].strip()
                    kpoint_list.append([float(num) for num in line.split()])
                else:
                    break
        return kpoint_list[self.skip_no_of_kpoints:]
    
    def get_kpath(self):
        kpath = [0]
        kpoint_list = np.array(self.kpoint_list)
        rec_basis = np.array(self.get_recp_latt())
        for prev_k, cur_k in zip(kpoint_list[0:-1], kpoint_list[1:]):
            d_k_frac_vec = cur_k - prev_k
            d_k_vec = d_k_frac_vec[0] * rec_basis[0] + d_k_frac_vec[1] * rec_basis[1]  + d_k_frac_vec[2] * rec_basis[2]
            d_k = pow(sum([dk**2 for dk in d_k_vec]), 0.5)
            kpath.append(d_k+kpath[-1])
        return kpath
    
    def get_high_symmetry_kpoints_list(self):
        skip_no_of_kpoints = self.skip_no_of_kpoints
        kpoint_list = self.kpoint_list[skip_no_of_kpoints:]
        high_sym_kpoint_list = [kpoint_list[0]]
        for prev_kpoint, cur_kpoint in zip(kpoint_list[1:-2], kpoint_list[2:-1]):
            if prev_kpoint == cur_kpoint:
                high_sym_kpoint_list.append(prev_kpoint)
        high_sym_kpoint_list.append(kpoint_list[-1])
        return high_sym_kpoint_list
    
    def get_high_symmetry_kpoint_index_list(self):
        skip_no_of_kpoints = self.skip_no_of_kpoints
        high_sym_kpoint_index_list = []
        kpoint_list = self.kpoint_list
        high_sym_kpoint_list = self.get_high_symmetry_kpoints_list()
        for kpoint_ind, kpoint in enumerate(kpoint_list):
            if kpoint in high_sym_kpoint_list:
                high_sym_kpoint_index_list.append(kpoint_ind)
        return high_sym_kpoint_index_list
    
    def get_fermi_level(self):
        with open(os.path.join(self.cal_loc, "vasprun.xml"), "r") as f:
            for line in f:
                if "<i name=\"efermi\">" in line:
                    return float(line.split(">")[1].split("<")[0])
                
    def get_eigenvalues(self, set_fermi_to_zero=True):
        if set_fermi_to_zero:
            E_fermi = self.get_fermi_level()
        else:
            E_fermi = 0
            
        eigenvalues = [[], []]
        spin_index = 0
        with open(os.path.join(self.cal_loc, "vasprun.xml"), "r") as f:
            is_eigenvalue_block = False
            eigenvalue_block = []
            for line in f:
                if "<eigenvalues>" in line:
                    is_eigenvalue_block = True
                    eigenvalue_block.append([])
                elif "</eigenvalues>" in line:
                    is_eigenvalue_block = False
                elif is_eigenvalue_block:
                    eigenvalue_block[-1].append(line)
        eigenvalues, occupations = [], []
        for eigen_bloc in eigenvalue_block:
            eigenvalues.append([])
            occupations.append([])
            for line in eigen_bloc:
                if "<set comment=\"spin" in line:
                    eigenvalues[-1].append([])
                    occupations[-1].append([])
                    continue
                if "<set comment=\"kpoint" in line:
                    eigenvalues[-1][-1].append([])
                    occupations[-1][-1].append([])
                if "<r>" in line:
                    eig , occ = [float(num) for num in line.split(">")[1].split("<")[0].strip().split()]
                    eig = eig-E_fermi
                    eigenvalues[-1][-1][-1].append(eig)
                    occupations[-1][-1][-1].append(occ)
        return eigenvalues[0], occupations[0]
    
    def get_band_str_data(self):
        eigenvalues = self.eigenvalues
        energies = [[], []] #spin up & spin down
        no_of_bands = len(eigenvalues[0][0])
        for spin_ind, spin_resolved_eigen in enumerate(eigenvalues):
            for band_ind in range(no_of_bands):
                energies[spin_ind].append([])
                for eigen_list in spin_resolved_eigen:
                    energies[spin_ind][-1].append(eigen_list[band_ind])
        return self.kpath, energies
    
    def save_band_str_data(self, filename="band_str.dat"):
        kpath, energies= self.get_band_str_data()
        with open(filename, "w") as f:
            for spin_resolved_energies in energies:
                for sub_band in spin_resolved_energies:
                    for k_, e_ in zip(kpath, sub_band):
                        f.write("%f    %f\n" % (k_, e_))
                    f.write("\n")
                    
                    
    def get_edge_band(self):
        eigenvalues, occupations = self.eigenvalues, self.occupations
        lowest_conduction_band, highest_valence_band = [], []
        for spin_resolved_occ, spin_resolved_eigen in zip(occupations, eigenvalues):
            lowest_conduction_band.append([])
            highest_valence_band.append([])
            for occ_list, eigen_list in zip(spin_resolved_occ, spin_resolved_eigen):
                #the occupation list and eigenvalue list of a given kpiont
                lcb_band_ind = sum([0 if occ<0.5 else 1 for occ in occ_list])
                lowest_conduction_band[-1].append(eigen_list[lcb_band_ind])
                highest_valence_band[-1].append(eigen_list[lcb_band_ind-1])
        return highest_valence_band, lowest_conduction_band
    
    def find_local_extrema_of_edge_band(self):
        hvb, lcb = self.get_edge_band() # hvb: highest_valence_band; lcb: lowest_conduction_band
        hvb_local_maximums_ind, lcb_local_minimums_ind = [], []
        hvb_varitaion, lcb_variation = [], []
        for spin_resolved_hvb, spin_resolved_lcb in zip(hvb, lcb):
            hvb_variation_list, hvb_local_maximums_ind_list = VasprunXml.find_local_maximums_of_a_band(band=spin_resolved_hvb, k_window=self.k_window)
            lcb_variation_list, lcb_local_minimums_ind_list = VasprunXml.find_local_minimums_of_a_band(band=spin_resolved_lcb, k_window=self.k_window)
            hvb_local_maximums_ind.append(hvb_local_maximums_ind_list)
            lcb_local_minimums_ind.append(lcb_local_minimums_ind_list)
            hvb_varitaion.append(hvb_variation_list)
            lcb_variation.append(lcb_variation_list)
        return hvb_varitaion, hvb_local_maximums_ind, lcb_variation, lcb_local_minimums_ind
    
    @classmethod
    def find_local_maximums_of_a_band(cls, band, k_window):
        band_variation_list = [0]
        for prev_e, cur_e in zip(band[:-1], band[1:]):
            if prev_e > cur_e:
                band_variation_list.append(-1)
            elif prev_e < cur_e:
                band_variation_list.append(1)
            else:
                band_variation_list.append(0)

        #print(band)
        #print("variation", band_variation_list)
        local_maximums_ind_list = []
        right_threshold = len(band)-k_window-1

        
        for band_var_ind, band_var in enumerate(band_variation_list):
            if band_var_ind < k_window:
                left_points = band_variation_list[:band_var_ind]
                right_points = band_variation_list[band_var_ind+1:band_var_ind+1+k_window]
                if sum(left_points) in [len(left_points)-1, len(left_points)] and sum(right_points) in [-len(right_points), -len(right_points)+1]:
                    if band[band_var_ind] == max(band[:band_var_ind+1+k_window]):
                        local_maximums_ind_list.append(band_var_ind)
            elif band_var_ind > right_threshold:
                left_points = band_variation_list[band_var_ind-k_window:band_var_ind]
                right_points = band_variation_list[band_var_ind:]
                if right_points:
                    right_points = right_points[1:]
                if sum(left_points) in [k_window-1, k_window] and sum(right_points) in [-len(right_points), -len(right_points)+1]:
                    if band[band_var_ind] == max(band[band_var_ind-k_window:]):
                        local_maximums_ind_list.append(band_var_ind)
            else:
                left_points = band_variation_list[band_var_ind-k_window:band_var_ind]
                right_points = band_variation_list[band_var_ind+1:band_var_ind+1+k_window]
                if sum(left_points) in [k_window-1, k_window] and sum(right_points) in [-k_window, -k_window+1]:
                    if band[band_var_ind] == max(band[band_var_ind-k_window:band_var_ind+1+k_window]):
                        local_maximums_ind_list.append(band_var_ind)
        
        
        local_maximums_ind_list_ = []
        for local_maximums_ind in sorted(local_maximums_ind_list):
            if local_maximums_ind-1 not in local_maximums_ind_list_:
                local_maximums_ind_list_.append(local_maximums_ind)
        #print(sorted(local_maximums_ind_list), local_maximums_ind_list_)
        local_maximums_ind_list = local_maximums_ind_list_
        
        local_maximums_ind_list = sorted(local_maximums_ind_list, key=lambda k_ind: band[k_ind], reverse=True)
        #print(band_variation_list)
        
        return band_variation_list, local_maximums_ind_list

    @classmethod
    def find_local_minimums_of_a_band(cls, band, k_window):
        band_variation_list = [0]
        for prev_e, cur_e in zip(band[0:-1], band[1:]):
            if prev_e > cur_e:
                band_variation_list.append(-1)
            elif prev_e < cur_e:
                band_variation_list.append(1)
            else:
                band_variation_list.append(0)
                
        local_minimums_ind_list = []
        right_threshold = len(band)-k_window-1
        for band_var_ind, band_var in enumerate(band_variation_list):
            if band_var_ind < k_window:
                left_points = band_variation_list[:band_var_ind]
                right_points = band_variation_list[band_var_ind+1:band_var_ind+1+k_window]
                #print(len(left_points), cb_var, len(right_points), [-len(left_points), -len(left_points)+1], [k_window-1, k_window])
                if sum(left_points) in [-len(left_points), -len(left_points)+1] and sum(right_points) in [k_window-1, k_window]:
                    if band[band_var_ind] == min(band[:band_var_ind+1+k_window]):
                        local_minimums_ind_list.append(band_var_ind)
            elif band_var_ind > right_threshold:
                left_points = band_variation_list[band_var_ind-k_window:band_var_ind]
                right_points = band_variation_list[band_var_ind:]
                if right_points:
                    right_points = right_points[1:]
                if sum(left_points) in [-k_window, -k_window+1] and sum(right_points) in [len(right_points)-1, len(right_points)]:
                    if band[band_var_ind] == min(band[band_var_ind-k_window:]):
                        local_minimums_ind_list.append(band_var_ind)
            else:
                left_points = band_variation_list[band_var_ind-k_window:band_var_ind]
                right_points = band_variation_list[band_var_ind+1:band_var_ind+1+k_window]
                if sum(left_points) in [-k_window, -k_window+1] and sum(right_points) in [k_window-1, k_window]:
                    if band[band_var_ind] == min(band[band_var_ind-k_window:band_var_ind+1+k_window]):
                        local_minimums_ind_list.append(band_var_ind)
        
        
        local_minimums_ind_list_ = []
        for local_minimum_ind in sorted(local_minimums_ind_list):
            if local_minimum_ind -1 not in local_minimums_ind_list_:
                local_minimums_ind_list_.append(local_minimum_ind)
        #print(sorted(local_minimums_ind_list), local_minimums_ind_list_)
        local_minimums_ind_list = local_minimums_ind_list_
        
        local_minimums_ind_list = sorted(local_minimums_ind_list, key=lambda k_ind: band[k_ind])
        
        return band_variation_list, local_minimums_ind_list

    
    
    
    #def calculate_effective_mass(self, no_of_local_extrema=2, k_window=5):
    #    pass
    
    def calculate_effective_mass(self,k_ind, band, append_k_data_to_output=False):
        k_window = self.k_window
        high_sym_kpoint_index_list = self.get_high_symmetry_kpoint_index_list()
        polynomial_fit_list = []
        if k_ind < k_window:
            #only fit the right half parabola
            k_list, y_list = VasprunXml.complement_half_parabola(self.kpath[k_ind:k_ind+1+k_window], band[k_ind:k_ind+1+k_window])
            p0, p1, p2 = np.polyfit(k_list, y_list, deg=2)
            polynomial_fit_list.append([(p0, p1, p2), self._effective_mass_from_p0(p0)])
            if append_k_data_to_output:
                polynomial_fit_list[-1].append(k_list)
        elif k_ind > len(self.kpath) - k_window-1:
            #only fit the left half parabola
            k_list, y_list = VasprunXml.complement_half_parabola(self.kpath[k_ind-k_window:k_ind+1], band[k_ind-k_window:k_ind+1])
            p0, p1, p2 = np.polyfit(k_list, y_list, deg=2)
            polynomial_fit_list.append([(p0, p1, p2), self._effective_mass_from_p0(p0)])
            if append_k_data_to_output:
                polynomial_fit_list[-1].append(k_list)
        else:
            if k_ind in high_sym_kpoint_index_list:
                left_k_list, left_y_list = VasprunXml.complement_half_parabola(self.kpath[k_ind-k_window:k_ind+1], band[k_ind-k_window:k_ind+1], left_is_missing=False)
                left_p0, left_p1, left_p2 = np.polyfit(left_k_list, left_y_list, deg=2)
                right_k_list, right_y_list = VasprunXml.complement_half_parabola(self.kpath[k_ind+1:k_ind+2+k_window], band[k_ind+1:k_ind+2+k_window])
                right_p0, right_p1, right_p2 = np.polyfit(right_k_list, right_y_list, deg=2)
                polynomial_fit_list.append([(left_p0, left_p1, left_p2), self._effective_mass_from_p0(left_p0)])
                if append_k_data_to_output:
                    polynomial_fit_list[-1].append(left_k_list)
                polynomial_fit_list.append([(right_p0, right_p1, right_p2), self._effective_mass_from_p0(right_p0)])
                if append_k_data_to_output:
                    polynomial_fit_list[-1].append(right_k_list)
            else:
                p0, p1, p2 = np.polyfit(self.kpath[k_ind-k_window:k_ind+1+k_window], band[k_ind-k_window:k_ind+1+k_window], deg=2)
                polynomial_fit_list.append([(p0, p1, p2), self._effective_mass_from_p0(p0)])
                if append_k_data_to_output:
                    polynomial_fit_list[-1].append(self.kpath[k_ind-k_window:k_ind+1+k_window])
        
        return polynomial_fit_list
    
    def plot_fitted_parabola_in_band_str(self, k_ind_list=[], band_list=[], y_lim=[-3, 3], figure_name="total.png", which_spin="both"):
        kpath, energies = self.get_band_str_data()
        x_lim = [min(kpath), max(kpath)]
        if not self.set_fermi_to_zero:
            fermi_level = self.get_fermi_level()
            y_lim = [y_-fermi_level for y_ in y_lim]
        
        plt.cla()
        
        for spin_ind in range(len(energies)):
            for sub_band in energies[spin_ind]:
                if spin_ind == 0:
                    if which_spin in ["both", "up"]:
                        plt.plot(kpath, sub_band, color="red")
                else:
                    if which_spin in ["both", "down"]:
                        plt.plot(kpath, sub_band, color="blue")
                
        #high-symmetry k-point
        high_symmetry_kpoint_ind_list = self.get_high_symmetry_kpoint_index_list()        
        x_ticks_coords, x_ticks_labels = [], []
        for high_sym_kpoint_ind in high_symmetry_kpoint_ind_list:
            plt.plot([self.kpath[high_sym_kpoint_ind], self.kpath[high_sym_kpoint_ind]], y_lim, color='grey', alpha=0.5)
            x_ticks_coords.append(self.kpath[high_sym_kpoint_ind])
            x_label = "(%.3f,%.3f,%.3f)" % tuple(self.kpoint_list[high_sym_kpoint_ind])
            x_ticks_labels.append(x_label.replace(".000", ""))
        plt.xticks(x_ticks_coords, x_ticks_labels, rotation=45, fontsize="x-small")
        if self.set_fermi_to_zero:
            plt.plot(x_lim, [0, 0], color='grey', alpha=0.5, linestyle="-.")
        else:
            plt.plot(x_lim, [fermi_level, fermi_level], color='grey', alpha=0.5, linestyle="-.")
            
        for k_ind, band in zip(k_ind_list, band_list):
            polynomial_fit_list = self.calculate_effective_mass(k_ind=k_ind, band=band, append_k_data_to_output=True)
            if len(polynomial_fit_list) >= 1:
                polynomial_fit = polynomial_fit_list[0]
                x_list = polynomial_fit[-1]
                y_list = [self.quadratic_formula(polynomial_fit[0], x_) for x_ in x_list]
                plt.plot(x_list, y_list, linewidth=0.8, color="green")
                if len(polynomial_fit_list) == 1:
                    plt.text(x=sum(x_list)/len(x_list), y=min(y_list)-0.1*abs(min(y_list)), s="%.3fm$_e$" % polynomial_fit[1], rotation=90,fontsize="small", alpha=0.8, color="green", verticalalignment="top", horizontalalignment="center")
                else:
                    plt.text(x=min(x_list), y=min(y_list)-0.1*abs(min(y_list)), s="%.3fm$_e$" % polynomial_fit[1], rotation=90,fontsize="small", alpha=0.8, color="green", verticalalignment="top", horizontalalignment="right")
            if len(polynomial_fit_list) >= 2:
                polynomial_fit = polynomial_fit_list[1]
                x_list = polynomial_fit[-1]
                y_list = [self.quadratic_formula(polynomial_fit[0], x_) for x_ in x_list]
                plt.plot(x_list, y_list, linewidth=0.4, color="orange")
                plt.text(x=max(x_list), y=min(y_list)-0.1*abs(min(y_list)), s="%.3fm$_e$" % polynomial_fit[1], rotation=90, fontsize="small", alpha=0.8, color="orange", verticalalignment="top", horizontalalignment="left")
                
                
                
        plt.legend()
        plt.xlim(x_lim)
        plt.ylim(y_lim)
        plt.tight_layout()
        plt.tick_params(axis="both", which="both", direction="in")
        plt.ylabel("Energy (eV)")
        plt.savefig(figure_name, format="png", dpi=1000)
    
    def quadratic_formula(self, coefficients, x):
        return coefficients[0]*x**2 + coefficients[1]*x + coefficients[2]
    
    def _effective_mass_from_p0(self, p0):
        """
        Suppose a local extremum has a energy E0 at k0, the band dispersion around the local extremum can harmonically 
        approximated as 
                        E(k)=E0+[h_bar*(k-k0)]**2/2m, 
        where h_bar is the reduced Planck constant and m is the effective
        mass of holes or electrons around the local extremum. We normally compare the effective mass with the free electron
        mass m_e. So the above equation can be modified as:
                        E(k)=E0+[h_bar*(k-k0)]**2/[2*m_e*(m/m_e)]
                            =E0+[h_bar*c*(k-k0)]**2/[2*m_e_c**2*(m/m_e)]
        Note that in VASP, 2*pi is factored out of the reciprocal lattice vector, i.e. (a_i, b_i)=1. So h_bar here is 
        actually h. And note that hc = 1.23984193 eV*miu_m and m_e*c**2=0.5109989461 MeV. The equation can be smplified as
                        E(k)=E0+1.50412e-18 * (k-k0)**2 / (m/m_e)
        Because in VASP, the unit of k is angstrom, let's change it to m:
                        E(k)=E0+150.412*(k-k0)**2/(m/m_e)
        The input p0 is the coefficient of k**2 term through the parabola fitting around a local extremum. So we may have
                        m/m_e = 150.412/p0
        """
        return 150.412/p0
        
    @classmethod
    def complement_half_parabola(cls, x_list, y_list, left_is_missing=True):
        if left_is_missing:
            central_x = x_list[0]
            missing_half_x_list = list(reversed([2*central_x - x for x in x_list[1:]]))
            missing_half_y_list = list(reversed(y_list[1:]))
            return missing_half_x_list + x_list, missing_half_y_list + y_list
        else:
            central_x = x_list[-1]
            missing_half_x_list = list(reversed([2*central_x - x for x in x_list[:-1]]))
            missing_half_y_list = list(reversed(y_list[:-1]))
            return x_list + missing_half_x_list, y_list + missing_half_y_list
                
                
    def get_NBANDS(self):
        with open(os.path.join(self.cal_loc, "vasprun.xml"), "r") as f:
            for line in f:
                if "<i type=\"int\" name=\"NBANDS\">" in line:
                    return int(line.split(">")[1].split("<")[0])
    
    def get_conduction_band(self):
        pass
    
def polynomial(kpath, k_ind, k_window, p0, p1, p2):
    if k_ind < k_window:
        k_list = kpath[:k_window+k_ind+1]
    elif k_ind + k_window > len(kpath):
        k_list = kpath[k_ind-k_window:]
    else:
        k_list = kpath[k_ind-k_window:k_ind+1+k_window]
    return k_list, [p0*k**2+ p1*k + p2 for k in k_list]


# In[4]:


if __name__ == "__main__":
    vasprun = VasprunXml("./", k_window=int(input("How many k-pionts should be included in the half parabola next to the local extremum for effective mass estimation: ")))
    y_lim = [-3, 3]
    vb, cb = vasprun.get_edge_band()
    vb_var, vb_extrema_k_ind_list, cb_var, cb_extrema_k_ind_list = vasprun.find_local_extrema_of_edge_band()
    k_ind_list, band_list = [], []
    
    for k_ind in vb_extrema_k_ind_list[0]:
        k_ind_list.append(k_ind)
        band_list.append(vb[0])
    for k_ind in cb_extrema_k_ind_list[0]:
        k_ind_list.append(k_ind)
        band_list.append(cb[0])
        
    if len(cb) == 1:
        vasprun.plot_fitted_parabola_in_band_str(k_ind_list=k_ind_list, band_list=band_list, y_lim=y_lim, figure_name="total.png", which_spin="both")
    else:
        vasprun.plot_fitted_parabola_in_band_str(k_ind_list=k_ind_list, band_list=band_list, y_lim=y_lim, figure_name="spin_up.png", which_spin="up")
        
        down_k_ind_list, down_band_list = [], []
        for k_ind in vb_extrema_k_ind_list[1]:
            down_k_ind_list.append(k_ind)
            down_band_list.append(vb[1])
        for k_ind in cb_extrema_k_ind_list[1]:
            down_k_ind_list.append(k_ind)
            down_band_list.append(cb[1])
            
        vasprun.plot_fitted_parabola_in_band_str(k_ind_list=down_k_ind_list, band_list=down_band_list, y_lim=y_lim, figure_name="spin_down.png", which_spin="down")

        total_k_ind_list = k_ind_list + down_k_ind_list
        total_band_list = band_list + down_band_list
        vasprun.plot_fitted_parabola_in_band_str(k_ind_list=total_k_ind_list, band_list=total_band_list, y_lim=y_lim, figure_name="total.png", which_spin="both")

