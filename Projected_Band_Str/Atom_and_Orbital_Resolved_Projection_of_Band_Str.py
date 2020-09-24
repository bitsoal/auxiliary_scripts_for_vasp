#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, re, math

#import matplotlib.pyplot as plt
#%matplotlib inline


# In[2]:


class PROCAR():

    """
        This class aims to extract the calculated band structures from PROCAR and batch process projections on either selected atoms or orbitals.
        OUTCAR is needed to extract additional information:
                    1. the reciprocal lattice vectors and the high symmetry k-path for scalar k-path construction
                    2. E-fermi --> It will be set to zero in the output band structure
                    3. ISPIN & LOSOBITAL
        This class should be able to process both collinear and non-collinear calculations. The former could be either non-spin-polarized or spin-polarized.
    """
    
    def __init__(self, cal_loc='.', procar_filename="PROCAR", outcar_filename="OUTCAR"):
        self.cal_loc = cal_loc
        self.procar_filename = procar_filename
        self.outcar_filename = outcar_filename
        self.cal_summary_dict = PROCAR.obtain_cal_summary_from_outcar(cal_loc=cal_loc, outcar_filename=outcar_filename)
        
    @classmethod
    def obtain_cal_summary_from_outcar(cls, cal_loc=".", outcar_filename="OUTCAR"):
        outcar_filename = os.path.join(cal_loc, outcar_filename)
        
        cal_summary_dict = {}
        with open(outcar_filename, "r") as outcar_f:
            for line in outcar_f:
                line = line.strip()
                if line.startswith("ISPIN"):
                    cal_summary_dict["ispin"] = int(line.split()[2])
                elif line.startswith("LORBIT"):
                    cal_summary_dict["lorbit"] = int(line.split()[2])
                elif line.startswith("LSORBIT"):
                    cal_summary_dict["lsorbit"] = True if "t" in line.split()[2].lower() else False
                elif "NBANDS" in line and "number of bands" in line:
                    cal_summary_dict["nbands"] = int(line.split("NBANDS=")[-1].strip())
                    cal_summary_dict["nkpoints"] = int(line.split("k-points in BZ")[0].split("=")[1].strip())
                elif "NIONS" in line:
                    cal_summary_dict["nions"] = int(line.split("NIONS =")[-1].strip())
                elif line.startswith("E-fermi"):
                    cal_summary_dict["e_fermi"] = float(line.split()[2])
                
        
        with open(outcar_filename, "r") as outcar_f:
            #Retrieve the reciprocal lattice vectors
            for line in outcar_f:
                if "direct lattice vectors                 reciprocal lattice vectors" in line:
                    break
            direct_latt_vectors, rec_latt_vectors = [], []
            for i in range(3):
                latt_vec = next(outcar_f).replace("-", " -").split()
                direct_latt_vectors.append([float(dir_) for dir_ in latt_vec[:3]])
                rec_latt_vectors.append([float(rec_) * 2 * math.pi for rec_ in latt_vec[-3:]])
            cal_summary_dict["rec_latt_vectors"] = rec_latt_vectors
            cal_summary_dict["dir_latt_vectors"] = direct_latt_vectors
            print("Pls note that the reciproal lattice vectors in OURCAR omit 2*pi. We multiply them by 2*pi.")
            
            #Retieve k-point lines
            for line in outcar_f:
                line = line.strip()
                if line.startswith("k-point  1 :") and "plane waves:" in line:
                    break
            kpoint_lines = [line]
            for i in range(cal_summary_dict["nkpoints"]-1):
                kpoint_lines.append(next(outcar_f).strip())
        cal_summary_dict["kpoint_lines"] = kpoint_lines
        
        #Parse k-points from the retrieved k-point lines
        kpoint_vectors = []
        for kpoint_line in cal_summary_dict["kpoint_lines"]:
            kpoint_line = kpoint_line.split(":")[1]
            kpoint_vec = re.findall("[-]?[0-9\.]+", kpoint_line)
            kpoint_vectors.append([float(k_) for k_ in kpoint_vec])
        cal_summary_dict["kpoint_vectors"] = kpoint_vectors
        
        #Convert k-points to a scalar k-path
        cal_summary_dict["scalar_kpath"] = PROCAR.cal_scalar_kpath_from_kpoint_vectors(kpoint_vectors=kpoint_vectors, rec_latt_vectors=rec_latt_vectors)
        
        avaiable_orbitals = []
        with open(os.path.join(cal_loc, "PROCAR"), "r") as procar_f:
            for line in procar_f:
                if line.startswith("ion"):
                    avaiable_orbitals = line.strip().split()[1:]
                    break
        cal_summary_dict["available_orbitals"] = avaiable_orbitals
        print("Projection on the following orbitals of each atom is available by parsing PROCAR: \n" + "\t".join(avaiable_orbitals))
        
        return cal_summary_dict
        
        
        
    @classmethod
    def cal_scalar_kpath_from_kpoint_vectors(cls, kpoint_vectors, rec_latt_vectors):
        kpath = [0]
        for k0, k1 in zip(kpoint_vectors[0:-1], kpoint_vectors[1:]):
            frac_diff_kpoint_vec = [k1_ - k0_ for k1_, k0_ in zip(k1, k0)]
            cartesian_diff_kpoint_vec = [0, 0, 0]
            for dk_i, rec_latt_vec_i in zip(frac_diff_kpoint_vec, rec_latt_vectors):
                for rec_latt_vec_ij_ind, rec_latt_vec_ij in enumerate(rec_latt_vec_i):
                    cartesian_diff_kpoint_vec[rec_latt_vec_ij_ind] += dk_i * rec_latt_vec_ij
            dkpath = pow(sum([dk_ij * dk_ij for dk_ij in cartesian_diff_kpoint_vec]), 0.5)
            kpath.append(dkpath+kpath[-1])
        return kpath
    
    def read_target_projections(self):
        """
        Read and parse target projections from a file named "target_projections" under cal_loc.
        
        In target_projections, lines starting with "target_proj=" define projections on selected atoms and orbitals. Note that each line only defines ONE target projection.
            Format: "projection_type|ions:orbitals|short_name"
                -- "projection_type": "tot", "mx", "my" or "mz". CAN ONLY select one of them.
                -- "ions": 1-based ion index. Separate them with comma if it is a combination.
                -- "orbitals": any orbital appearing in PROCAR ("s", "py", "pz", "px", ..., "tot"). Separate them with comma if it is a combination.
                -- "short_name": the short name of the target projection. Will be used in the output file of the projected band structure.
            White spaces and what follows "#" will be ignored.
            e.g.
                Example 1: target_proj = tot|1:s|atom_1_s <--> project on the s orbital of ion/atom 1
                Example 2: target_proj = tot|1,2,4:s|atoms_1_2_4_s <--> project on the s orbital of ion/atom 1, 2 and 4
                Example 3: target_proj = tot|1:s + 1:px,py,pz|atom_1_s_p <--> project on s, px, py and pz orbitals of ion/atom 1
                Example 4: target_proj = mx|1-3:pz,py|atom_1_to_3_pz_py <--> project on pz and py orbitals of ion/atom 1, 2 and 3 along mx direction.
        Return a 2-element tuple:
            -- 1st element: a list of target projections. The same format as the input arguement target_projections of method get_projections_on_atom_and_orbital
            -- 2nd element: a list of short names associated with the target projections in the 1st element.
        """
        
        target_projections_filename = os.path.join(self.cal_loc, "target_projections")
        if not os.path.isfile(target_projections_filename):
            print("\n" + self.read_target_projections.__doc__ + "\n")
            with open(target_projections_filename, "w") as target_projection_f:
                target_projection_f.write("#Projection on the following orbitals of each atom is available by parsing PROCAR:\n")
                target_projection_f.write("#"+"\t".join(self.cal_summary_dict["available_orbitals"])+"\n")
            raise Exception("Cannot find file 'target_projections' defining the target projections. We create an empty 'target_projections' for you .Please follow the format guidelines above and define your target projections in it")
            
        target_projections_str = []
        proj_ind_list = []
        with open(target_projections_filename, "r") as f:
            for line in f:
                line = line.split("#")[0].strip().replace(" ", "")
                if line.startswith("target_proj="):
                    assert line.count("=") == 1, "Each line should define ONLY ONE target projection. Try to change the line below\n%s" % line
                    target_projections_str.append(line.split("=")[1])
        
        target_projections = []
        short_projection_names = []
        for proj_str in target_projections_str:
            proj_type, ion_orbital_pairs, short_name = proj_str.split("|")
            if proj_type in ["mx", "my", "mz"]:
                short_projection_names.append("+" + short_name)
                short_projection_names.append("-" + short_name)
            else:
                short_projection_names.append(short_name)
            target_proj = []
            for ion_orbital in ion_orbital_pairs.split("+"):
                try:
                    ion_part, orbital_part = ion_orbital.split(":")
                    ion_ind_list = []
                    for ion_ in ion_part.split(","):
                        if "-" in ion_:
                            start_ion_ind, end_ion_end = [int(ion_ind) for ion_ind in ion_.split("-")]
                            ion_ind_list.extend(range(start_ion_ind, end_ion_end+1))
                        else:
                            ion_ind_list.append(ion_)
                    target_proj.extend([str(ion_ind) + "-" + orbital for ion_ind in ion_ind_list for orbital in orbital_part.split(",")])
                except:
                    print("\nFail to parse the target projection below")
                    print(proj_str)
                    print("Try to correct it by refereing to the guidelines below")
                    print(self.read_target_projections.__doc__)
                    raise
            new_target_proj = []
            for ion_orbital_ in target_proj:
                if ion_orbital_ in new_target_proj:
                    print("Warning: We found you set %s more than once when reading the line below from file 'target_projections'. We will just it once." % ion_orbital_)
                    print(proj_str)
                else:
                    new_target_proj.append(ion_orbital_)
            target_projections.append(new_target_proj + [proj_type])
            
            
        return target_projections, short_projection_names
    
    def write_projected_band_str(self, energies, projections, projection_details, short_projection_names):
        """
        Input arguements:
            energies, projections, projection_details: exactly the output of get_projections_on_atom_and_orbital
            short_projection_names: a list of short names associated with the projections. 
                                    It is also the 2nd element of the 2-element output tuple of method read_target_projections
        """
        if len(energies) == 2: #collinear spin-polarized calculation
            output_filename_list = ["projected_band_str_spin_up_channel.dat", "projected_band_str_spin_down_channel.dat"]
        else:
            output_filename_list = ["projected_band_str.dat"]
        output_filename_list = [os.path.join(self.cal_loc, filename) for filename in output_filename_list]
        
        no_of_projections = len(projection_details)
        for spin_ind in range(len(energies)):
            with open(output_filename_list[spin_ind], "w") as output_f:
                #Write header
                title_line = "#1_based_band_ind\t1_based_k_ind\tscalar_k_path\tenergy\t"
                for short_name, projection_detail in zip(short_projection_names, projection_details):
                    output_f.write("#%s: %s\n" % (short_name, projection_detail))
                    title_line += "%s\t%s_weight\t" % (short_name, short_name)
                output_f.write(title_line + "\n")
                #Start to write band structures and target projections
                for band_ind in range(self.cal_summary_dict["nbands"]):
                    for k_ind in range(self.cal_summary_dict["nkpoints"]):
                        output_f.write("%d\t%d\t" % (band_ind+1, k_ind+1))
                        output_f.write("%f\t%f\t" % (self.cal_summary_dict["scalar_kpath"][k_ind], energies[spin_ind][band_ind][k_ind]))
                        for proj_ind in range(no_of_projections):
                            output_f.write("%f\t%f\t" % (energies[spin_ind][band_ind][k_ind], projections[spin_ind][band_ind][k_ind][proj_ind]))
                        output_f.write("\n")
                    output_f.write("\n")
        print("Done")
        
    def get_projections_on_atom_and_orbital(self, target_projections):
        """
        target_projections is a list of lists. Each list entry of target_projections is an individual projection 
            on selected ions' orbitals. 
        Each list entry of target_projections has a format below:
            Let's denote the projection on "orbital" of the n-th atom as P("orbital"@n), where n is 1-based and "orbital" could be any ortibal appearing in PROCAR, inclusive of "tot"
            example for collinear cal: ["1-s", "1-py, "tot"] <--> P(s@1) + P(py@1)
            example for collinear cal: ["1-x2-y2", "1-dz2", "tot"] <--> P(t2g states@1)
            example for collinear cal: ["1-tot", "2-tot", "tot"] <--> P(all orbitals@1) + P(all orbitals@2)
            example for collinear cal: ["tot-s", "tot"] <--> P(s@all atoms)
            example for non-collinear cal: ["1-s", , "2-dz2", "mx"] <--> P(s@1) + P(py@1) + (x2-y2@2) along the x-axis direction.
            Note that the last entry denotes the projection type:
                For collinear cal, it is always "tot";
                For non-collinear cal, it can be "tot", "mx", "my" or "mz"
                
        Return a 3-element tuple:
            - 1st element:
                    If ispin = 1 or lsorbit = True:
                        a 1-element tuple.
                        The only entry is 2d eigenvalues w.r.t. E-fermi: E[band_ind][k_ind]
                    else if it is a collinear spin-polarized calculation:
                        a 2-element tuple of 2d eigenvalues: (spin-up 2d eigenvalues, spin-down 2d eigenvalues)
            - 2nd element: 
                    if ispin = 1 or lsorbit = True:
                        a 1-element tuple.
                        The only entry is 3d projections: projections[band_ind][k_ind][proj_ind]
                        k_ind, band_ind and proj_ind are 0-based integer indexes
                    else if it is a collinear spin-polarized calculatoin:
                        a 2-element tuple of 3d projectoins: (spin-up 3d projections, spin-down 3d projections)
            - 3rd element:
                    a list of projection details in accordace with the 2nd element projections
        """
        
        e_fermi = self.cal_summary_dict["e_fermi"]
        ispin = self.cal_summary_dict["ispin"]
        lsorbit = self.cal_summary_dict["lsorbit"] # True for soc calculations; False otherwise
        if lsorbit:
            for target_proj in target_projections:
                assert target_proj[-1] in ["tot", "mx", "my", "mz"], "The last entry of each target projection list should be 'tot', 'mx', 'my' or 'mz' for non-collinear calculations"
        else:
            for target_proj in target_projections:
                assert target_proj[-1] == "tot", "The last entry of each target projection list shoud be 'tot' for collinear calculations"
        
        projection_details = []
        for target_proj in target_projections:
            if target_proj[-1] == "tot":
                projection_details.append("tot: " + ",".join(target_proj[:-1]))
            else:
                projection_details.append("+" + target_proj[-1] + ": " + ",".join(target_proj[:-1]))
                projection_details.append("-" + target_proj[-1] + ": " + ",".join(target_proj[:-1]))
        
        #projections[band_ind][k_ind][proj_ind]
        projections = [[] for i in range(self.cal_summary_dict["nbands"])]
        #energies[band_ind][k_ind]. 
        energies = [[] for i in range(self.cal_summary_dict["nbands"])]
        
        no_of_kpoint_0_and_band_0 = 0
        for bloch_state_dict in self.get_bloch_state_iterator():
            #convert 1-based indexes to 0-based indexes
            k_ind = bloch_state_dict["k_ind"] - 1
            band_ind = bloch_state_dict["band_ind"] - 1
            if k_ind == 0 and band_ind == 0:
                no_of_kpoint_0_and_band_0 += 1
            if no_of_kpoint_0_and_band_0 == 2:
                no_of_kpoint_0_and_band_0 += 1 #deactivate this tag
                spin_up_projections = projections
                spin_up_energies = energies
                #projections and energies below essentially correspond to the spin-down channel
                projections = [[] for i in range(self.cal_summary_dict["nbands"])]
                energies = [[] for i in range(self.cal_summary_dict["nbands"])]
                
            #set E_fermi to 0, i.e. the reference energy
            energy = bloch_state_dict["energy"] - e_fermi
            energies[band_ind].append(energy)
            #For double check ONLY
            #print(k_ind, band_ind, no_of_kpoint_0_and_band_0)
            assert energies[band_ind][k_ind] == energy, "energies is not appended to the right place (band_ind, k_ind)=(%d, %d)" % (band_ind, k_ind)
            
            projection_list = []
            for target_proj in target_projections:
                proj_type = target_proj[-1]
                proj = [bloch_state_dict["bloch_state_dict"][proj_type][ion_orbital_key] for ion_orbital_key in target_proj[:-1]]
                if proj_type == "tot":
                    projection_list.append(sum(proj))
                else:
                    projection_list.append(sum([proj_ for proj_ in proj if proj_>=0]))
                    projection_list.append(sum([proj_ for proj_ in proj if proj_<0]))
            
            #append it to projections and energies
            projections[band_ind].append(projection_list)
            #For double check ONLY
            assert projections[band_ind][k_ind] == projection_list, "projection_list is not appended to the right place (band_ind, k_ind)=(%d, %d)" % (band_ind, k_ind)
            
            
        if no_of_kpoint_0_and_band_0 == 1: #collinear non-spin-polarized cal or non-collinear cal
            return (energies, ), (projections, ), projection_details
        else: #collinear spin-polarized cal
            #Here energies and projections correspond to the spin-down channel
            return (spin_up_energies, energies), (spin_up_projections, projections), projection_details
    

    def get_bloch_state_of(self, k_ind, band_ind, spin_channel="up"):
        """
        k_ind and band_ind are the 1-based index of k-point and band
        spin_channel could be either "up" or "down". "down" is meaningless for non-spin-polarized or non-collinear calculations
        
        return the bloch state of k_ind and band_ind.
        
        The returned bloch state has the exactly same format of the entry of get_bloch_state_iterator()
        """
        assert k_ind <= self.cal_summary_dict["nkpoints"], "the 1-based k_ind should be in between 1 and %d, inclusive of both ends" % self.cal_summary_dict["nkpoints"]
        assert band_ind <= self.cal_summary_dict["nbands"], "the 1-based band_ind should be in between 1 and %d, inclusive of both ends" % self.cal_summary_dict["nbands"]
        assert spin_channel in ["up", "down"], "spin_channel must be either 'up' or 'down'"
        if self.cal_summary_dict["ispin"] == 2 and not self.cal_summary_dict["lsorbit"]:
            target_spin_channel = spin_channel
        else:
            target_spin_channel = "up"
        
        current_spin_channel = "up"
        for bloch_state in self.get_bloch_state_iterator(parse_bloch_state_block=False):
            if bloch_state["k_ind"] == k_ind and bloch_state["band_ind"] == band_ind:
                if current_spin_channel == target_spin_channel:
                    bloch_state["bloch_state_dict"] = self.parse_a_bloch_state_block(bloch_state["bloch_state_block_lines"])
                    return bloch_state
                else:
                    current_spin_channel = "down"
                
        
    def get_bloch_state_iterator(self, parse_bloch_state_block=True):
        """
        Read PROCAR and return a bloch state iterator.
        Each entry of the iterator is a dictionary:
            - "procar_line_ind" (intger): the 1-based index of the last line of the bloch state
            - "k_ind" (integer): 1-based k-point index as shown in PROCAR
            - "band_ind" (integer): 1-based band index as shown in PROCAR
            - "energy" (float): energy of the bloch state
            - "energy_str" (str): a str version of energy of the bloch state
            - "bloch_state_block_lines" (list): the bloch state block
            - "bloch_state_dict" (dict; only for parse_bloch_state_block=True): 
                    if the calculatoin takes into account the spin-orbital coupling:
                        the dictionary contains four sub-directories storing the projection of the total magnetization on each ion and orbital, 
                        and the projectoin of the magnetization along the x, y and z axis.
                        The keys referring to the sub-directionaries: "tot", "mx", "my" and "mz"
                    if the calculation does not take into account the spin-orbital coupling:
                        the dictionary contains only one sub-directory storing the projection on each ion and orbital.
                        The key referring to the sub-directory: "tot"
                    The format of the sub-directory/sub-directories:
                        key: "ion-orbital", where "ion" is the 1-based ion index, which is the same as the first column in PROCAR.
                                "orbital" is the orbital name, e.g. s, py, pz, px, dxy, dyz. The orbital name must be the same as that shown in PROCAR
                        value: the projection of the bloch state on "orbital" of "ion"    
        """
        procar_filename = os.path.join(self.cal_loc, self.procar_filename)
        lsorbit = self.cal_summary_dict["lsorbit"] # True for soc calculations; False otherwise
        
        #for soc cal, the bloch state block has four sub-block, i.e. projection on tot, mx, my and mz.
        #     each sub-block ends with a line which starts with "tot"
        #for non-sco cal, there is only one sub-block, namely "tot"
        max_no_of_tot = 4 if lsorbit else 1
        
        with open(procar_filename, "r") as procar_f:
            is_bloch_state_block = False
            bloch_state_block_lines = []
            
            for line_ind, line in enumerate(procar_f):
                line = line.strip()
                if line == "":
                    continue
                elif is_bloch_state_block:
                    bloch_state_block_lines.append(line)
                    if line.startswith("tot"):
                        no_of_tot += 1
                    if no_of_tot == max_no_of_tot:
                        is_bloch_state_block = False
                        bloch_state_dict["procar_line_ind"] = line_ind
                        bloch_state_dict["bloch_state_block_lines"] = bloch_state_block_lines
                        if parse_bloch_state_block:
                            bloch_state_dict["bloch_state_dict"] = self.parse_a_bloch_state_block(bloch_state_block_lines)
                        yield bloch_state_dict
                elif line.startswith("k-point") and "weight" in line:
                    k_ind = int(line.split(":")[0].split("k-point")[1].strip())                     
                elif line.startswith("band") and "energy" in line:
                    band_ind = int(line.split("#")[0].split("band")[1].strip())
                    energy_str = line.split("#")[1].split("energy")[1].strip()
                    energy = float(energy_str)
                    
                    bloch_state_dict = {"k_ind": k_ind, "band_ind": band_ind, "all_index_is_1_based": ""}
                    bloch_state_dict["energy"] = energy
                    bloch_state_dict["energy_str"] = energy_str
                    #What follows this line is the corresponding bloch state block
                    is_bloch_state_block = True
                    no_of_tot = 0 
                    bloch_state_block_lines = []
    

        
    def parse_a_bloch_state_block(self, bloch_state_block_lines):
        assert bloch_state_block_lines[0].startswith("ion"), "the bloch state block must start with a line 'ion      s     py     pz ...'"
        orbital_name_list = bloch_state_block_lines[0].split()[1:]
        nions = self.cal_summary_dict["nions"]
        #print(bloch_state_block_lines)
        #print()
        
        bloch_state_dict = {}
        if self.cal_summary_dict["lsorbit"]:
            #line index range
            #orbital name list: the 0th line
            #tot: from 1 to 1+nions, inclusive of both ends
            #mx: from 2+ions to 2+2*nions, inclusive of both ends
            #my: from 3+2*nions to 3+3*nions, inclusive of both ends
            #mz: from 4+3*nions to 4+4*nions, inclusive of both ends
            bloch_state_dict["tot"] = self._parse_a_sub_bloch_state_block(bloch_state_block_lines[1:2+nions], orbital_name_list)
            bloch_state_dict["mx"] = self._parse_a_sub_bloch_state_block(bloch_state_block_lines[2+nions:3+2*nions], orbital_name_list)
            bloch_state_dict["my"] = self._parse_a_sub_bloch_state_block(bloch_state_block_lines[3+2*nions:4+3*nions], orbital_name_list)
            bloch_state_dict["mz"] = self._parse_a_sub_bloch_state_block(bloch_state_block_lines[4+3*nions:], orbital_name_list)
        else:
            bloch_state_dict["tot"] = self._parse_a_sub_bloch_state_block(bloch_state_block_lines[1:], orbital_name_list)

        return bloch_state_dict
            
    def _parse_a_sub_bloch_state_block(self, sub_bloch_state_block_lines, orbital_name_list):
        sub_block_dict = {}
        for line in sub_bloch_state_block_lines:
            items = re.findall("[-]?[0-9\.to]+", line)
            assert len(items)-1 == len(orbital_name_list), "Fail to correctly parse the line below %s" % line
            ion_ind = items[0] #This is a 1-based ion index
            for weight, orbital_name in zip(items[1:], orbital_name_list):
                ion_orbital_key = ion_ind + "-" + orbital_name
                sub_block_dict[ion_orbital_key] = float(weight)
        return sub_block_dict


# In[3]:


if __name__ == "__main__":
    procar = PROCAR(cal_loc=".")
    target_projections, short_projection_names = procar.read_target_projections()
    print(target_projections, short_projection_names)
    energies, projections, projection_details = procar.get_projections_on_atom_and_orbital(target_projections)
    procar.write_projected_band_str(energies, projections, projection_details, short_projection_names)


# plt.figure(figsize=(16, 8))
# scalar_kpath = procar.cal_summary_dict["scalar_kpath"]
# plt.subplot(1, 2, 1)
# for band_ind in range(procar.cal_summary_dict["nbands"]):
#     plt.plot(scalar_kpath, energies[0][band_ind], color="gray")
#     plt.scatter(scalar_kpath, energies[0][band_ind], s=[proj[4]*100 for proj in projections[0][band_ind]], color="red")
# plt.ylim([-5, 5])
# plt.xlim([0, max(scalar_kpath)])    
# 
# plt.subplot(1, 2, 2)
# for band_ind in range(procar.cal_summary_dict["nbands"]):
#     plt.plot(scalar_kpath, energies[0][band_ind], color="gray")
#     plt.scatter(scalar_kpath, energies[0][band_ind], s=[abs(proj[5])*100 for proj in projections[0][band_ind]], color="orange")
# plt.ylim([-5, 5])
# plt.xlim([0, max(scalar_kpath)])
# projection_names

# band_ind = 55
# plt.figure(figsize=[16, 6])
# plt.subplot(1, 2, 1)
# plt.plot(scalar_kpath, [proj[0] for proj in projections[0][band_ind]], color="blue", label="Ag 1 4")
# plt.plot(scalar_kpath, [proj[1] for proj in projections[0][band_ind]], color="red", label="Ag 2 3")
# plt.plot(scalar_kpath, [proj[2] for proj in projections[0][band_ind]], color="orange", label="S")
# plt.scatter([scalar_kpath[81], scalar_kpath[161]], [0.3, 0.3])
# plt.xlim([min(scalar_kpath), max(scalar_kpath)])
# plt.legend()
# 
# plt.subplot(1, 2, 2)
# plt.plot(scalar_kpath, [sum(proj) for proj in projections[0][band_ind]], color="black", label="tot")
# plt.scatter([scalar_kpath[81], scalar_kpath[161]], [0.8, 0.8])
# plt.xlim([min(scalar_kpath), max(scalar_kpath)])
