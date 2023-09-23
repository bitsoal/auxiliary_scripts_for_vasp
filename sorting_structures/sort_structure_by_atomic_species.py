#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pymatgen.core import Structure

import sys


# In[3]:


filename = sys.argv[1]


# In[8]:


struct = Structure.from_file(filename)

sorted_sites = sorted(struct, key=lambda site: site.species_string)
sorted_struct = Structure(lattice=struct.lattice.matrix, 
                          species=[site.species_string for site in sorted_sites], 
                          coords=[site.frac_coords for site in sorted_sites],
                          coords_are_cartesian=False)


# In[12]:


sorted_struct.to(filename="sorted_POSCAR", fmt="poscar")
print("Sort {} and save the sorted structure into file sorted_POSCAR".format(filename))

