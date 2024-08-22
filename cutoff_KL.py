#!/usr/bin/env python
# coding: utf-8

# In[2]:


#Install Required Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from IPython.display import display, Math
import matplotlib.colors


# Define Model Equations
#Power Law
def pl(x,a, c, d,m, b):
    return d + (a * x - d) / ((1 + (x / c)**b)**m)
# Black Body equation
def bb(x, g, f):
    exponent = f * x
    y = 2*(g)**4 * x**3 / (np.exp(exponent) - 1)
    return y
#Black body with fixed f_col
def bb_fix(x, h):
    exponent = h * x
    y = 2*(1.7)**4 * x**3 / (np.exp(exponent) - 1)
    return y

#Blackbody co-ordinates
x_bb = np.array((3.58, 5, 7,9, 10,11, 13))
y_bb = np.array((10.8, 5.2, 1.2,0.2, 0.095,0.038, 0.01))
#Power Law co-ordinates
x_pl = np.array((3.58, 7, 10, 20, 30, 40, 50))
y_pl = np.array((0.058, 0.067,0.07 ,0.06, 0.052, 0.04, 0.0295))



# In[3]:


from scipy.stats import entropy

def calculate_kl_divergence(cut_off):
    # Segment the data into two parts based on the cut-off
    x_bb_segment = x_bb[x_bb <= cut_off]
    y_bb_segment = y_bb[x_bb <= cut_off]
    
    x_pl_segment = x_pl[x_pl > cut_off]
    y_pl_segment = y_pl[x_pl > cut_off]
    
    # Skip calculation if there are too few points in any segment
    if len(x_bb_segment) < 2 or len(x_pl_segment) < 5:
        return -np.inf  # Assign a very low KL value to avoid choosing this cutoff
    
    # Normalize y-values to create probability distributions
    p_bb = y_bb_segment / np.sum(y_bb_segment)
    p_pl = y_pl_segment / np.sum(y_pl_segment)
    
    # Match the lengths of the segments by interpolation or re-binning if necessary
    if len(p_bb) != len(p_pl):
        p_bb_interp = np.interp(np.linspace(min(x_bb_segment), max(x_bb_segment), len(p_pl)),
                                x_bb_segment, p_bb)
        p_pl_interp = p_pl
    else:
        p_bb_interp = p_bb
        p_pl_interp = p_pl
    
    # Calculate KL divergence
    kl_divergence = entropy(p_bb_interp, p_pl_interp)
    return kl_divergence


# In[10]:


# Range of cut-off points
cut_off_points = np.linspace(min(x_bb), max(x_pl), 10000)
kl_divergence_values = [calculate_kl_divergence(cut_off) for cut_off in cut_off_points]

# Find the optimal cut-off point that maximizes KL divergence
optimal_cut_off_kl = cut_off_points[np.argmax(kl_divergence_values)]

# Plot KL Divergence vs cut-off points
fig, ax = plt.subplots(figsize=(10, 6))

plt.plot(cut_off_points, kl_divergence_values, label='KL Divergence')
plt.axvline(optimal_cut_off_kl, color='red', linestyle='--', label=f'Optimal Cut-off: {optimal_cut_off_kl:.3f} keV')
plt.xlabel('Cut-off Point (keV)', fontsize=14)
plt.ylabel('KL Divergence', fontsize=14)
ax.xaxis.set_tick_params(which='both', labeltop=False)
ax.yaxis.set_tick_params(which='both', labelright=False)
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
plt.legend(fontsize=14)
plt.savefig('cutoff_kl', dpi=300, bbox_inches = 'tight')
plt.show()


# In[ ]:




