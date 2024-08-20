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

def calculate_ssr(cut_off):
    # Fit the blackbody model to the segment before the cut-off point
    x_bb_segment = x_bb[x_bb <= cut_off]
    y_bb_segment = y_bb[x_bb <= cut_off]
    
    # Skip the fit if there are fewer than 2 points
    if len(x_bb_segment) < 2:
        return np.inf  # Assign a very high SSR value
    
    popt_bb, _ = curve_fit(bb, x_bb_segment, y_bb_segment)
    g_fit, f_fit = popt_bb
    
    # Fit the power law model to the segment after the cut-off point
    x_pl_segment = x_pl[x_pl > cut_off]
    y_pl_segment = y_pl[x_pl > cut_off]
    
    # Skip the fit if there are fewer than 5 points
    if len(x_pl_segment) < 5:
        return np.inf  # Assign a very high SSR value
    
    initial_guess_pl = [-0.0116895, 2.8322e10, 0.840665, 94726.6, 0.596613]
    popt_pl, _ = curve_fit(pl, x_pl_segment, y_pl_segment, p0=initial_guess_pl)
    
    # Calculate residuals and SSR
    residuals_bb = y_bb_segment - bb(x_bb_segment, *popt_bb)
    residuals_pl = y_pl_segment - pl(x_pl_segment, *popt_pl)
    ssr = np.sum(residuals_bb**2) + np.sum(residuals_pl**2)
    return ssr

# Range of cut-off points
cut_off_points = np.linspace(min(x_bb), max(x_pl), 10000)
ssr_values = [calculate_ssr(cut_off) for cut_off in cut_off_points]

# Find the optimal cut-off point
optimal_cut_off = cut_off_points[np.argmin(ssr_values)]

# Plot SSR vs cut-off points
plt.figure(figsize=(10, 6))
plt.plot(cut_off_points, ssr_values, label='SSR vs Cut-off Point')
plt.axvline(optimal_cut_off, color='red', linestyle='--', label=f'Optimal Cut-off: {optimal_cut_off:.3f} keV')
plt.xlabel('Cut-off Point (keV)', fontsize=14)
plt.ylabel('Sum of Squared Residuals (SSR)', fontsize=14)
plt.legend(fontsize=13)
plt.show()


# In[ ]:




