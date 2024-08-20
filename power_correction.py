#!/usr/bin/env python
# coding: utf-8

# In[6]:


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
x_bb = (3.58, 5, 7,9, 10,11, 13)
y_bb = (10.8, 5.2, 1.2,0.2, 0.095,0.038, 0.01)
#Power Law co-ordinates
x_pl = (3.58, 7, 10, 20, 30, 40, 50)
y_pl = (0.058, 0.067,0.07 ,0.06, 0.052, 0.04, 0.0295)

#curve fit for bbody
popt_bb1, pcov_bb = curve_fit(bb, x_bb, y_bb)
g_fit, f_fit = popt_bb1
perr = np.sqrt(np.diag(pcov_bb))
g_fit_unc = perr[0]
#curve fit for bbody with fixed cc factor
popt_bb_fix, pcov_bb_fix = curve_fit(bb_fix, x_bb, y_bb)
h_fit = popt_bb_fix

#curve fit for power law
initial_guess_pl = [-0.0116895, 2.8322e10,0.840665 ,94726.6 , 0.596613]
popt_pl1, pcov_pl = curve_fit(pl, x_pl, y_pl, p0=initial_guess_pl)
a_fit, c_fit, d_fit, m_fit, b_fit = popt_pl1

x_fit = np.linspace(3.58, 60, 1000)

def total_fit(x, d, a, c, b, m, f, g):
    pl = d + (a * x - d) / ((1 + (x / c)**b)**m)
    exponent = f * x
    bb = 2*(g)**4 * x**3 / (np.exp(exponent) - 1)
    return pl + bb

x_data = np.array((3.58, 5,7, 10, 13, 20, 30, 40, 50))
y_data = np.array((11, 5.3,1.3, 0.16, 0.07, 0.06, 0.052, 0.04, 0.0295))
popt_t, pcov_t = curve_fit(total_fit, x_data, y_data, p0 =[0.840665, -0.0116895, 2.8322e10, 0.596613, 94726.6, 1.19662891, 1.7])
d_fit, a_fit, c_fit, b_fit, m_fit, f_fit, g_fit = popt_t

perr = np.sqrt(np.diag(pcov_t))
g_fit_unc = perr[6]
x_fit = np.linspace(3.58, 60, 1000)
y_fit_t = total_fit(x_fit, d_fit, a_fit, c_fit, b_fit, m_fit, f_fit, g_fit)

y_fit_bb = bb(x_fit, g_fit, f_fit)

#Blackbody with fixed cc factor plot & Uncertainty
y_fit_bb_fix = bb_fix(x_fit, h_fit)

#Power law Plot
y_fit_pl = pl(x_fit, a_fit, c_fit, d_fit, m_fit, b_fit)

# Axis settings
fig, ax = plt.subplots(figsize=(10, 6))

ax.set_ylim(0.01, 15)
ax.set_xlim(3.58, 60)
ax.set_xscale('log')
ax.set_yscale('log')
ax.scatter(x_data, y_data, label = 'Model data', color = 'black')
ax.plot(x_fit, y_fit_t, linestyle = 'dashdot', label = 'Original model', color = 'red')
ax.plot(x_fit, y_fit_pl+y_fit_bb, linestyle = 'dashed', label = 'Power-corrected Model', color = 'blue')
    
# Adding grid lines
ax.axvline(x=10, color='gray', linestyle='-', linewidth=0.5)
ax.axhline(y=0.1, color='gray', linestyle='-', linewidth=0.5)
ax.axhline(y=1, color='gray', linestyle='-', linewidth=0.5)
ax.axhline(y=10, color='gray', linestyle='-', linewidth=0.5)

ax.xaxis.set_tick_params(which='both', labeltop=False)
ax.yaxis.set_tick_params(which='both', labelright=False)
plt.xscale('log')
plt.yscale('log')
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
plt.text(27, 2, fr'$f_{{col}} = {g_fit:.3f} \pm {g_fit_unc:.3f}$', fontsize=14)
plt.xlabel('Energy (keV)', fontsize=14)
ax.set_ylabel(r'$\mathrm{keV^2}\left[\mathrm{Photons\ s^{-1}\ keV^{-1}}\right]$', fontsize=14)
plt.legend(loc = 'upper right', fontsize=13)
plt.show()


# In[ ]:




