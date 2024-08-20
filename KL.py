#!/usr/bin/env python
# coding: utf-8

# In[3]:


#KL Div
import numpy as np
from scipy.integrate import quad
from scipy.special import rel_entr
import matplotlib.pyplot as plt

norm_bbfix = 18.6657
norm_bb = 18.5315
norm_pl = 2.70862

###########
def total_norm(x):
    d, a, c, b, m, f, g = 1.13501074e+00, -7.83291086e-03,  2.67747616e+10,  6.10924332e-01, 6.92494602e+04,  1.21742085,  1.74206201
    pl = d + (a * x - d) / ((1 + (x / c)**b)**m)
    exponent = f * x
    bb = 2*(g)**4 * x**3 / (np.exp(exponent) - 1)
    return (pl + bb)/21.1005

def separate_norm(x):
    a, c, d, m, b = -1.40065322e-02,  2.91701373e+10,  6.89663407e-01,  1.03757338e+05, 5.86783018e-01
    pl = d + (a * x - d) / ((1 + (x / c)**b)**m)
    g, f = 1.73623036, 1.21749901
    exponent = f * x
    bb = 2*(g)**4 * x**3 / (np.exp(exponent) - 1)
    return (pl + bb)/(20.986)

x_fit = np.linspace(3.58, 60, 1000)

# Calculate y values for the normalized models
y_fit_total_norm = total_norm(x_fit)
y_fit_separate_norm = separate_norm(x_fit)

# Calculate the KL divergence
kl_divergence = np.sum(rel_entr(y_fit_separate_norm, y_fit_total_norm))

print(f"KL Divergence: {kl_divergence}")

