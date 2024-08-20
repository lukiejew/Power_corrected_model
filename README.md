# Power_corrected_model
#his Repository contains Python code used for the implementation of "Power-Law correction to a continuum-fitting emission model" by L.E. Jew.

The curve_modelling folder contains two files:

 model_curve.py: code for the implementation of the model approximation for the power-law segment curve, the blackbody curves, both with and without a fixed spectral hardening factor (Figure2).
 
 power_correction.py: code implementing the power-corrected model curve alongside the standard total model (Figure 5)

The Statistical_Analysis folder contains three files:

 cutoff.py: code to plot the modified 'total model' curve considering the cut-off point suggested by SSR analysis (Figure 4).

 KL.py: code for the calculation of the Kullback-Leibler divergence between the Power-corrected and original model curves.

 ssr_new.py: code for the calculation of SSR and the optimal cut-off Energy point, and the plotting of this value.
