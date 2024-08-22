# Power_corrected_model
This Repository contains Python code used for the implementation of "Power-Law correction to a continuum-fitting emission model for Black Hole X-ray Binaries" by L.E. Jew.

The repository contains the following files:

  model_curve.py: code for the implementation of the model approximation for the power-law segment curve, the blackbody 
  curves, both with and without a fixed spectral hardening factor (Figure2).

  ssr_new.py: code for the calculation of SSR and the optimal cut-off Energy point, and the plotting of this value. (Figure 
  3)

  cutoff.py: code to plot the modified 'total model' curve considering the cut-off point suggested by SSR analysis (Figure 
  4).

  cutoff_KL.py: code for the calculation and plotting of K-L divergence against cut-off Energy (Figure 5).

  comp_kl.py: code for the plotting of the modified total curve considering the K-L-optimized cut-off Energy, alongside 
  the residuals in a subplot (Figure 6).

  power_correction.py: code implementing the power-corrected model curve alongside the standard total model (Figure 7), 
  with a subplot of the residuals.
  
  KL.py: code for the calculation of the Kullback-Leibler divergence between the Power-corrected and original model curves.
  
  

  
