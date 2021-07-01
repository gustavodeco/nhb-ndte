# nhb-ndte
MATLAB code for Deco, Vidaurre and Kringelbach paper in Nature Human Behaviour:
Deco, G., Vidaurre, D. & Kringelbach, M.L. Revisiting the global workspace orchestrating the hierarchical organization of the human brain. Nat Hum Behav 5, 497–511 (2021). https://doi.org/10.1038/s41562-020-01003-6

In the sub-directory NDTE, you will find the MATLAB code for computing the NDTE for the HCP data (rest and seven tasks). The data is freely available from https://www.humanconnectome.org/ in the young adults HCP cohort and has been preprocessed as specified in the paper (i.e. pre-processed timeseries for the DK80 parcellation). 
This produces the results used for Figure 2 in the published paper. More specifically, running e.g. ndte_all(1) will produce the results for the resting state data.  
In the MODEL sub-directory, you will find MATLAB code for the whole-brain model which produces the data for Figure 5. 
The main program is fitt_hopf_a_particleswarm.m, which will generate the results for the whole-brain model.
