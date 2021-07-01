# nhb-ndte
Code for Deco, Vidaurre and Kringelbach paper in Nature Human Behaviour

Deco, G., Vidaurre, D. & Kringelbach, M.L. Revisiting the global workspace orchestrating the hierarchical organization of the human brain. Nat Hum Behav 5, 497–511 (2021). https://doi.org/10.1038/s41562-020-01003-6

Within this directory, we have created the sub-directory NDTE, where we have uploaded the MATLAB code for computing the NDTE for the HCP data (rest and seven tasks). The data is freely available from https://www.humanconnectome.org/ in the young adults HCP cohort and has been preprocessed as specified in the paper (i.e. pre-processed timeseries for the DK80 parcellation). 

This produces the results used for Figure 2 in the published paper. More specifically, running e.g. ndte_all(1) will produce the results for the resting state data.  
In addition, in the MODEL sub-directory, we have uploaded the whole-brain model which produces the data for Figure 5. 
