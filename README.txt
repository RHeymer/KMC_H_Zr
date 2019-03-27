chooseVector chooses a vector with probability proportional to its rate, given an array of rates.
clearSummary clears Summary.mat in the Results folder, and saves it as an archive.
clearSummary_gradient does the same, but with extra columns for gradient related data.
KMC runs kinetic monte carlo code for vacancies and H in alpha zirconium.
KMC_Hm does the same, but for a given migration enthalpy gradient.
KMC_T does the same, but for a given temperature gradient.
latticeCoords contains arrays of coordinates of lattice sites and interstitial sites.  Multiplied by 24 for integer usage instead of decimals.
latticeCreate creates the lattices, making latticeCoords
NN creates the vectors between nearest neighbours in the lattice (for on-lattice and interstitials).
NNArrays contains these vector arrays.
postScriptAnalysis analyses bulk data from running KMC, calculates useful values and saves them in Results/Summary.mat
postScriptAnalysis_Hm does the same for KMC_Hm
postScriptAnalysis_T does the same for KMC_T.
rateCalc creates an array of rates of jumps between sites, used in KMC.
rateCalc_Hm does the same for KMC_Hm
rateCalc_T does the smae for KMC_T
siteID identifies the type of interstitial, used in all KMC files.
vec2dist converts a vector between sites (scaled to 1 per unit cell, not 24) into a distance.