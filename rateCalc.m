%% rateCalc simply calculates the rates for each movement before the main loop runs

rateOO = (23.32e12) *exp(-0.398*ec/kbT); %rates from Zhang et al.
rateOT = (23.32e12) *exp(-0.346*ec/kbT);
rateTO = (36.87e12) *exp(-0.406*ec/kbT);
rateTT = (36.87e12) *exp(-0.129*ec/kbT);
rateVAb = (1.61e13) *exp(-0.510*ec/kbT); %rates from Ruiz et al.
rateVAnb =(1.92e13) *exp(-0.530*ec/kbT);

%RatesINN is an array with rows corresponding to the current siteID, and
%columns as jump destinations, for interstitial H diffusion
ratesINN = [rateOO,rateOO,rateOT,rateOT,rateOT,rateOT,rateOT,rateOT;...
    rateOO,rateOO,rateOT,rateOT,rateOT,rateOT,rateOT,rateOT;...
    rateTT,rateTO,rateTO,rateTO,0,0,0,0;...
    rateTT,rateTO,rateTO,rateTO,0,0,0,0;...
    rateTT,rateTO,rateTO,rateTO,0,0,0,0;...
    rateTT,rateTO,rateTO,rateTO,0,0,0,0]; 

%RatesANN is similar for vacancies on the lattice.  Only one site type,
%with 12 destinations
ratesANN = [rateVAb rateVAb rateVAb rateVAb rateVAb rateVAb...
    rateVAnb rateVAnb rateVAnb rateVAnb rateVAnb rateVAnb]; 
sumRatesANN = sum(ratesANN);
