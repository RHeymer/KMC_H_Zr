%% rateCalc_T simply calculates the rates for each movement in a lattice 
% with a temperature gradient gradientT in the <a> direction
load('NNArrays.mat')
load('latticeCoords.mat')

%Calculate x-displacement for each jump
dxANN = a0/24*(vectorsANN(:,:,1)+(0.5*vectorsANN(:,:,2))); %dx = a0*(da+(0.5*db), where da&db scale from -1 to 1. We work on scale 24.
dxINN = a0/24*(vectorsINN(:,:,1)+(0.5*vectorsINN(:,:,2)));
%Calculate temperature difference dT across each jump
dTANN = gradientT * dxANN;
dTINN = gradientT * dxINN;

%Arrays below of jump frequencies nu and migration enthalpies Hm for each
%jump
nuANN = [1.61e13 1.61e13 1.61e13 1.61e13 1.61e13 1.61e13 1.92e13 1.92e13 1.92e13 1.92e13 1.92e13 1.92e13;...
    1.61e13 1.61e13 1.61e13 1.61e13 1.61e13 1.61e13 1.92e13 1.92e13 1.92e13 1.92e13 1.92e13 1.92e13];
nuINN = [23.32e12 23.32e12 23.32e12 23.32e12 23.32e12 23.32e12 23.32e12 23.32e12;...
    23.32e12 23.32e12 23.32e12 23.32e12 23.32e12 23.32e12 23.32e12 23.32e12;...
    36.87e12 36.87e12 36.87e12 36.87e12 0 0 0 0;...
    36.87e12 36.87e12 36.87e12 36.87e12 0 0 0 0;...
    36.87e12 36.87e12 36.87e12 36.87e12 0 0 0 0;...
    36.87e12 36.87e12 36.87e12 36.87e12 0 0 0 0];

HmANN = [0.51 0.51 0.51 0.51 0.51 0.51 0.53 0.53 0.53 0.53 0.53 0.53;...
    0.51 0.51 0.51 0.51 0.51 0.51 0.53 0.53 0.53 0.53 0.53 0.53];
HmINN = [0.398 0.398 0.346 0.346 0.346 0.346 0.346 0.346;...
    0.398 0.398 0.346 0.346 0.346 0.346 0.346 0.346;...
    0.129 0.406 0.406 0.406 99999 99999 99999 99999;...
    0.129 0.406 0.406 0.406 99999 99999 99999 99999;...
    0.129 0.406 0.406 0.406 99999 99999 99999 99999;...
    0.129 0.406 0.406 0.406 99999 99999 99999 99999]; %99999 to represent impossible jumps; unused

%Calculate the rates for each jump, accounting for temperature difference
%dT across the jump distance
ratesANN = nuANN.*exp((-HmANN.*ec)./(kb*(T+dTANN)));
ratesINN = nuINN.*exp((-HmINN.*ec)./(kb*(T+dTINN)));

sumRatesANN = sum(ratesANN(1,:)); %sum of rates from corners is same as sum of rates into corners; this is still always true    

%Calculate RateTT and RateTO from each tetrahedral site (siteID 3-6), for
%use in acceleration
allRateTTRateTO = NaN(6,2);
allRateTTRateTO(3:6,1) = ratesINN(3:6,1);
allRateTTRateTO(3:6,2) = mean(ratesINN(3:6,2:4),2);
%Calculate acceleration probabilities from each tetrahedral site 
%(siteID 3-6)
allPlPrTe = NaN(6,3);
allPlPrTe(:,1) = (allRateTTRateTO(:,1)+allRateTTRateTO(:,2))./((2*allRateTTRateTO(:,1))+allRateTTRateTO(:,2)); %Pl, where each row is current siteType
allPlPrTe(:,2) = allRateTTRateTO(:,1)./((2*allRateTTRateTO(:,1))+allRateTTRateTO(:,2)); %Pr, where each row corresponds to current siteType
allPlPrTe(:,3) = 1./allRateTTRateTO(:,2); %tExit