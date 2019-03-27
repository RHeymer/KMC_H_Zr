%% siteID is a function that takes interstitial vector coords [a,b,c] and outputs the type of interstitial
% 1=O, 2=O', 3=T, 4=T', 5=T'', 6=T'''
% These are the same numbers as indices for ratesINN and vectorsINN

function siteType = siteID(x)
modCoords = mod(x,24);
if abs(modCoords-[16,16,6]) == 0 %approximately equal; accounts for float mismatch. Might be a better way.
    siteType = 1;
elseif abs(modCoords - [16,16,18])== 0
    siteType = 2;
elseif abs(modCoords - [0,0,9])== 0
    siteType = 3;
elseif abs(modCoords - [0,0,15])== 0
    siteType = 4;
elseif abs(modCoords - [8,8,3])== 0
    siteType = 5;
elseif abs(modCoords - [8,8,21])== 0
    siteType = 6;
else
    siteType = 0;
    error(['siteID did not find site for coords [',num2str(x),']']) %throw error if it doesn't match any of these.  
end


