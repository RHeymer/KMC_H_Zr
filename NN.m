a0 = 3.232e-10; %lattice parameter in metres
%% vectorsINN is a 6 (current interstitial site)* 8 (nearest neighbour interstitial site) * 3 (directions a,b,c) array
% with values of the vectors [da, db, dc] to get from current to NN site.
% Current site is in order of O, O', T, T', T'', T''
% For current site O,    NN in order of O',   O', T,  T,  T,  T', T', T'
% For current site O',   NN in order of O,    O,  T', T', T', T,  T,  T
% For current site T,    NN in order of T',   O,  O,  O
% For current site T',   NN in order of T,    O', O', O'
% For current site T'',  NN in order of T''', O,  O,  O
% For current site T''', NN in order of T'',  O', O', O'
% Empty spaces are filled with 0.

vectorsINN = zeros(6,8,3);
vectorsINN(1,:,:) = [0 0 +12; 0 0 -12; +8 +8 +3; +8 -16 +3; -16 +8 +3; -8 -8 -3; -8 +16 -3; +16 -8 -3];
vectorsINN(2,:,:) = [0 0 -12; 0 0 +12; +8 +8 -3; +8 -16 -3; -16 +8 -3; -8 -8 +3; -8 +16 +3; +16 -8 +3];

vectorsINN(3,1:4,:) = [0 0 +6; -8 -8 -3; -8 +16 -3; +16 -8 -3];
vectorsINN(4,1:4,:) = [0 0 -6; -8 -8 +3; -8 +16 +3; +16 -8 +3];
vectorsINN(5,1:4,:) = [0 0 -6; +8 +8 +3; +8 -16 +3; -16 +8 +3];
vectorsINN(6,1:4,:) = [0 0 +6; +8 +8 -3; +8 -16 -3; -16 +8 -3];



%% vectorsANN is a 2(middle/corner unit cell site) * 12 (nearest neighbour atom site) * 3 (directions a,b,c) array
% with values of the vectors [da, db, dc] to get from current to NN lattice site.
% All sites are energetically identical; rate of movement between them is
% based on rateVA only.

vectorsANN = zeros(2,12,3);
vectorsANN(1,:,:) = [24 0 0; 24 -24 0; 0 -24 0; -24 0 0; -24 24 0; 0 24 0;...
    8 8 12; 8 -16 12; -16 8 12; 8 8 -12; 8 -16 -12; -16 8 -12];
vectorsANN(2,:,:) = -vectorsANN(1,:,:);

save('NNArrays','vectorsINN','vectorsANN')





