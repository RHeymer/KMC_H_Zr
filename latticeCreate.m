%% Create the lattice (if we're just using vectors, we might not need this?)
a0 = 3.232e-10; %lattice parameter in metres
scSize = 3; %number of unit cells in the repeated supercell

%% atomCoords is a 54 (atom sites) * 3 (coords a,b,c) array
% with vector coords of atoms in lattice
% odd indices refer to [0,0,0] initial atom repeats
% even indices refer to [8,8,12] initial atom repeats

atomCoords = zeros(54,3);
atomInitCoords = [0 0 0;8 8 12];
i = 1; %counting index for dimension 1 of atomCoords

for xRep = 0:2 
    for yRep = 0:2
        for zRep = 0:2
            for atomInit = 1:2
                atomInitCoord = atomInitCoords(atomInit,:);
                atomCoords(i,:) = atomInitCoord + 24.*[xRep yRep zRep];
                i = i+1;
            end
        end
    end
end


%% inCoords is a 162 (total interstitial sites (6 types*27cells)) * 3 (coords a,b,c) array

inCoords = zeros(162,3);
inInitCoords = [16 16 6; 16 16 18; 0,0,9; 0,0,15; 8,8,3; 8 8 21];
i = 1; %counting index for dimension 1 of inCoords

for xRep = 0:2 
    for yRep = 0:2
        for zRep = 0:2
            for inInit = 1:6
                inInitCoord = inInitCoords(inInit,:);
                inCoords(i,:) = inInitCoord + 24*[xRep yRep zRep];
                i = i+1;
            end
        end
    end
end

% if working on coords 0 to 1 rather than 0 to 3, uncomment:
% inCoords = inCoords / 3;

save('latticeCoords','atomCoords','inCoords')

%% Define initial vacancy and hydrogen positions

initVPos = atomCoords(randi(54),:);
initHPos = inCoords(randi(162),:);

    

