%% Main loop of code; runs Kinetic Monte Carlo Simulation for vacancies and hydorgen diffusion in alpha-Zr at a given temperature

%% Constants
kb = 1.38064852e-23; %Boltzmann constant
ec = 1.60217662e-19; %elementary charge
a0 = 3.232e-10; %lattice parameter in metres

%% Inputs
for Tc = [0,50,100,150,200,210,220,230,240,250,260,270,280,290] %in degrees Celsius
    T = Tc + 273.15; %Kelvin
    kbT = kb*T; %used frequently in rate calculations
    numEvents = 100000; %number of steps to run each KMC for
    numIterations = 10000; %number of times to run KMC at this temperature
    
    %% Load Arrays
    load('latticeCoords') %load coords of atoms and interstitials in the 3*3*3 supercell
    rateCalc %calculates rates
    
    %% Probability calculations (For TT trapping)
    Pr = double(rateTT)/(2*double(rateTT)+double(rateTO));
    Pl = double((rateTT+rateTO))/double((2*rateTT+rateTO));
    tExit = 1/rateTO;
    accelerationCompleted = false;
    
    %% Empty Arrays
    allFinalVacVectors = NaN(numIterations,3);
    allFinalVacTimes = NaN(numIterations,1);
    allFinalHVectors = NaN(numIterations,3);
    allFinalHTimes = NaN(numIterations,1);
    
    %% Main loop
    for iteration = 1:numIterations
        
        initVacPos = atomCoords(randi(54),:); %random start positions within the supercell
        initHPos = inCoords(randi(162),:);
        
        % FOR VACANCY DIFFUSION
        
        totalRatesTimesVac = zeros(numEvents+1,2);%first column rates (unused, included for future study), 2nd time
        storePosVac = NaN(numEvents+1,3);
        storePosVac(1,:) = initVacPos;
        
        for i=2:numEvents+1
            deltaTimeVA = -1*(log(rand(1,1))/sumRatesANN); %sumRatesANN is sum of all rates, see rateCalc.m
            totalRatesTimesVac(i,2) = totalRatesTimesVac((i-1),2)+deltaTimeVA; %increment time
            
            siteRatesVac = ratesANN; % vector containing all rates from any vacancy pos (could put outside loop for now, but if there's anything site dependent later we may want it here
            chosenVectorIDVac = chooseVector(siteRatesVac); % gives vector ID, based on random selection proportional to the rate
            
            chosenVectorVac = vectorsANN(chosenVectorIDVac,:); %gives corresponding vector (from corner site)
            if all(abs(mod(storePosVac(i-1,:),24))) ~= 0 %If the vacancy is at the "middle" site in the unit cell (not the corners which all have coords multiples of 24)
                chosenVectorVac = -1.*chosenVectorVac; %Make the vector negative (-1 * vectors from corner sites = 1* vecs to corner sites = vecs from middle sites)
            end
            % The above method is the most efficient way of managing this
            % with identical rates in all directions
            
            storePosVac(i,:) = storePosVac(i-1,:) + chosenVectorVac; %Increment position
            deltaTimeVac = -1*log(rand(1))./sum(siteRatesVac);            
        end
        
        % FOR HYDROGEN DIFFUSION
        totalRatesTimesINN = zeros(numEvents+1,2); %first row is rates, second is times
        storePosINN = NaN(numEvents+1,3);
        storePosINN(1,:) = initHPos;
        
        for i=2:numEvents+1
            
            siteType=siteID(storePosINN(i-1,:));
            if (siteType == 1 || siteType == 2) % If currently occupying octahedral site
                numPossMoves = 8;
                siteRatesINN = ratesINN(siteType,1:numPossMoves); % vector containing all rates from current site
                chosenVectorID = chooseVector(siteRatesINN); % gives vector ID, based on random selection proportional to the rate
                chosenVectorINN = vectorsINN(siteType,chosenVectorID,:);
                chosenVectorINNa = chosenVectorINN(1);
                chosenVectorINNb = chosenVectorINN(2);
                chosenVectorINNc = chosenVectorINN(3);
                storePosINN(i,:) = storePosINN(i-1,:) + [chosenVectorINNa, chosenVectorINNb, chosenVectorINNc]; % increment position
                deltaTimeINN = -1*log(rand(1))./sum(siteRatesINN);
                totalRatesTimesINN(i,2) = totalRatesTimesINN(i-1,2) + deltaTimeINN;
                
            else % If currently occupying tetrahedral site
                numPossMoves = 4;
                siteRatesINN = ratesINN(siteType,1:numPossMoves); % vector containing all rates from current site
                chosenVectorID = chooseVector(siteRatesINN); % gives vector ID, based on random selection proportional to the rate
                chosenVectorINN = vectorsINN(siteType,chosenVectorID,:);
                
                if chosenVectorID == 1 %TT movement
                    %If the H atom is in a tetrahedral site and wants to move
                    %to the next tetrahedral site, apply acceleration
                    accelerationCompleted = true;
                    leftOrRightSelection = rand;
                    
                    %PL > PR should always hold therefore
                    if leftOrRightSelection <= Pr
                        %exits right, so from other T site
                        %hydrogen moves to other T site and then exits
                        chosenVectorINNa = chosenVectorINN(1);
                        chosenVectorINNb = chosenVectorINN(2);
                        chosenVectorINNc = chosenVectorINN(3);
                        accelerationPosINN = storePosINN(i-1,:) + [chosenVectorINNa, chosenVectorINNb, chosenVectorINNc];
                        siteType = siteID(accelerationPosINN);
                        accelerationSiteRatesINN = ratesINN(siteType,(2:numPossMoves));
                        accelerationChosenVectorID = chooseVector(accelerationSiteRatesINN);
                        accelerationChosenVectorINN = vectorsINN(siteType,accelerationChosenVectorID,:);
                        chosenVectorINN = chosenVectorINN + accelerationChosenVectorINN;
                    else
                        %exits left, so from where it enters the basin
                        %hydrogen exits from current T site
                        accelerationPosINN = storePosINN(i-1,:);
                        siteType = siteID(accelerationPosINN);
                        accelerationSiteRatesINN = ratesINN(siteType,(2:numPossMoves));
                        accelerationChosenVectorID = chooseVector(accelerationSiteRatesINN);
                        accelerationChosenVectorINN = vectorsINN(siteType,accelerationChosenVectorID,:);
                        chosenVectorINN = accelerationChosenVectorINN;
                    end
                end
                
                chosenVectorINNa = chosenVectorINN(1);
                chosenVectorINNb = chosenVectorINN(2);
                chosenVectorINNc = chosenVectorINN(3);
                storePosINN(i,:) = storePosINN(i-1,:) + [chosenVectorINNa, chosenVectorINNb, chosenVectorINNc];
                
                if accelerationCompleted == false % if we don't accelerate
                    deltaTimeINN = -1*log(rand(1))./sum(siteRatesINN);
                    totalRatesTimesINN(i,2) = totalRatesTimesINN(i-1,2) + deltaTimeINN; %increment everything normally
                else
                    deltaTimeINN = -1*tExit*log(rand(1));
                    totalRatesTimesINN(i,1) = rateTT;
                    totalRatesTimesINN(i,2) = totalRatesTimesINN(i-1,2) + deltaTimeINN;
                    accelerationCompleted = false;
                end
            end
        end
        
        %After all events have occurred per run
        totalVectorVac = (storePosVac(end,:) - storePosVac(1,:)) / 24; %Divide by 24 here to get true vector based on 1*1*1 unit cell
        totalTimeVac = totalRatesTimesVac(end,2);
        allFinalVacVectors(iteration,:) = totalVectorVac;
        allFinalVacTimes(iteration,:) = totalTimeVac;
        
        totalVectorH = (storePosINN(end,:) - storePosINN(1,:)) / 24;
        totalTimeH = totalRatesTimesINN(end,2);
        allFinalHVectors(iteration,:) = totalVectorH;
        allFinalHTimes(iteration,:) = totalTimeH;
    end
    postScriptAnalysis % See separate file
end