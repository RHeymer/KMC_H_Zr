function chosenVectorID = chooseVector(rates)

    numRates = size(rates,2); %gets the length of the longest dimension of rates
    sumRates = sum(rates);
    cumuRates = NaN(numRates,1);
    chosenVectorID = 0;
  
    % Generate cumulative rates
    cumuRates(1) = rates(1);
    for i = 2:numRates
        cumuRates(i) = cumuRates(i-1) + rates(i);
    end
    cumuRates = cumuRates./sumRates; % Normalise to 1
    randNum = rand;
    truthArray = zeros(numRates,1);
    % For each value in cumuRates, compare it to the random number; store a
    % 1 if it is larger
    for j = 1:numRates
            if randNum < cumuRates(j)
                truthArray(j) = 1;
            end
    end
    %The chosen vector is the lowest cumulative rate where it is greater than the random number    
    chosenVectorID = find(truthArray,1);
    
    if chosenVectorID == 0
        error('Vector Was Not Chosen; check rates')
    end
    
    