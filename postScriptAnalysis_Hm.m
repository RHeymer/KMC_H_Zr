allFinalVacDists = NaN(numIterations,1);
allFinalVacDistsBasal = NaN(numIterations,1);
allFinalVacDistsNonBasal = NaN(numIterations,1);
allFinalHDists = NaN(numIterations,1);
allFinalHDistsBasal = NaN(numIterations,1);
allFinalHDistsNonBasal = NaN(numIterations,1);

% Calculate important distances
for iteration = 1:numIterations
    allFinalVacDists(iteration) = vec2dist(allFinalVacVectors(iteration,:)); %total distance
    allFinalHDists(iteration) = vec2dist(allFinalHVectors(iteration,:));
    
    allFinalVacDistsBasal(iteration) = vec2dist([allFinalVacVectors(iteration,1) 0 0]); %just distance on a plane
    allFinalHDistsBasal(iteration) = vec2dist([allFinalHVectors(iteration,1) 0 0]);
    
    allFinalVacDistsNonBasal(iteration) = vec2dist([0 0 allFinalVacVectors(iteration,3)]); %just distance on c plane
    allFinalHDistsNonBasal(iteration) = vec2dist([0 0 allFinalHVectors(iteration,3)]);
end

% Calculate overall diffusivities
allVacDiffusivities = (allFinalVacDists.^2)./(6*allFinalVacTimes);
meanVacDiffusivity = mean(allVacDiffusivities);
stdVacDiffusivity = std(allVacDiffusivities)/sqrt(numIterations); %I think this underestimates the variance; may need to think more

allHDiffusivities = (allFinalHDists.^2)./(6*allFinalHTimes);
meanHDiffusivity = mean(allHDiffusivities);
stdHDiffusivity = std(allHDiffusivities)/sqrt(numIterations);

% Calculate basal and non-basal diffusivities
allVacBasalDiffusivities = (allFinalVacDistsBasal.^2)./(6*allFinalVacTimes);
allVacNonBasalDiffusivities = (allFinalVacDistsNonBasal.^2)./(6*allFinalVacTimes);
meanVacBasalDiffusivity = mean(allVacBasalDiffusivities);
meanVacNonBasalDiffusivity = mean(allVacNonBasalDiffusivities);
stdVacBasalDiffusivity = std(allVacBasalDiffusivities)/sqrt(numIterations);
stdVacNonBasalDiffusivity = std(allVacNonBasalDiffusivities)/sqrt(numIterations);

allHBasalDiffusivities = (allFinalHDistsBasal.^2)./(6*allFinalHTimes);
allHNonBasalDiffusivities = (allFinalHDistsNonBasal.^2)./(6*allFinalHTimes);
meanHBasalDiffusivity = mean(allHBasalDiffusivities);
meanHNonBasalDiffusivity = mean(allHNonBasalDiffusivities);
stdHBasalDiffusivity = std(allHBasalDiffusivities)/sqrt(numIterations);
stdHNonBasalDiffusivity = std(allHNonBasalDiffusivities)/sqrt(numIterations);

% Calculate anisotropies
if (meanVacBasalDiffusivity ~= 0 || meanHBasalDiffusivity ~= 0)
    meanVacAnisotropy = meanVacNonBasalDiffusivity/meanVacBasalDiffusivity;
    stdVacAnisotropy = meanVacAnisotropy*sqrt((stdVacBasalDiffusivity/meanVacBasalDiffusivity)^2+(stdVacNonBasalDiffusivity/meanVacNonBasalDiffusivity)^2); %standard error propagation

    meanHAnisotropy = meanHNonBasalDiffusivity/meanHBasalDiffusivity;
    stdHAnisotropy = meanHAnisotropy*sqrt((stdHBasalDiffusivity/meanHBasalDiffusivity)^2+(stdHNonBasalDiffusivity/meanHNonBasalDiffusivity)^2);
else
    error('One of the mean diffusivities was zero, and you don`t want to divide by it (do more iterations next time)')
end

% Calculate diffusion speeds
allVacXVelocities = allFinalVacXDisplacements./allFinalVacTimes;
meanVacXVelocity = mean(allVacXVelocities);
stdVacXVelocity = std(allVacXVelocities)/sqrt(numIterations);

allHXVelocities = allFinalHXDisplacements./allFinalHTimes;
meanHXVelocity = mean(allHXVelocities);
stdHXVelocity = std(allHXVelocities)/sqrt(numIterations);

% Print results for testing
fprintf( 'Vacancy diffusivity = %0.6e +- %0.6e m^2 s^-1 with anisotropy %0.2f +- %0.2f and x velocity %0.2e +- %0.2e\n', meanVacDiffusivity, stdVacDiffusivity, meanVacAnisotropy, stdVacAnisotropy, meanVacXVelocity, stdVacXVelocity);
fprintf('Hydrogen diffusivity = %0.6e +- %0.6e m^2 s^-1 with anisotropy %0.2f +- %0.2f and x velocity %0.2e +- %0.2e\n', meanHDiffusivity, stdHDiffusivity, meanHAnisotropy, stdHAnisotropy, meanHXVelocity, stdHXVelocity);

%% Save the data (probably worth sticking in as a function once we're happy)

now = datetime('now');
dateNow = yyyymmdd(now);
timeNow = sprintf('%02d%02d%02d',hour(now),minute(now),floor(second(now)));
savePath = fullfile(pwd,'Results');


%Save the summary of the data in Summary.mat
saveDataSummary = [dateNow str2num(timeNow) T gradientHm numIterations numEvents...
    meanVacDiffusivity stdVacDiffusivity meanVacAnisotropy stdVacAnisotropy meanVacXVelocity stdVacXVelocity...
    meanHDiffusivity   stdHDiffusivity   meanHAnisotropy   stdHAnisotropy   meanHXVelocity   stdHXVelocity];
load(fullfile(savePath,'Summary.mat'));
Summary = [Summary; saveDataSummary];
save(fullfile(savePath,'Summary'),'Summary');
clear Summary;
