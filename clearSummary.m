savePath = fullfile(pwd,'Results');
now = datetime('now');
dateNow = yyyymmdd(datetime('now'));
timeNow = sprintf('%02d%02d%02d',hour(now),minute(now),floor(second(now)));

load(fullfile(savePath,'Summary.mat'));
archiveSaveName = [num2str(dateNow) '_' timeNow '_archiveSummary'];
save(fullfile(savePath,archiveSaveName),'Summary');

Summary = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
save(fullfile(savePath,'Summary'),'Summary');
clear Summary;