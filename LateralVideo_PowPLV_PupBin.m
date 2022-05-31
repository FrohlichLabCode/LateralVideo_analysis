clear
tic


cluster = 0;
animalCode = '0169';
%excludeWin = [-0.1 0.6]; % in seconds

if cluster == 0
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/'];

elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];    %code for initialising parallel computing
    
%     c = parcluster('local'); % build the 'local' cluster object
%     numCore = c.NumWorkers;        % get the number of workers
%     saveProfile(c);  
% 
%     display(num2str(numCore))
%     
%     parpool('local',numCore,'SpmdEnabled',false);
end

fileInfo   = dir([AnalysisDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto

theseRecs = [1:100]; %29:100 are 2eyes sessions
numRecs = numel(theseRecs);

for irec =  31%:numRecs%1:numel(fileInfo) %12:15%
    
    rec2Analyze = theseRecs(irec);
    
    recName = fileInfo(rec2Analyze).name;
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/PupilBin/'];
    
    %is_detectSaccadesEphys_V3(rootPreprocessDir)    
    %is_pupilState(rootPreprocessDir)
    
    % find epochs to exclude
%     % if use trigger data:
%     [triggerData,Fs] = is_load([rootPreprocessDir '/triggerData'],'triggerData','Fs');
%     % load trigger data
%     ttlInd = 1;
%     excludeTimes = find(diff(triggerData(ttlInd,:))==1)./Fs;
     
    % better use processed event time data: (took photodiode data when possible)
    [evtTimes,condNames,condID] = is_load([rootPreprocessDir 'eventTimes.mat'], 'evtTimes', 'condNames','condID');
%     for iCond = 1:numel(condID)
%         excludeTimes = evtTimes{iCond};
%         condName = condNames{condID(iCond)};
%     end
    for iAnalysis = 1:3
        if iAnalysis == 1
            excludeTimes = [0.01]; excludeWin = [0 10]; % exclude first 10 sec of recording,can't use 0, error in indexing
            condName = 'all';
        elseif iAnalysis == 2
            excludeTimes = evtTimes{2}; %only leave left video trials [-1,7]
            excludeWin   = [-3 8];
            condName = 'lVideo';
        elseif iAnalysis == 3
            excludeTimes = evtTimes{1}; %only leave right video trials
            excludeWin   = [-3 8];
            condName = 'rVideo';
        end
    is_pupilReanalysis_2eyes(rootPreprocessDir,rootAnalysisDir,condName,excludeTimes,excludeWin)
    end
end