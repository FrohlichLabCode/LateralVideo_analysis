% Select all visually responsive channels across sessions

%% prepare diSessiontory
clear all

cluster = 0;
skipRec = 0;
animalCode = '0172';
analysisType = 'Saccade';
folderSuffix = '_validChns_new';
doPlot = 1;
lfpFs = 1000;


if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/' analysisType '/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];  
end

fileInfo   = dir([AnalysisDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto

% animal info
regionNames = {'lEye','rEye'};
condNames   = {'lVideo','rVideo'};
numRegions  = numel(regionNames);
numConds    = numel(condNames);
%fullVideoSessions = {'039','042','045','048','051','054','058','062','064','069','073','077','083','086','092','100','136'};
%fullVideoSessions = {'039','042','045','048','051','054','058','061','062','064','067','068','069','073','074',...
%                     '077','078','083','084','086','087','091','092','097','100'};
% 
% %if find(contains(fullVideoSessions,sessionName))>0; fullVideo = 1; else fullVideo = 0; end
% session2analyze = [1:100];
% % -----------
% monocularIndex = [2,5,10,15,16,21,23,24,26,27,76,82,90,93,96]; %only analyze some sessions
% oneEyeSessions = [1:28]; %not code to process this yet
% twoEyeSessions = [29:100];
% %excludeSessions= [34,
% session2analyze = twoEyeSessions;
% numSessions     = numel(session2analyze);
% -----------
medial  = [40,42,43,45,46,48,49,51,52,54,55,57,59,60,62,63,65,66,68,70];
lateral = [1,2,4,6,8,10,12,14,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39];
half    = [3,5,7,9,11,13,15,17,20,23,26,29,32,35,38,41,44,47,53,56,58,61,64,69]; %50

sessionTypes =  {medial, lateral, half};
sessionTypeNames = {'medial', 'lateral', 'half'};



% parameteres
twin = [-2 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-2,8];
binSize = 0.5;    % in seconds
baseTwin = [-2 0]; %<<<--- baseline to subtract from spectrogram 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Saccade PSTH %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop through each recording
if ~exist(join(GroupAnalysisDir),'dir'); mkdir(join(GroupAnalysisDir));end % to save figures for each session

% allSaccades = cell(numSessions, numRegions); % have to pre-allocate
% allEvents   = cell(numSessions, numConds);
for i = 1:numel(sessionTypes)
    sessionIDs = sessionTypes{i};
    sessionType = sessionTypeNames{i};
    irec = 0;
    sacPSTH{numRegions,numConds} = [];
    sacPSTHTrial{numRegions,numConds} = [];
    sacPSTHNormed{numRegions,numConds} = [];
    sessionNames = {};
for iSession = 1:numel(sessionIDs)
    sessionID = sessionIDs(iSession);
    sessionName = sprintf('%03d',sessionID); % add leading zeros        
    fileInfo = dir([AnalysisDir animalCode '_LateralVideo_' sessionName '*']); % detect files to load/convert  '_LateralVideo*'
    irec = irec +1;
    recName = fileInfo.name;
    rootPreprocessDir     = [PreprocessDir recName '/'];
    
%     if ismember(str2num(sessionName),session2analyze)
%     if find(contains(fullVideoSessions,sessionName))>0; continue;end % skip fullVideoSessions
%     sessionNames{end+1} = {sessionName};    
    if exist(join([rootPreprocessDir, 'saccades.mat']),'file')
       sacTimes = is_load([rootPreprocessDir ,'saccades'],'sacTime'); % directly load saved pupil data
    else
       [sacTimes,~,~,~] = is_detectSaccadesEphysFun_V3(rootPreprocessDir);
    end
    

% sacSamp = is_load([rootPreprocessDir 'saccades.mat'], 'sacSamp'); %numCond x numChn x tvecPupil
    [evtTimes,lfpFs] = is_load([rootPreprocessDir 'eventTimes.mat'],'evtTimes','lfpFs');
    condID = [1,2];
%     allSaccades(iSession, :) = flip(sacTime,2); % get Left and Right order
%     allEvents(iSession, :) = is_load([rootPreprocessDir 'eventTimes'],'evtTimes');
    for iEye  = 1:numRegions
        for iCond = 1:numConds
            sacTime = sacTimes{iEye};
            evtTime = evtTimes{condID(iCond)};
            try % some sessions don't have saccades detected
            [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,sacTime,twin,binSize); % CZ: PSTHrate's 1st timept is time of saccade
            catch
            end
            %data2Average = psthTrial(any(psthTrial,2),:); %get number of trials with saccade data
            data2Average = psthTrial;
            sacPSTH{iEye,iCond} = [sacPSTH{iEye,iCond}; sum(psthTrial,1)/size(data2Average,1)/binSize];
            sacPSTHTrial{iEye,iCond} = [sacPSTHTrial{iEye,iCond}; psthTrial];
            [temp_normed,~,~] = baselineNorm(nansum(psthTrial,1)/size(data2Average,1)/binSize, timePSTH, baseTwin, 1);
            sacPSTHNormed{iEye,iCond} = [sacPSTHNormed{iEye,iCond}; temp_normed];
        end
    end
    end




save([GroupAnalysisDir 'SaccadeAll_' sessionType '.mat'],'timePSTH','sacPSTH','sacPSTHTrial','sacPSTHNormed','regionNames', 'condNames','sessionNames', '-v7.3');

%%%%%%%%%%%%%%%%%%
%% Average Pupil %%
%%%%%%%%%%%%%%%%%%

% median
fig = figure('name','Session mean saccade PSTH');
iPanel = 1;
for iCond = 1:numConds
    for iRegion = 1:numRegions
        subplot(numConds,numRegions,iPanel)
        iPanel = iPanel + 1;
        mean = nanmean(sacPSTH{iRegion,iCond},1);
        median = nanmedian(sacPSTH{iRegion,iCond},1); 
        sem = nanstd(sacPSTH{iRegion,iCond},[],1)/sqrt(size(sacPSTHTrial{iRegion,iCond},1)); 
        shadedErrorBar(timePSTH,mean,sem); %ylim([0,1.6])%,@std});         
        title(regionNames{iRegion});ylabel([condNames{condID(iCond)} ': [saccades/sec]']);
        xlim(xLim);
    end %median moves baseline down, so use mean!
end
savefig(fig, [GroupAnalysisDir 'Session_mean_saccadePSTH_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'Session_mean_saccadePSTH_' sessionType '.png']);
 

% normalized
fig = figure('name','baselineNormed session mean saccade PSTH');
iPanel = 1;
for iCond = 1:numConds
    for iRegion = 1:numRegions
        subplot(numConds,numRegions,iPanel)
        iPanel = iPanel + 1;
        mean = nanmean(sacPSTHNormed{iRegion,iCond},1);
        median = nanmedian(sacPSTHNormed{iRegion,iCond},1); 
        sem = nanstd(sacPSTHNormed{iRegion,iCond},[],1)/sqrt(size(sacPSTHNormed{iRegion,iCond},1)); 
        shadedErrorBar(timePSTH,mean,sem); %,@std});        
        title(regionNames{iRegion});ylabel(condNames{condID(iCond)});
        xlim(xLim);
    end %median moves baseline down, so use mean!
end
        
savefig(fig, [GroupAnalysisDir 'Session_mean_saccadePSTH_baselineNormed_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'Session_mean_saccadePSTH_baselineNormed_' sessionType '.png']);


end