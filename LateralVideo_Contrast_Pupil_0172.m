% Select all visually responsive channels across sessions

%% prepare directory
clear all

cluster = 0;
skipRec = 1;
animalCode = '0172';
analysisType = 'Pupil';
folderSuffix = '_validChns';
doPlot = 1;



if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/Pupil/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];  
end

fileInfo   = dir([AnalysisDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto

% animal info
regionNames = {'lEye','rEye'};
condNames   = {'lVideo', 'rVideo'};
numRegion   = numel(regionNames);
numConds    = numel(condNames);
%fullVideoSessions = {'039','042','045','048','051','054','058','062','064','069','073','077','083','086','092','100','136'};
%fullVideoSessions = {'039','042','045','048','051','054','058','061','062','064','067','068','069','073','074',...
%                      '077','078','083','084','086','087','091','092','097','100'};

% %if find(contains(fullVideoSessions,sessionID))>0; fullVideo = 1; else fullVideo = 0; end
% session2analyze = [1:100];
% % -----------
% monocularIndex = [2,5,10,15,16,21,23,24,26,27,76,82,90,93,96]; %only analyze some sessions
% oneEyeSessions = [1:28]; %not code to process this yet
% twoEyeSessions = [29:100];
% session2analyze = twoEyeSessions;
% -----------
medial  = [40,42,43,45,46,48,49,51,52,54,55,57,59,60,62,63,65,66,68,70];
lateral = [1,2,4,6,8,10,12,14,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39];
half    = [3,5,7,9,11,13,15,17,20,23,26,29,32,35,38,41,44,47,53,56,58,61,64,69]; %50

sessionTypes =  {medial, lateral, half};
sessionTypeNames = {'medial', 'lateral', 'half'};

% parameteres
twin = [-2 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-2,8];
baseTwin = [-2 0]; %<<<--- baseline to subtract from spectrogram 
tvecPupil = is_load('D:\FerretData\0172\Analyzed\0172_LateralVideo_001_20181218\Pupil\Pupil_baselineNormed.mat','timeVec');
eyeChn = {'leftX','rightX','leftY','rightY','leftD','rightD'};

%%%%%%%%%%%%%%
%% Contrast %%
%%%%%%%%%%%%%%

% loop through each recording
if ~exist(join([GroupAnalysisDir 'contrastFigures/']),'dir'); mkdir(join([GroupAnalysisDir 'contrastFigures/']));end % to save figures for each session


for i = 1:numel(sessionTypes)
    sessionIDs = sessionTypes{i};
    sessionType = sessionTypeNames{i};
    irec = 0;
    allPupil = cell(numConds, numRegion); % have to pre-allocate

    for iSession = 1:numel(sessionIDs)
        sessionID = sessionIDs(iSession);
        sessionName = sprintf('%03d',sessionID); % add leading zeros        
        fileInfo = dir([AnalysisDir animalCode '_LateralVideo_' sessionName '*']); % detect files to load/convert  '_LateralVideo*'
        irec = irec +1;
        recName = fileInfo.name;
    %if ismember(str2num(sessionName),session2analyze)
    %if find(contains(fullVideoSessions,sessionName))>0; continue;end % skip fullVideoSessions
    rootPupilDir     = [AnalysisDir recName '/Pupil/'];
    if exist([rootPupilDir 'Pupil_baselineNormed.mat'])
        pupil = is_load([rootPupilDir 'Pupil_baselineNormed.mat'], 'evtDat_trialAvg'); %numCond x numChn x tvecPupil
    else
        pupil = is_load([rootPupilDir 'Pupil_baselineNormed_median.mat'], 'evtDat_trialAvg'); %numCond x numChn x tvecPupil
    end
    display(['processing ' rootPupilDir])
    for iCond = 1:numConds
        for iRegion = 1:numRegion
            sliceData = reshape(pupil(iCond,4+iRegion,:),1,[]); %numChn x tvecPupil
            data2Average = sliceData; % delete channels (2nd dimension) without spikes
            allPupil{iCond,iRegion} = [allPupil{iCond,iRegion}; data2Average];  % append the visual channel
        end
    end
    end

             
save([GroupAnalysisDir 'PupilAll_' sessionType '.mat'],'tvecPupil','allPupil', '-v7.3');

%%%%%%%%%%%%%%%%%%
%% Average Pupil %%
%%%%%%%%%%%%%%%%%%

% median
for iCond = 1:numConds
    for iRegion = 1:numRegion
        avgPupil(iCond,iRegion,:) = nanmedian(allPupil{iCond,iRegion},1); % 1st dimension is channel
    end 
end

Pupil.contraAvg(1,:) = squeeze(avgPupil(2,1,:)); %left eye
Pupil.contraAvg(2,:) = squeeze(avgPupil(1,2,:)); %right eye
Pupil.ipsiAvg(1,:)   = squeeze(avgPupil(1,1,:));
Pupil.ipsiAvg(2,:)   = squeeze(avgPupil(2,2,:));
Pupil.contrast{1}    = allPupil{1,1} - allPupil{2,1}; %leye=lvideo-rvideo
Pupil.contrast{2}    = allPupil{2,2} - allPupil{1,2}; %reye=rvideo-lvideo
Pupil.contrastAvg    = Pupil.ipsiAvg - Pupil.contraAvg;
Pupil.leftAvg        = squeeze(avgPupil(1,:,:));
Pupil.rightAvg       = squeeze(avgPupil(2,:,:));
Pupil.l_r{1}         = allPupil{1,1} - allPupil{2,1}; % left-right
Pupil.l_r{2}         = allPupil{1,2} - allPupil{2,2}; % left-right
Pupil.l_rAvg         = Pupil.leftAvg - Pupil.rightAvg;
save([GroupAnalysisDir 'PupilAvg_median.mat'],'Pupil','tvecPupil','avgPupil', '-v7.3');


% plot Pupil
%-------- plot Pupil for all channels
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/3 (screensize(4)-150)*2/3]); %(x,y,width,height) screensize(3)-100

for iRegion = 1:numRegion
    if iRegion == 1;yLim = [-1.5,1.5];
    elseif iRegion == 2; yLim = [-1.5,1.5];end
    
    subplot(3,numRegion,iRegion)
    mean = squeeze(Pupil.leftAvg(iRegion,:));
    sem  = nanstd(allPupil{1,iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.leftAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);
    title(regionNames{iRegion});
    if iRegion == 1;ylabel('Left normed amplitude');end
    
    subplot(3,numRegion,iRegion+numRegion)
    mean = squeeze(Pupil.rightAvg(iRegion,:));
    sem  = nanstd(allPupil{2,iRegion},[],1)/sqrt(size(allPupil{2,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.rightAvg(iRegion,:)),'k',  'LineWidth', 1);
    xlim(xLim); ylim(yLim);%ylim([-10,5]);
    if iRegion == 1;ylabel('Right normed amplitude');end
    
    
    
    subplot(3,numRegion,iRegion+2*numRegion)
    mean = squeeze(nanmedian(Pupil.contrast{iRegion},1));
    sem  = nanstd(Pupil.contrast{iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.l_rAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);%ylim([-1,2]); %ylim([-0.8,0.8]);
    xlabel('Time [s]');
    if iRegion == 1; ylabel('Ipsi-Contra normed amplitude');end
    
%     subplot(3,numRegion,iRegion+2*numRegion)
%     mean = squeeze(nanmedian(Pupil.l_r{iRegion},1));
%     sem  = nanstd(Pupil.l_r{iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
%     shadedErrorBar(tvecPupil,mean,sem);
%     %plot(tvecPupil, squeeze(Pupil.l_rAvg(iRegion,:)),'k', 'LineWidth', 1);
%     xlim(xLim); ylim(yLim);%ylim([-1,2]); %ylim([-0.8,0.8]);
%     xlabel('Time [s]');
%     if iRegion == 1; ylabel('Left-right normed amplitude');end
end

savefig(fig, [GroupAnalysisDir 'contrastAvg_Pupil_median_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_Pupil_median_' sessionType '.png']);

%% channel mean (vs. median)
for iCond = 1:numConds
    for iRegion = 1:numRegion
        avgPupil(iCond,iRegion,:) = nanmean(allPupil{iCond,iRegion},1); % 1st dimension is channel
    end %median moves baseline down, so use mean!
end

% session average
Pupil.contraAvg(1,:) = squeeze(avgPupil(2,1,:)); %left eye
Pupil.contraAvg(2,:) = squeeze(avgPupil(1,2,:)); %right eye
Pupil.ipsiAvg(1,:)   = squeeze(avgPupil(1,1,:));
Pupil.ipsiAvg(2,:)   = squeeze(avgPupil(2,2,:));
Pupil.contrastAvg    = Pupil.contraAvg - Pupil.ipsiAvg;
Pupil.leftAvg        = squeeze(avgPupil(1,:,:));
Pupil.rightAvg       = squeeze(avgPupil(2,:,:));
Pupil.l_r{1}         = allPupil{1,1} - allPupil{2,1}; % left-right
Pupil.l_r{2}         = allPupil{1,2} - allPupil{2,2}; % left-right
Pupil.l_rAvg         = Pupil.leftAvg - Pupil.rightAvg;
save([GroupAnalysisDir 'PupilAvg_mean_' sessionType '.mat'],'Pupil','tvecPupil','avgPupil', '-v7.3');


% plot Pupil mean
%-------- plot Pupil for all channels
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/3 (screensize(4)-150)*2/3]); %(x,y,width,height) screensize(3)-100

for iRegion = 1:numRegion
%     if iRegion == 1; ipsiSide = 1; contraSide = 2;
%     else iRegion == 2; ipsiSide = 2; contraSide = 1;
%     end
    if iRegion == 1;yLim = [-1.5,1.5];
    elseif iRegion == 2; yLim = [-1.5,1.5];end
    subplot(3,numRegion,iRegion)
    mean = squeeze(Pupil.leftAvg(iRegion,:));
    sem  = nanstd(allPupil{1,iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.leftAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);
    title(regionNames{iRegion});
    if iRegion == 1;ylabel('Left normed amplitude');end
    
    subplot(3,numRegion,iRegion+numRegion)
    mean = squeeze(Pupil.rightAvg(iRegion,:));
    sem  = nanstd(allPupil{2,iRegion},[],1)/sqrt(size(allPupil{2,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.rightAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);
    if iRegion == 1;ylabel('Right normed amplitude');end
    
    subplot(3,numRegion,iRegion+2*numRegion)
    mean = squeeze(nanmean(Pupil.contrast{iRegion}));
    sem  = nanstd(Pupil.contrast{iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.l_rAvg(iRegion,:)), 'k','LineWidth', 1);
    xlim(xLim); ylim(yLim); %ylim([-0.8,0.8]);
    xlabel('Time [s]');
    if iRegion == 1; ylabel('Ipsi-Contra normed amplitude');end
  
    
%     subplot(3,numRegion,iRegion+2*numRegion)
%     mean = squeeze(nanmean(Pupil.l_r{iRegion}));
%     sem  = nanstd(Pupil.l_r{iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
%     shadedErrorBar(tvecPupil,mean,sem);
%     %plot(tvecPupil, squeeze(Pupil.l_rAvg(iRegion,:)), 'k','LineWidth', 1);
%     xlim(xLim); ylim(yLim); %ylim([-0.8,0.8]);
%     xlabel('Time [s]');
%     if iRegion == 1; ylabel('Left-Right normed amplitude');end
end

savefig(fig, [GroupAnalysisDir 'contrastAvg_Pupil_mean_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_Pupil_mean_' sessionType '.png']);
clear allPupil avgPupil Pupil
end