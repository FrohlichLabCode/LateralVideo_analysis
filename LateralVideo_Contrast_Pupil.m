% Select all visually responsive channels across sessions

%% prepare directory
clear all

cluster = 0;
skipRec = 0;
animalCode = '0172';
analysisType = 'Pupil';
folderSuffix = '_validChns';
doPlot = 1;



if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/Pupil_clean/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];  
end

fileInfo   = dir([AnalysisDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto

% animal info
regionNames = {'lEye','rEye'};
numRegion   = numel(regionNames);
numConds = 2;
%fullVideoSessions = {'039','042','045','048','051','054','058','062','064','069','073','077','083','086','092','100','136'};
fullVideoSessions = {'039','042','045','048','051','054','058','061','062','064','067','068','069','073','074',...
                      '077','078','083','084','086','087','091','092','097','100'};

%if find(contains(fullVideoSessions,sessionID))>0; fullVideo = 1; else fullVideo = 0; end
session2analyze = [1:100];
% -----------
monocularIndex = [2,5,10,15,16,21,23,24,26,27,76,82,90,93,96]; %only analyze some sessions
oneEyeSessions = [1:28]; %not code to process this yet
twoEyeSessions = [29:100];
session2analyze = twoEyeSessions;
% -----------

% parameteres
twin = [-2 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-2,7];
baseTwin = [-2 0]; %<<<--- baseline to subtract from spectrogram 
tvecPupil = is_load('D:\FerretData\0169\Analyzed\0169_LateralVideo_001_20180713\Pupil\Pupil_contrast.mat','timeVec');
eyeChn = {'leftX','rightX','leftY','rightY','leftD','rightD'};

%%%%%%%%%%%%%%
%% Contrast %%
%%%%%%%%%%%%%%

% loop through each recording
if ~exist(join([GroupAnalysisDir 'contrastFigures/']),'dir'); mkdir(join([GroupAnalysisDir 'contrastFigures/']));end % to save figures for each session

allPupil = cell(numConds, numRegion); % have to pre-allocate
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName = strsplit(recName,'_');
    sessionID = splitName{3}; 
    if ismember(str2num(sessionID),session2analyze)
    if find(contains(fullVideoSessions,sessionID))>0; continue;end % skip fullVideoSessions
    rootPupilDir     = [AnalysisDir recName '/Pupil_clean/'];
    pupil = is_load([rootPupilDir 'Pupil_baselineNormed.mat'], 'evtDat_trialAvg'); %numCond x numChn x tvecPupil
    display(['processing ' rootPupilDir])
    for iCond = 1:numConds
        for iRegion = 1:numRegion
            sliceData = reshape(pupil(iCond,4+iRegion,:),1,[]); %numChn x tvecPupil
            data2Average = sliceData; % delete channels (2nd dimension) without spikes
            allPupil{iCond,iRegion} = [allPupil{iCond,iRegion}; data2Average];  % append the visual channel
        end
    end
    end
end
             
save([GroupAnalysisDir 'PupilAll.mat'],'tvecPupil','allPupil', '-v7.3');

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
Pupil.contrastAvg    = Pupil.contraAvg - Pupil.ipsiAvg;
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

xLim = [-2, 7];
for iRegion = 1:numRegion
    if iRegion == 1;yLim = [-5,15];
    elseif iRegion == 2; yLim = [-2.5,0.8];end
    
    subplot(3,numRegion,iRegion)
    mean = squeeze(Pupil.leftAvg(iRegion,:));
    sem  = nanstd(allPupil{1,iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.leftAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);
    title(regionNames{iRegion});
    if iRegion == 1;ylabel('Left normed amplitude');end
    
    subplot(3,numRegion,iRegion+numRegion)
    mean = squeeze(Pupil.leftAvg(iRegion,:));
    sem  = nanstd(allPupil{2,iRegion},[],1)/sqrt(size(allPupil{2,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.rightAvg(iRegion,:)),'k',  'LineWidth', 1);
    xlim(xLim); ylim(yLim);%ylim([-10,5]);
    if iRegion == 1;ylabel('Right normed amplitude');end
    
    subplot(3,numRegion,iRegion+2*numRegion)
    mean = squeeze(nanmedian(Pupil.l_r{iRegion},1));
    sem  = nanstd(Pupil.l_r{iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.l_rAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);%ylim([-1,2]); %ylim([-0.8,0.8]);
    xlabel('Time [s]');
    if iRegion == 1; ylabel('Left-right normed amplitude');end
end

savefig(fig, [GroupAnalysisDir 'contrastAvg_Pupil_median.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_Pupil_median.png']);

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
save([GroupAnalysisDir 'PupilAvg_mean.mat'],'Pupil','tvecPupil','avgPupil', '-v7.3');


% plot Pupil mean
%-------- plot Pupil for all channels
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/3 (screensize(4)-150)*2/3]); %(x,y,width,height) screensize(3)-100

xLim = [-2, 7];
for iRegion = 1:numRegion
%     if iRegion == 1; ipsiSide = 1; contraSide = 2;
%     else iRegion == 2; ipsiSide = 2; contraSide = 1;
%     end
    if iRegion == 1;yLim = [-5,15];
    elseif iRegion == 2; yLim = [-2.5,0.8];end
    subplot(3,numRegion,iRegion)
    mean = squeeze(Pupil.leftAvg(iRegion,:));
    sem  = nanstd(allPupil{1,iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.leftAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);
    title(regionNames{iRegion});
    if iRegion == 1;ylabel('Left normed amplitude');end
    
    subplot(3,numRegion,iRegion+numRegion)
    mean = squeeze(Pupil.leftAvg(iRegion,:));
    sem  = nanstd(allPupil{2,iRegion},[],1)/sqrt(size(allPupil{2,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.rightAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim(yLim);
    if iRegion == 1;ylabel('Right normed amplitude');end
    
    subplot(3,numRegion,iRegion+2*numRegion)
    mean = squeeze(nanmean(Pupil.l_r{iRegion}));
    sem  = nanstd(Pupil.l_r{iRegion},[],1)/sqrt(size(allPupil{1,iRegion},1));
    shadedErrorBar(tvecPupil,mean,sem);
    %plot(tvecPupil, squeeze(Pupil.l_rAvg(iRegion,:)), 'k','LineWidth', 1);
    xlim(xLim); ylim(yLim); %ylim([-0.8,0.8]);
    xlabel('Time [s]');
    if iRegion == 1; ylabel('Left-Right normed amplitude');end
end

savefig(fig, [GroupAnalysisDir 'contrastAvg_Pupil_mean.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_Pupil_mean.png']);
