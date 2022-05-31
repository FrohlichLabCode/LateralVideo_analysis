% Select all visually responsive channels across sessions

%% prepare directory
clear all

cluster = 0;
skipRec = 1;
animalCode = '0172';
analysisType = 'PSTH';
folderSuffix = '';%'_validChns_new';
doPlot = 1;



if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];  
end

fileInfo   = dir([AnalysisDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto

% animal info
switch animalCode
    case '0169'
        regionNames = {'lLPl','lPPC','rPPC','rLPl','lVCEEG','rVCEEG'};
        numRegion   = numel(regionNames);
        numSpkRegion = numel(regionNames) - 2;
        numConds = 2;
        fullVideoSessions = {'039','042','045','048','051','054','058','062','064','069','073','077','083','086','092','100'};
        %if find(contains(fullVideoSessions,sessionID))>0; fullVideo = 1; else fullVideo = 0; end
        allIndex = [1:100];
        % -----------
        monocularIndex = [2,5,10,15,16,21,23,24,26,27,76,82,90,93,96]; %only analyze some sessions
        session2analyze = allIndex; %monocularIndex
    
    case '0172'
        regionNames = {'FC','LPl','PPC','VC'};
        numRegion   = numel(regionNames);
        numSpkRegion = numel(regionNames);
        numConds = 2;
        allIndex = [1:24];
        bigIndex = [3,5,7,9,11,13,15,17,20,23,26];
        smallIndex = setdiff(session2analyze,bigIndex);
        session2analyze = allIndex;
end
% -----------

% parameteres
twin = [-3 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-2,7];
baseTwin = [-2.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
tvecPSTH = is_load('D:\FerretData\0169\Analyzed\0169_LateralVideo_001_20180713\PSTH\zPSTH_mean.mat','timePSTH');
validChn = is_load('D:\FerretData\0169\Analyzed\0169_LateralVideo_001_20180713\PSTH\zPSTH_mean.mat','validChn'); %1 x numRegion cell
         
%%%%%%%%%%%%%%
%% Contrast %%
%%%%%%%%%%%%%%

% loop through each recording
if ~exist(join([GroupAnalysisDir 'contrastFigures/']),'dir'); mkdir(join([GroupAnalysisDir 'contrastFigures/']));end % to save figures for each session

allVisualPSTH = cell(numConds, numSpkRegion); % have to pre-allocate
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName = strsplit(recName,'_');
    sessionID = splitName{3}; 
    if ismember(str2num(sessionID),session2analyze)
    if find(contains(fullVideoSessions,sessionID))>0; continue;end % skip fullVideoSessions
    rootPSTHDir     = [AnalysisDir recName '/PSTH/'];
    frZ = is_load([rootPSTHDir 'zPSTH_mean_visual.mat'], 'frZ'); %numCond x numChn x tvecPSTH
    
    for iCond = 1:numConds
        for iRegion = 1:numSpkRegion
            sliceData = reshape(frZ(iCond, validChn{iRegion},:),length(validChn{iRegion}),[]); %numChn x tvecPSTH
            data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
            allVisualPSTH{iCond,iRegion} = [allVisualPSTH{iCond,iRegion}; data2Average];  % append the visual channel
        end
    end
    end
end
             
save([GroupAnalysisDir 'PSTHAll_visual.mat'],'tvecPSTH','allVisualPSTH','validChn', '-v7.3');

%%%%%%%%%%%%%%%%%%
%% Average PSTH %%
%%%%%%%%%%%%%%%%%%

for iCond = 1:numConds
    for iRegion = 1:numSpkRegion
        avgVisualPSTH(iCond,iRegion,:) = nanmean(allVisualPSTH{iCond,iRegion},1); % 1st dimension is channel
    end %median moves baseline down, so use mean!
end

PSTH.contraAvg(1:2,:) = squeeze(avgVisualPSTH(2,1:2,:));
PSTH.contraAvg(3:4,:) = squeeze(avgVisualPSTH(1,3:4,:));
PSTH.ipsiAvg(1:2,:)   = squeeze(avgVisualPSTH(1,1:2,:));
PSTH.ipsiAvg(3:4,:)   = squeeze(avgVisualPSTH(2,3:4,:));
PSTH.contrastAvg      = PSTH.contraAvg - PSTH.ipsiAvg;
save([GroupAnalysisDir 'PSTHAvg_visual.mat'],'PSTH','tvecPSTH','avgVisualPSTH','validChn', '-v7.3');


% plot PSTH
%-------- plot PSTH for all channels
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/3 (screensize(4)-150)*2/3]); %(x,y,width,height) screensize(3)-100
yLim = [-0.8,2];
for iRegion = 1:numSpkRegion
    subplot(3,4,iRegion)
    plot(tvecPSTH, squeeze(PSTH.contraAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim([-0.5,1]);
    title(regionNames{iRegion});
    xlim(xLim);
    if iRegion == 1;ylabel('Contra normed FR');end
    
    subplot(3,4,iRegion+numSpkRegion)
    plot(tvecPSTH, squeeze(PSTH.ipsiAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim([-0.5,1]);
    xlim(xLim);
    if iRegion == 1;ylabel('Ipsi normed FR');end
    
    subplot(3,4,iRegion+2*numSpkRegion)
    plot(tvecPSTH, squeeze(PSTH.contrastAvg(iRegion,:)),'k', 'LineWidth', 1);
    xlim(xLim); ylim([-0.5,1]); %ylim([-0.8,0.8]);
    xlabel('Time [s]');xlim(xLim);
    if iRegion == 1; ylabel('Contra-ipsi normed FR');end
end

savefig(fig, [GroupAnalysisDir 'contrastAvg_PSTH_visual.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_PSTH_visual.png']);


        


