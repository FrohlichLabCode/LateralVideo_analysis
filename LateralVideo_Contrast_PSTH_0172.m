% Select all visually responsive channels across sessions

%% prepare directory
clear all

cluster = 0;
skipRec = 0;
animalCode = '0172';
analysisType = 'PSTH_combineCondition';
visualUnits = 1;
doPlot = 1;

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

if visualUnits == 1
    suffix = '_visual';%'_validChns_new';
else
    suffix = '';
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
        %allIndex = [1:24];
        %bigIndex = [3,5,7,9,11,13,15,17,20,23,26];
        %smallIndex = setdiff(session2analyze,bigIndex);
        lateral = [1,2,4,6,8,10,12,14,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39];
        half    = [3,5,7,9,11,13,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,58,61,64,69];
        medial  = [40,42,43,45,46,48,49,51,52,54,55,57,59,60,62,63,65,66,68,70];
        sessionTypes =  {lateral, half, medial};
        sessionTypeNames = {'lateral', 'half', 'medial'};
end
% -----------

% parameteres
twin = [-3 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-2,7];
baseTwin = [-2.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
tvecPSTH = is_load('D:\FerretData\0169\Analyzed\0169_LateralVideo_001_20180713\PSTH\zPSTH_mean.mat','timePSTH');
         
%%%%%%%%%%%%%%
%% Contrast %%
%%%%%%%%%%%%%%

% loop through each recording
if ~exist(join([GroupAnalysisDir]),'dir'); mkdir(join([GroupAnalysisDir]));end % to save figures for each session
for i = 1:numel(sessionTypes)
    sessionIDs = sessionTypes{i};
    sessionType = sessionTypeNames{i};
    irec = 0;
    allVisualPSTH = cell(numConds, numSpkRegion); % have to pre-allocate
    
    for iSession = 1:numel(sessionIDs)
        sessionID = sessionIDs(iSession);
        sessionName = sprintf('%03d',sessionID); % add leading zeros        
        fileInfo = dir([AnalysisDir animalCode '_LateralVideo_' sessionName '*']); % detect files to load/convert  '_LateralVideo*'
        irec = irec +1;
        recName = fileInfo.name;
    %if ismember(str2num(sessionID),session2analyze)
    %if find(contains(fullVideoSessions,sessionID))>0; continue;end % skip fullVideoSessions
        rootPSTHDir     = [AnalysisDir recName '/' analysisType '/'];
        [frZ, validChn] = is_load([rootPSTHDir 'zPSTH_median' suffix '.mat'], 'frZ','validChn'); %numCond x numChn x tvecPSTH

    for iCond = 1:numConds
        for iRegion = 1:numSpkRegion
            sliceData = reshape(frZ(iCond, validChn{iRegion},:),length(validChn{iRegion}),[]); %numChn x tvecPSTH
            data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
            allVisualPSTH{iCond,iRegion} = [allVisualPSTH{iCond,iRegion}; data2Average];  % append the visual channel if exist
        end
    end
    end

             
save([GroupAnalysisDir 'PSTHAll_' sessionType suffix '.mat'],'tvecPSTH','allVisualPSTH','validChn', '-v7.3');

%%%%%%%%%%%%%%%%%%
%% Average PSTH %%
%%%%%%%%%%%%%%%%%%

for iRegion = 1:numSpkRegion
    for iCond = 1:numConds
        avgVisualPSTH(iCond,iRegion,:) = nanmean(allVisualPSTH{iCond,iRegion},1); % 1st dimension is channel
        semVisualPSTH(iCond,iRegion,:) = nanstd(allVisualPSTH{iCond,iRegion},1)/sqrt(size(allVisualPSTH{iCond,iRegion},1));
    end %median moves baseline down, so use mean!
    allVisualContrast{iRegion} = allVisualPSTH{2,iRegion} - allVisualPSTH{1,iRegion};
    avgVisualContrast(iRegion,:) = nanmean(allVisualContrast{iRegion},1);
    semVisualContrast(iRegion,:) = nanstd(allVisualContrast{iRegion},1)/sqrt(size(allVisualContrast{iRegion},1));
end

% PSTH.contraAvg(1:2,:) = squeeze(avgVisualPSTH(2,1:2,:));
% PSTH.contraAvg(3:4,:) = squeeze(avgVisualPSTH(1,3:4,:));
% PSTH.contraSem(1:2,:) = squeeze(semVisualPSTH(2,1:2,:));
% PSTH.contraSem(3:4,:) = squeeze(semVisualPSTH(1,3:4,:));
% 
% PSTH.ipsiAvg(1:2,:)   = squeeze(avgVisualPSTH(1,1:2,:));
% PSTH.ipsiAvg(3:4,:)   = squeeze(avgVisualPSTH(2,3:4,:));
% PSTH.ipsiSem(1:2,:)   = squeeze(semVisualPSTH(1,1:2,:));
% PSTH.ipsiSem(3:4,:)   = squeeze(semVisualPSTH(2,3:4,:));
% PSTH.contrastAvg      = PSTH.contraAvg - PSTH.ipsiAvg;

save([GroupAnalysisDir 'PSTHAvg_' sessionType suffix '.mat'],'tvecPSTH',...
    'avgVisualPSTH','semVisualPSTH','allVisualContrast','avgVisualContrast','semVisualContrast','validChn', '-v7.3');
clear allVisualPSTH

% plot PSTH
%-------- plot PSTH for all channels
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/3 (screensize(4)-150)*2/3]); %(x,y,width,height) screensize(3)-100
yLim = {[-0.7,0.7],[-1,1.5],[-1,1]}; % ylim for each session type
for iRegion = 1:numSpkRegion
    subplot(3,4,iRegion)
    toPlot = squeeze(avgVisualPSTH(2,iRegion,:));
    toPlot = smoothdata(toPlot,'gaussian',10);
    %plot(tvecPSTH, toPlot,'k', 'LineWidth', 1);
    shadedErrorBar(tvecPSTH, toPlot, squeeze(semVisualPSTH(2,iRegion,:)), '-k',4);
    xlim(xLim); ylim(yLim{i});
    title(regionNames{iRegion});
    xlim(xLim);
    if iRegion == 1;ylabel('Contra normed FR');end
    
    subplot(3,4,iRegion+numSpkRegion)
    toPlot = squeeze(avgVisualPSTH(1,iRegion,:));
    toPlot = smoothdata(toPlot,'gaussian',10);
    %plot(tvecPSTH, toPlot,'k', 'LineWidth', 1);
    shadedErrorBar(tvecPSTH, toPlot, squeeze(semVisualPSTH(1,iRegion,:)),'-k',2);
    xlim(xLim); ylim(yLim{i});
    xlim(xLim);
    if iRegion == 1;ylabel('Ipsi normed FR');end
    
    subplot(3,4,iRegion+2*numSpkRegion)
    toPlot = squeeze(avgVisualContrast(iRegion,:));
    toPlot = smoothdata(toPlot,'gaussian',10);
    %plot(tvecPSTH, toPlot,'k', 'LineWidth', 1);
    shadedErrorBar(tvecPSTH, toPlot, squeeze(semVisualContrast(iRegion,:)),'-k',2);
    xlim(xLim); ylim(yLim{i}); %ylim([-0.8,0.8]);
    xlabel('Time [s]');xlim(xLim);
    if iRegion == 1; ylabel('Contra-ipsi normed FR');end
end

savefig(fig, [GroupAnalysisDir 'contrastAvg_' sessionType '_PSTH' suffix '_sem.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_' sessionType '_PSTH' suffix '_sem.png']);

end
        


