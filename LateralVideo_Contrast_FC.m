%% This script is used after LateralVideo_FC_Normalization
% Input: plvAvgNormed.mat, GCAvgNormed.mat
% AH 5/2018: 0169 and 0172 are complete, need to adjust for 0168

%% prepare directory
clear
tic

cluster = 0;
skipRec = 0;
animalCode = '0169';
analysisType = 'FC';
folderSuffix = '_validChns'; %'_validChns_new'
doPlot = 1;


if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    baseDir = ['Z:\Individual\Angel\FerretData\'];
    PreprocessDir = [baseDir animalCode '/Preprocessed/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    BehavDatDir   = [baseDir animalCode '/behav/'];
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];  
end

fileInfo   = dir([AnalysisDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto

fullVideoSessions = {};
% animal info
switch animalCode
    case '0168'
        regionNames = {'LPl','PPC','VC'};
        regionPairs = {[1,2],[1,3],[2,3]}; 
    case '0169'
        regionNames = {'lLPl','lPPC','rPPC','rLPl','lVCEEG','rVCEEG'};
        regionPairs = {[1,2],[4,3],[1,4],[2,3]}; %{[2,5],[1,5],[1,2],[4,6],[3,6],[4,3],[1,4],[2,3],[5,6]};
        %{'lLPl->lPPC','lPPC->lLPl','rLPl->rPPC','rPPC->rLPl','lLPl->rLPl','rLPl->lLPl','lPPC->rPPC','rPPC->lPPC'};
        fullVideoSessions = {'039','042','045','048','051','054','058','061','062','064','067','068','069','073','074',...
                             '077','078','083','084','086','087','091','092','097','100'};
    case '0172'
        regionNames = {'FC','LPl','PPC','VC'};
        regionPairs = {[1,3],[2,3],[2,4],[3,4]};
        
end

numRegion   = numel(regionNames);
for iRegion = 1:numel(regionPairs)
    regionPairNames{iRegion}  = [regionNames{regionPairs{iRegion}(1)} '-' regionNames{regionPairs{iRegion}(2)}];
    regionPairNames1{iRegion} = [regionNames{regionPairs{iRegion}(1)} '->' regionNames{regionPairs{iRegion}(2)}];
    regionPairNames2{iRegion} = [regionNames{regionPairs{iRegion}(2)} '->' regionNames{regionPairs{iRegion}(1)}];
    regionPairNamesGC{2*iRegion-1} = regionPairNames1{iRegion};
    regionPairNamesGC{2*iRegion}   = regionPairNames2{iRegion};
    
end % eg. lLPl-lPPC

numAllSessions     = numel(fileInfo);
numFullSessions    = numel(fullVideoSessions);
numLateralSessions = numAllSessions - numFullSessions;

% parameteres
twin = [-3 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-2,7];
baseTwin = [-2.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
tvec = is_load([baseDir '0169\Analyzed\0169_LateralVideo_001_20180713\FC_validChns_new\lLPl-lVCEEG\specAll_Left.mat'],'tvec');
tvecGC = is_load([baseDir '0169\Analyzed\0169_LateralVideo_001_20180713\FC_validChns_new\lLPl-lVCEEG\GCAll_Left.mat'],'tvecGC');
foi  = is_load([baseDir '0169\Analyzed\0169_LateralVideo_001_20180713\FC_validChns_new\lLPl-lVCEEG\specAll_Left.mat'],'foi');
basetvecMask = tvec>= baseTwin(1) & tvec<= baseTwin(2);
basetvecMaskGC = tvecGC>= baseTwin(1) & tvecGC<= baseTwin(2);
tvecPSTH = is_load([baseDir '0169\Analyzed\0169_LateralVideo_001_20180713\PSTH\zPSTH_mean.mat'],'timePSTH');
%stimtvecMask = tvec>= stimTwin(1) & tvec<= stimTwin(2);
numfoi = length(foi);
numtvec = length(tvec);
numtvecGC = length(tvecGC);

linORlog = 2; lowFreq = 2; highFreq = 128;

if linORlog == 1
    fois = [2, 5:5:highFreq];
    tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
elseif linORlog == 2
    fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
    tickLabel = string(fois);
end
for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end



%%%%%%%%%%%%%%
%% Contrast %%
%%%%%%%%%%%%%%

% loop through each recording
if ~exist(join([GroupAnalysisDir 'contrastFigures/']),'dir'); mkdir(join([GroupAnalysisDir 'contrastFigures/']));end % to save figures for each session

for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName = strsplit(recName,'_');
    sessionID = splitName{3}; 
    if find(contains(fullVideoSessions,sessionID))>0; continue;end % skip fullVideoSessions
    rootAnalysisDir = [AnalysisDir recName '/' analysisType folderSuffix '/'];
    rootPSTHDir     = [AnalysisDir recName '/PSTH/'];
    
    % check if already processed, skip if wanted
    if length(dir([GroupAnalysisDir '*.fig'])) >= 1 %each condition generate a fig file
    fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
    if skipRec == 1; continue; end;end  % to skip already analyzed records
    fprintf('Normalizing %s \n',rootAnalysisDir');
    
    contrastAll.regionNames = regionNames; %{'lLPl','lPPC','rPPC','rLPl'}
    contrastAll.regionPairNames = regionPairNames; %{'lLPl-lPPC','rLPl-rPPC','lLPl-rLPl','lPPC-rPPC'}
    contrastAll.regionPairNamesGC = regionPairNamesGC;
    contrastAll.Spec(irec,1,:,:) = -(is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Left.mat'],'avgXNormed') - is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Right.mat'],'avgXNormed'));
    contrastAll.Spec(irec,2,:,:) = -(is_load([rootAnalysisDir 'lPPC-lVCEEG/funcCon_avg_Left.mat'],'avgXNormed') - is_load([rootAnalysisDir 'lPPC-lVCEEG/funcCon_avg_Right.mat'],'avgXNormed'));
    contrastAll.Spec(irec,3,:,:) = -(is_load([rootAnalysisDir 'rPPC-rVCEEG/funcCon_avg_Right.mat'],'avgXNormed') - is_load([rootAnalysisDir 'rPPC-rVCEEG/funcCon_avg_Left.mat'],'avgXNormed'));
    contrastAll.Spec(irec,4,:,:) = -(is_load([rootAnalysisDir 'rLPl-rVCEEG/funcCon_avg_Right.mat'],'avgXNormed') - is_load([rootAnalysisDir 'rLPl-rVCEEG/funcCon_avg_Left.mat'],'avgXNormed'));
%     contrast.lVCEEG = is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Left.mat'],'avgYNormed') - is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Right.mat'],'avgYNormed');
%     contrast.rVCEEG = is_load([rootAnalysisDir 'lLPl-rVCEEG/funcCon_avg_Right.mat'],'avgYNormed') - is_load([rootAnalysisDir 'lLPl-rVCEEG/funcCon_avg_Left.mat'],'avgYNormed');
    
    temp = is_load([rootPSTHDir 'zPSTH_mean.mat'],'toPlot');
    contrastAll.PSTH(irec,1:2,:) = squeeze(temp(2,1:2,:) - temp(1,1:2,:)); %Cond1-cond2 -> nRegion x tvecPSTH
    contrastAll.PSTH(irec,3:4,:) = squeeze(temp(1,3:4,:) - temp(2,3:4,:)); %Cond1-cond2 -> nRegion x tvecPSTH

    contrastAll.PLV(irec,1,:,:) = -(is_load([rootAnalysisDir 'lLPl-lPPC/plvAvgNormed.mat'],'plvAvgNormed_Left') - is_load([rootAnalysisDir 'lLPl-lPPC/plvAvgNormed.mat'],'plvAvgNormed_Right'));
    contrastAll.PLV(irec,2,:,:) = -(is_load([rootAnalysisDir 'rLPl-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Right') - is_load([rootAnalysisDir 'rLPl-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Left'));
    contrastAll.PLV(irec,3,:,:) = -(is_load([rootAnalysisDir 'lLPl-rLPl/plvAvgNormed.mat'],'plvAvgNormed_Left') - is_load([rootAnalysisDir 'lLPl-rLPl/plvAvgNormed.mat'],'plvAvgNormed_Right'));
    contrastAll.PLV(irec,4,:,:) = -(is_load([rootAnalysisDir 'lPPC-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Left') - is_load([rootAnalysisDir 'lPPC-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Right'));

    %{'lLPl->lPPC','lPPC->lLPl','rLPl->rPPC','rPPC->rLPl','rLPl->lLPl','lLPl->rLPl','lPPC->rPPC','rPPC->lPPC'}
    contrastAll.GC(irec,1,:,:) = -(is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left') - is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right'));
    contrastAll.GC(irec,2,:,:) = -(is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left') - is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right'));
    contrastAll.GC(irec,3,:,:) = -(is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left') - is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right'));
    contrastAll.GC(irec,4,:,:) = -(is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right') - is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left'));
    contrastAll.GC(irec,5,:,:) = -(is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right') - is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left'));
    contrastAll.GC(irec,6,:,:) = -(is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right') - is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left'));   
    contrastAll.GC(irec,7,:,:) = -(is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left') - is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right'));
    contrastAll.GC(irec,8,:,:) = -(is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right') - is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left'));

    if doPlot == 1
    %-------- plot contrast of PSTH for each session
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)*1/4]); %(x,y,width,height) screensize(3)-100

    for iRegion = 1:size(contrastAll.PSTH,2)
        subplot(1,4,iRegion)
        plot(tvecPSTH, squeeze(contrastAll.PSTH(irec,iRegion,:)), 'LineWidth', 1);
        xlim(xLim);ylim([-1.5,1.5]);
        title(regionNames{iRegion})
        xlabel('Time [s]');
        ylabel('Contra-ipsi normed FR');
    end
    %savefig(fig, [GroupAnalysisDir 'contrastFigures/' sessionID '_Normed_PSTH.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'contrastFigures/' sessionID '_Normed_PSTH.png']);
    
    %-------- plot contrast of spectrogram for each session    
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150) (screensize(4)-150)*1/3]); %(x,y,width,height) screensize(3)-100
    
    for iRegion = 1:size(contrastAll.PSTH,2)
        subplot(1,4,iRegion)
        imagesc(tvec,1:numel(foi),squeeze(contrastAll.Spec(irec,iRegion,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-2 2]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['Contra-ipsi normed power: ' regionNames{iRegion}],'FontSize',12)
    end
   
    colormap(jet)
    %savefig(fig, [GroupAnalysisDir 'contrastFigures/' sessionID '_Normed_Spec.fig'],'compact');    
    saveas(fig, [GroupAnalysisDir 'contrastFigures/' sessionID '_Normed_Spec.png']);
    
    %-------- plot contrast of FC for each session, first row PLV, 2,3 row GC    
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150)*numel(regionPairs)/4 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
    
    
    for iRegionPair = 1:numel(regionPairs)
        region_pair = [regionNames{regionPairs{iRegionPair}(1)} '_' regionNames{regionPairs{iRegionPair}(2)}];
        region_pair2 = [regionNames{regionPairs{iRegionPair}(2)} '_' regionNames{regionPairs{iRegionPair}(1)}];
        subplot(3,numel(regionPairs),iRegionPair)
        imagesc(tvec,1:numel(foi),squeeze(contrastAll.PLV(irec,iRegionPair,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['Contra-ipsi PLV: ' regionPairNames{iRegionPair}],'FontSize',12)

        subplot(3,numel(regionPairs),iRegionPair+numel(regionPairs))
        imagesc(tvecGC,1:numel(foi),squeeze(contrastAll.GC(irec,2*iRegionPair-1,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames1{iRegionPair}],'FontSize',12)

        subplot(3,numel(regionPairs),iRegionPair+2*numel(regionPairs))
        imagesc(tvecGC,1:numel(foi),squeeze(contrastAll.GC(irec,2*iRegionPair,:,:)));
        xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames2{iRegionPair}],'FontSize',12)

        colormap(jet)
        
    end
    %savefig(fig, [GroupAnalysisDir 'contrastFigures/' sessionID '_Normed_FC.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'contrastFigures/' sessionID '_Normed_FC.png']);
    close all
    end
end

save([GroupAnalysisDir 'contrastAll.mat'],'contrastAll','-v7.3');


%% 
%%%%%%%%%%%%
%% Contra %%
%%%%%%%%%%%%
% loop through each recording
if ~exist(join([GroupAnalysisDir 'contraFigures/']),'dir'); mkdir(join([GroupAnalysisDir 'contraFigures/']));end % to save figures for each session


for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName = strsplit(recName,'_');
    sessionID = splitName{3}; 
    if find(contains(fullVideoSessions,sessionID))>0; continue;end % skip fullVideoSessions
    rootAnalysisDir = [AnalysisDir recName '/' analysisType folderSuffix '/'];
    rootPSTHDir     = [AnalysisDir recName '/PSTH/'];
    
    % check if already processed, skip if wanted
    if length(dir([GroupAnalysisDir '*.fig'])) >= 1 %each condition generate a fig file
    fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
    if skipRec == 1; continue; end;end  % to skip already analyzed records
    fprintf('Normalizing %s \n',rootAnalysisDir');
    
    contraAll.regionNames = regionNames; %{'lLPl','lPPC','rPPC','rLPl'}
    contraAll.regionPairNames = regionPairNames; %{'lLPl-lPPC','rLPl-rPPC','lLPl-rLPl','lPPC-rPPC'}
    contraAll.regionPairNamesGC = regionPairNamesGC;
    contraAll.Spec(irec,1,:,:) = is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Right.mat'],'avgXNormed');
    contraAll.Spec(irec,2,:,:) = is_load([rootAnalysisDir 'lPPC-lVCEEG/funcCon_avg_Right.mat'],'avgXNormed');
    contraAll.Spec(irec,3,:,:) = is_load([rootAnalysisDir 'rPPC-rVCEEG/funcCon_avg_Left.mat'],'avgXNormed');
    contraAll.Spec(irec,4,:,:) = is_load([rootAnalysisDir 'rLPl-rVCEEG/funcCon_avg_Left.mat'],'avgXNormed');
%     contra.lVCEEG = is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Left.mat'],'avgYNormed') - is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Right.mat'],'avgYNormed');
%     contra.rVCEEG = is_load([rootAnalysisDir 'lLPl-rVCEEG/funcCon_avg_Right.mat'],'avgYNormed') - is_load([rootAnalysisDir 'lLPl-rVCEEG/funcCon_avg_Left.mat'],'avgYNormed');
    
    temp = is_load([rootPSTHDir 'zPSTH_mean.mat'],'toPlot');
    contraAll.PSTH(irec,1:2,:) = squeeze(temp(2,1:2,:)); %cond2 -> nRegion x tvecPSTH
    contraAll.PSTH(irec,3:4,:) = squeeze(temp(1,3:4,:)); %Cond1 -> nRegion x tvecPSTH

    contraAll.PLV(irec,1,:,:) = is_load([rootAnalysisDir 'lLPl-lPPC/plvAvgNormed.mat'],'plvAvgNormed_Right');
    contraAll.PLV(irec,2,:,:) = is_load([rootAnalysisDir 'rLPl-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Left');
    contraAll.PLV(irec,3,:,:) = is_load([rootAnalysisDir 'lLPl-rLPl/plvAvgNormed.mat'],'plvAvgNormed_Right');
    contraAll.PLV(irec,4,:,:) = is_load([rootAnalysisDir 'lPPC-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Right');

    %{'lLPl->lPPC','lPPC->lLPl','rLPl->rPPC','rPPC->rLPl','rLPl->lLPl','lLPl->rLPl','lPPC->rPPC','rPPC->lPPC'}
    contraAll.GC(irec,1,:,:) = is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right');
    contraAll.GC(irec,2,:,:) = is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right');
    contraAll.GC(irec,3,:,:) = is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right');
    contraAll.GC(irec,4,:,:) = is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left');
    contraAll.GC(irec,5,:,:) = is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left');
    contraAll.GC(irec,6,:,:) = is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left');   
    contraAll.GC(irec,7,:,:) = is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right');
    contraAll.GC(irec,8,:,:) = is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left');

    if doPlot == 1
    %-------- plot contra of PSTH for each session
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)*1/4]); %(x,y,width,height) screensize(3)-100

    for iRegion = 1:size(contraAll.PSTH,2)
        subplot(1,4,iRegion)
        plot(tvecPSTH, squeeze(contraAll.PSTH(irec,iRegion,:)), 'LineWidth', 1);
        xlim(xLim);ylim([-1.5,1.5]);
        title(regionNames{iRegion})
        xlabel('Time [s]');
        ylabel('Contra normed FR');
    end
    %savefig(fig, [GroupAnalysisDir 'contraFigures/' sessionID '_Normed_PSTH.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'contraFigures/' sessionID '_Normed_PSTH.png']);
    
    %-------- plot contra of spectrogram for each session    
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150) (screensize(4)-150)*1/3]); %(x,y,width,height) screensize(3)-100
    
    for iRegion = 1:size(contraAll.PSTH,2)
        subplot(1,4,iRegion)
        imagesc(tvec,1:numel(foi),squeeze(contraAll.Spec(irec,iRegion,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 2.5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['Contra normed power: ' regionNames{iRegion}],'FontSize',12)
    end
   
    colormap(jet)
    %savefig(fig, [GroupAnalysisDir 'contraFigures/' sessionID '_Normed_Spec.fig'],'compact');    
    saveas(fig, [GroupAnalysisDir 'contraFigures/' sessionID '_Normed_Spec.png']);
    
    %-------- plot contra of FC for each session, first row PLV, 2,3 row GC    
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150)*numel(regionPairs)/4 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
    
    
    for iRegionPair = 1:numel(regionPairs)
        region_pair = [regionNames{regionPairs{iRegionPair}(1)} '_' regionNames{regionPairs{iRegionPair}(2)}];
        region_pair2 = [regionNames{regionPairs{iRegionPair}(2)} '_' regionNames{regionPairs{iRegionPair}(1)}];
        subplot(3,numel(regionPairs),iRegionPair)
        imagesc(tvec,1:numel(foi),squeeze(contraAll.PLV(irec,iRegionPair,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['contra PLV: ' regionPairNames{iRegionPair}],'FontSize',12)

        subplot(3,numel(regionPairs),iRegionPair+numel(regionPairs))
        imagesc(tvecGC,1:numel(foi),squeeze(contraAll.GC(irec,2*iRegionPair-1,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames1{iRegionPair}],'FontSize',12)

        subplot(3,numel(regionPairs),iRegionPair+2*numel(regionPairs))
        imagesc(tvecGC,1:numel(foi),squeeze(contraAll.GC(irec,2*iRegionPair,:,:)));
        xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames2{iRegionPair}],'FontSize',12)

        colormap(jet)
        
    end
    %savefig(fig, [GroupAnalysisDir 'contraFigures/' sessionID '_Normed_FC.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'contraFigures/' sessionID '_Normed_FC.png']);
    close all
    end
end

save([GroupAnalysisDir 'contraAll.mat'],'contraAll','-v7.3');



%%
%%%%%%%%%%
%% Ipsi %%
%%%%%%%%%%

% loop through each recording
if ~exist(join([GroupAnalysisDir 'ipsiFigures/']),'dir'); mkdir(join([GroupAnalysisDir 'ipsiFigures/']));end % to save figures for each session


for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName = strsplit(recName,'_');
    sessionID = splitName{3}; 
    if find(contains(fullVideoSessions,sessionID))>0; continue;end % skip fullVideoSessions
    rootAnalysisDir = [AnalysisDir recName '/' analysisType folderSuffix '/'];
    rootPSTHDir     = [AnalysisDir recName '/PSTH/'];
    
    % check if already processed, skip if wanted
    if length(dir([GroupAnalysisDir '*.fig'])) >= 1 %each condition generate a fig file
    fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
    if skipRec == 1; continue; end;end  % to skip already analyzed records
    fprintf('Normalizing %s \n',rootAnalysisDir');
    
    ipsiAll.regionNames = regionNames; %{'lLPl','lPPC','rPPC','rLPl'}
    ipsiAll.regionPairNames = regionPairNames; %{'lLPl-lPPC','rLPl-rPPC','lLPl-rLPl','lPPC-rPPC'}
    ipsiAll.regionPairNamesGC = regionPairNamesGC;
    ipsiAll.Spec(irec,1,:,:) = is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Left.mat'],'avgXNormed');
    ipsiAll.Spec(irec,2,:,:) = is_load([rootAnalysisDir 'lPPC-lVCEEG/funcCon_avg_Left.mat'],'avgXNormed');
    ipsiAll.Spec(irec,3,:,:) = is_load([rootAnalysisDir 'rPPC-rVCEEG/funcCon_avg_Right.mat'],'avgXNormed');
    ipsiAll.Spec(irec,4,:,:) = is_load([rootAnalysisDir 'rLPl-rVCEEG/funcCon_avg_Right.mat'],'avgXNormed');
    
%     ipsi.lVCEEG = is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Left.mat'],'avgYNormed') - is_load([rootAnalysisDir 'lLPl-lVCEEG/funcCon_avg_Right.mat'],'avgYNormed');
%     ipsi.rVCEEG = is_load([rootAnalysisDir 'lLPl-rVCEEG/funcCon_avg_Right.mat'],'avgYNormed') - is_load([rootAnalysisDir 'lLPl-rVCEEG/funcCon_avg_Left.mat'],'avgYNormed');
    
    temp = is_load([rootPSTHDir 'zPSTH_mean.mat'],'toPlot');
    ipsiAll.PSTH(irec,1:2,:) = squeeze(temp(1,1:2,:)); %Cond1-cond2 -> nRegion x tvecPSTH
    ipsiAll.PSTH(irec,3:4,:) = squeeze(temp(2,3:4,:)); %Cond1-cond2 -> nRegion x tvecPSTH

    ipsiAll.PLV(irec,1,:,:) = is_load([rootAnalysisDir 'lLPl-lPPC/plvAvgNormed.mat'],'plvAvgNormed_Left');
    ipsiAll.PLV(irec,2,:,:) = is_load([rootAnalysisDir 'rLPl-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Right');
    ipsiAll.PLV(irec,3,:,:) = is_load([rootAnalysisDir 'lLPl-rLPl/plvAvgNormed.mat'],'plvAvgNormed_Left');
    ipsiAll.PLV(irec,4,:,:) = is_load([rootAnalysisDir 'lPPC-rPPC/plvAvgNormed.mat'],'plvAvgNormed_Left');

    %{'lLPl->lPPC','lPPC->lLPl','rLPl->rPPC','rPPC->rLPl','rLPl->lLPl','lLPl->rLPl','lPPC->rPPC','rPPC->lPPC'}
    ipsiAll.GC(irec,1,:,:) = is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left');
    ipsiAll.GC(irec,2,:,:) = is_load([rootAnalysisDir 'lLPl-lPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Left');
    ipsiAll.GC(irec,3,:,:) = is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left');
    ipsiAll.GC(irec,4,:,:) = is_load([rootAnalysisDir 'lLPl-rLPl/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right');
    ipsiAll.GC(irec,5,:,:) = is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Right');
    ipsiAll.GC(irec,6,:,:) = is_load([rootAnalysisDir 'rLPl-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right');   
    ipsiAll.GC(irec,7,:,:) = is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left');
    ipsiAll.GC(irec,8,:,:) = is_load([rootAnalysisDir 'lPPC-rPPC/GCAvgNormed.mat'],'normedGC_YtoX_Avg_Right');

    if doPlot == 1
    %-------- plot ipsi of PSTH for each session
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)*1/4]); %(x,y,width,height) screensize(3)-100

    for iRegion = 1:size(ipsiAll.PSTH,2)
        subplot(1,4,iRegion)
        plot(tvecPSTH, squeeze(ipsiAll.PSTH(irec,iRegion,:)), 'LineWidth', 1);
        xlim(xLim);ylim([-1.5,1.5]);
        title(regionNames{iRegion})
        xlabel('Time [s]');
        ylabel('Ipsi normed FR');
    end
    %savefig(fig, [GroupAnalysisDir 'ipsiFigures/' sessionID '_Normed_PSTH.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'ipsiFigures/' sessionID '_Normed_PSTH.png']);
    
    %-------- plot ipsi of spectrogram for each session    
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150) (screensize(4)-150)*1/3]); %(x,y,width,height) screensize(3)-100
    
    for iRegion = 1:size(ipsiAll.PSTH,2)
        subplot(1,4,iRegion)
        imagesc(tvec,1:numel(foi),squeeze(ipsiAll.Spec(irec,iRegion,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 2.5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['Ipsi normed power: ' regionNames{iRegion}],'FontSize',12)
    end
   
    colormap(jet)
    %savefig(fig, [GroupAnalysisDir 'ipsiFigures/' sessionID '_Normed_Spec.fig'],'compact');    
    saveas(fig, [GroupAnalysisDir 'ipsiFigures/' sessionID '_Normed_Spec.png']);
    
    %-------- plot ipsi of FC for each session, first row PLV, 2,3 row GC    
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-150)*numel(regionPairs)/4 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
    
    
    for iRegionPair = 1:numel(regionPairs)
        region_pair = [regionNames{regionPairs{iRegionPair}(1)} '_' regionNames{regionPairs{iRegionPair}(2)}];
        region_pair2 = [regionNames{regionPairs{iRegionPair}(2)} '_' regionNames{regionPairs{iRegionPair}(1)}];
        subplot(3,numel(regionPairs),iRegionPair)
        imagesc(tvec,1:numel(foi),squeeze(ipsiAll.PLV(irec,iRegionPair,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['Ipsi PLV: ' regionPairNames{iRegionPair}],'FontSize',12)

        subplot(3,numel(regionPairs),iRegionPair+numel(regionPairs))
        imagesc(tvecGC,1:numel(foi),squeeze(ipsiAll.GC(irec,2*iRegionPair-1,:,:)));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames1{iRegionPair}],'FontSize',12)

        subplot(3,numel(regionPairs),iRegionPair+2*numel(regionPairs))
        imagesc(tvecGC,1:numel(foi),squeeze(ipsiAll.GC(irec,2*iRegionPair,:,:)));
        xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);xlim(xLim);
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames2{iRegionPair}],'FontSize',12)

        colormap(jet)
        
    end
    %savefig(fig, [GroupAnalysisDir 'ipsiFigures/' sessionID '_Normed_FC.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'ipsiFigures/' sessionID '_Normed_FC.png']);
    close all
    end
end

save([GroupAnalysisDir 'ipsiAll.mat'],'ipsiAll','-v7.3');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%% Session average %% Contrast
%%%%%%%%%%%%%%%%%%%%%
% if needed, load contrastAll files
load('D:\FerretData\0169\GroupAnalysis\lateralVideoSessions\contraAll.mat')
load('D:\FerretData\0169\GroupAnalysis\lateralVideoSessions\contrastAll.mat')
load('D:\FerretData\0169\GroupAnalysis\lateralVideoSessions\ipsiAll.mat')


fullIndex = cellfun(@str2num,fullVideoSessions); % 25 sessions
allIndex  = [1:100];
%lateralIndex = setdiff(allIndex,fullIndex); % 75 sessions

smallMonocularIndex = [2,5,10,15,16,21,23,24,26,27,76,82,90,93,96];
allMonocularIndex = [2,4,5,10,15,16,21,23,24,26:44,46,47,49,50,53,...
                     55:57,60,63,65,70:72,75,76,80:82,85,88,89,90,93,99,99];
lateralIndex = allMonocularIndex; % binocular sessions

% calculate average of all sessions
contrastAvg.regionNames = regionNames;
contrastAvg.regionPairNames = regionPairNames;
contrastAvg.regionPairNamesGC = regionPairNamesGC;
contrastAvg.Spec = squeeze(nanmedian(contrastAll.Spec(lateralIndex,:,:,:),1)); % 1st dimension is session
contrastAvg.PSTH = squeeze(nanmedian(contrastAll.PSTH(lateralIndex,:,:,:),1));
contrastAvg.PLV  = squeeze(nanmedian(contrastAll.PLV(lateralIndex,:,:,:),1));
contrastAvg.GC   = squeeze(nanmedian(contrastAll.GC(lateralIndex,:,:,:),1));

save([GroupAnalysisDir 'contrastAvg.mat'],'contrastAvg','-v7.3');

% plot contrastAvg
%-------- plot contrast of PSTH for each session
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)*1/4]); %(x,y,width,height) screensize(3)-100

for iRegion = 1:size(contrastAvg.PSTH,1)
    subplot(1,4,iRegion)
    plot(tvecPSTH, squeeze(contrastAvg.PSTH(iRegion,:,:)), 'LineWidth', 1);
    xlim(twin); ylim([-0.2,0.3]);
    title(regionNames{iRegion})
    xlabel('Time [s]');xlim(xLim);
    ylabel('Contra-ipsi normed FR');
end

savefig(fig, [GroupAnalysisDir 'contrastAvg_Normed_PSTH.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_Normed_PSTH.png']);


%-------- plot contrast of spectrogram for each session    
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150) (screensize(4)-150)*1/3]); %(x,y,width,height) screensize(3)-100
for iRegion = 1:size(contrastAvg.Spec,1)
    
    subplot(1,4,iRegion)
    imagesc(tvec,1:numel(foi),squeeze(contrastAvg.Spec(iRegion,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.3 0.3]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['Contra-ipsi power: ' regionNames{iRegion}],'FontSize',12)

end
colormap(jet)

savefig(fig, [GroupAnalysisDir 'contrastAvg_Normed_Spec.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_Normed_Spec.png']);
    
%-------- plot contrast of FC for each session    
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)*numel(regionPairs)/4 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
    
    
for iRegionPair = 1:numel(regionPairs)
    region_pair = [regionNames{regionPairs{iRegionPair}(1)} '_' regionNames{regionPairs{iRegionPair}(2)}];
    region_pair2 = [regionNames{regionPairs{iRegionPair}(2)} '_' regionNames{regionPairs{iRegionPair}(1)}];
    subplot(3,numel(regionPairs),iRegionPair)
    imagesc(tvec,1:numel(foi),squeeze(contrastAvg.PLV(iRegionPair,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['Contra-ipsi PLV: ' regionPairNames{iRegionPair}],'FontSize',12)

    subplot(3,numel(regionPairs),iRegionPair+numel(regionPairs))
    imagesc(tvecGC,1:numel(foi),squeeze(contrastAvg.GC(2*iRegionPair-1,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames1{iRegionPair}],'FontSize',12)

    subplot(3,numel(regionPairs),iRegionPair+2*numel(regionPairs))
    imagesc(tvecGC,1:numel(foi),squeeze(contrastAvg.GC(2*iRegionPair,:,:)));
    xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames2{iRegionPair}],'FontSize',12)

    colormap(jet)

end
savefig(fig, [GroupAnalysisDir 'contrastAvg_Normed_FC.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contrastAvg_Normed_FC.png']);


%%%%%%%%%%%%%%%%%%%%%
%% Session average %% Contra
%%%%%%%%%%%%%%%%%%%%%

% fullIndex = cellfun(@str2num,fullVideoSessions); % 25 sessions
% allIndex  = [1:100];
% lateralIndex = setdiff(allIndex,fullIndex); % 75 sessions

% calculate average of all sessions
contraAvg.regionNames = regionNames;
contraAvg.regionPairNames = regionPairNames;
contraAvg.regionPairNamesGC = regionPairNamesGC;
contraAvg.Spec = squeeze(nanmedian(contraAll.Spec(lateralIndex,:,:,:),1)); % 1st dimension is session
contraAvg.PSTH = squeeze(nanmedian(contraAll.PSTH(lateralIndex,:,:,:),1));
contraAvg.PLV  = squeeze(nanmedian(contraAll.PLV(lateralIndex,:,:,:),1));
contraAvg.GC   = squeeze(nanmedian(contraAll.GC(lateralIndex,:,:,:),1));

save([GroupAnalysisDir 'contraAvg.mat'],'contraAvg','-v7.3');

% plot contraAvg
%-------- plot contra of PSTH for each session
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)*1/4]); %(x,y,width,height) screensize(3)-100

for iRegion = 1:size(contraAvg.PSTH,1)
    subplot(1,4,iRegion)
    plot(tvecPSTH, squeeze(contraAvg.PSTH(iRegion,:,:)), 'LineWidth', 1);
    xlim(twin); ylim([-0.2,0.3]);
    title(regionNames{iRegion})
    xlabel('Time [s]');xlim(xLim);
    ylabel('contra normed FR');
end

savefig(fig, [GroupAnalysisDir 'contraAvg_Normed_PSTH.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contraAvg_Normed_PSTH.png']);


%-------- plot contra of spectrogram for each session    
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150) (screensize(4)-150)*1/3]); %(x,y,width,height) screensize(3)-100
for iRegion = 1:size(contraAvg.Spec,1)
    
    subplot(1,4,iRegion)
    imagesc(tvec,1:numel(foi),squeeze(contraAvg.Spec(iRegion,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([0.7 1.2]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['contra power: ' regionNames{iRegion}],'FontSize',12)

end
colormap(jet)

savefig(fig, [GroupAnalysisDir 'contraAvg_Normed_Spec.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contraAvg_Normed_Spec.png']);
    
%-------- plot contra of FC for each session    
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)*numel(regionPairs)/4 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
    
    
for iRegionPair = 1:numel(regionPairs)
    region_pair = [regionNames{regionPairs{iRegionPair}(1)} '_' regionNames{regionPairs{iRegionPair}(2)}];
    region_pair2 = [regionNames{regionPairs{iRegionPair}(2)} '_' regionNames{regionPairs{iRegionPair}(1)}];
    subplot(3,numel(regionPairs),iRegionPair)
    imagesc(tvec,1:numel(foi),squeeze(contraAvg.PLV(iRegionPair,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['contra PLV: ' regionPairNames{iRegionPair}],'FontSize',12)

    subplot(3,numel(regionPairs),iRegionPair+numel(regionPairs))
    imagesc(tvecGC,1:numel(foi),squeeze(contraAvg.GC(2*iRegionPair-1,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames1{iRegionPair}],'FontSize',12)

    subplot(3,numel(regionPairs),iRegionPair+2*numel(regionPairs))
    imagesc(tvecGC,1:numel(foi),squeeze(contraAvg.GC(2*iRegionPair,:,:)));
    xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames2{iRegionPair}],'FontSize',12)

    colormap(jet)

end
savefig(fig, [GroupAnalysisDir 'contraAvg_Normed_FC.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'contraAvg_Normed_FC.png']);



%%%%%%%%%%%%%%%%%%%%%
%% Session average %% Ipsi
%%%%%%%%%%%%%%%%%%%%%

% fullIndex = cellfun(@str2num,fullVideoSessions); % 25 sessions
% allIndex  = [1:100];
% lateralIndex = setdiff(allIndex,fullIndex); % 75 sessions

% calculate average of all sessions
ipsiAvg.regionNames = regionNames;
ipsiAvg.regionPairNames = regionPairNames;
ipsiAvg.regionPairNamesGC = regionPairNamesGC;
ipsiAvg.Spec = squeeze(nanmedian(ipsiAll.Spec(lateralIndex,:,:,:),1)); % 1st dimension is session
ipsiAvg.PSTH = squeeze(nanmedian(ipsiAll.PSTH(lateralIndex,:,:,:),1));
ipsiAvg.PLV  = squeeze(nanmedian(ipsiAll.PLV(lateralIndex,:,:,:),1));
ipsiAvg.GC   = squeeze(nanmedian(ipsiAll.GC(lateralIndex,:,:,:),1));

save([GroupAnalysisDir 'ipsiAvg.mat'],'ipsiAvg','-v7.3');

% plot ipsiAvg
%-------- plot ipsi of PSTH for each session
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)*1/4]); %(x,y,width,height) screensize(3)-100

for iRegion = 1:size(ipsiAvg.PSTH,1)
    subplot(1,4,iRegion)
    plot(tvecPSTH, squeeze(ipsiAvg.PSTH(iRegion,:,:)), 'LineWidth', 1);
    xlim(twin); ylim([-0.2,0.3]);
    title(regionNames{iRegion})
    xlabel('Time [s]');xlim(xLim);
    ylabel('Ipsi normed FR');
end

savefig(fig, [GroupAnalysisDir 'ipsiAvg_Normed_PSTH.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'ipsiAvg_Normed_PSTH.png']);


%-------- plot ipsi of spectrogram for each session    
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150) (screensize(4)-150)*1/3]); %(x,y,width,height) screensize(3)-100
for iRegion = 1:size(ipsiAvg.Spec,1)
    
    subplot(1,4,iRegion)
    imagesc(tvec,1:numel(foi),squeeze(ipsiAvg.Spec(iRegion,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([0.7 1.2]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['Ipsi power: ' regionNames{iRegion}],'FontSize',12)

end
colormap(jet)

savefig(fig, [GroupAnalysisDir 'ipsiAvg_Normed_Spec.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'ipsiAvg_Normed_Spec.png']);
    
%-------- plot ipsi of FC for each session    
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)*numel(regionPairs)/4 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
    
    
for iRegionPair = 1:numel(regionPairs)
    region_pair = [regionNames{regionPairs{iRegionPair}(1)} '_' regionNames{regionPairs{iRegionPair}(2)}];
    region_pair2 = [regionNames{regionPairs{iRegionPair}(2)} '_' regionNames{regionPairs{iRegionPair}(1)}];
    subplot(3,numel(regionPairs),iRegionPair)
    imagesc(tvec,1:numel(foi),squeeze(ipsiAvg.PLV(iRegionPair,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['Ipsi PLV: ' regionPairNames{iRegionPair}],'FontSize',12)

    subplot(3,numel(regionPairs),iRegionPair+numel(regionPairs))
    imagesc(tvecGC,1:numel(foi),squeeze(ipsiAvg.GC(2*iRegionPair-1,:,:)));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames1{iRegionPair}],'FontSize',12)

    subplot(3,numel(regionPairs),iRegionPair+2*numel(regionPairs))
    imagesc(tvecGC,1:numel(foi),squeeze(ipsiAvg.GC(2*iRegionPair,:,:)));
    xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    caxis([-0.5 0.5]);xlim(xLim);
    cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames2{iRegionPair}],'FontSize',12)

    colormap(jet)

end
savefig(fig, [GroupAnalysisDir 'ipsiAvg_Normed_FC.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'ipsiAvg_Normed_FC.png']);















% %% selecting half screen sessions
% load('D:\FerretData\0169\GroupAnalysis\contraAll.mat')
% %halfIndex = [1,22,25,61,67,68,74,78,84,87,91,97]; %waiting for cluster
% binocularIndex = [2,5,7,10,14,15,16,18,19,21,23,24,26,27,76,82,90,93,96];
% halfIndex = binocularIndex;
% % calculate average of all sessions
% contrastAvg2.regionNames = regionNames;
% contrastAvg2.regionPairNames = regionPairNames;
% contrastAvg2.regionPairNamesGC = regionPairNamesGC;
% contrastAvg2.PSTH = squeeze(nanmedian(contrastAll.PSTH(halfIndex,:,:,:),1)); % 1st dimension is session
% contrastAvg2.Spec = squeeze(nanmedian(contrastAll.Spec(halfIndex,:,:,:),1));
% contrastAvg2.PLV  = squeeze(nanmedian(contrastAll.PLV(halfIndex,:,:,:),1));
% contrastAvg2.GC   = squeeze(nanmedian(contrastAll.GC(halfIndex,:,:,:),1));
% 
% save([GroupAnalysisDir 'contrastAvg2.mat'],'contrastAvg2','-v7.3');
% 
% % plot contrastAvg2
% %-------- plot contrast of PSTH for each session
% screensize = get( groot, 'Screensize' );
% fig = figure('Position',[10 50 (screensize(3)-150)/2 (screensize(4)-150)*1/4]); %(x,y,width,height) screensize(3)-100
% 
% for iRegion = 1:size(contrastAvg2.PSTH,1)
%     subplot(1,4,iRegion)
%     plot(tvecPSTH, squeeze(contrastAvg2.PSTH(iRegion,:,:)), 'LineWidth', 1);
%     xlim(twin); ylim([-0.2,0.3]);
%     title(regionNames{iRegion})
%     xlabel('Time [s]');xlim(xLim);
%     ylabel('Contra-ipsi normed FR');
% end
% 
% savefig(fig, [GroupAnalysisDir 'contrastAvg2_Normed_PSTH.fig'],'compact');
% saveas(fig, [GroupAnalysisDir 'contrastAvg2_Normed_PSTH.png']);
% 
% 
% %-------- plot contrast of spectrogram for each session    
% screensize = get( groot, 'Screensize' );
% fig = figure('Position',[10 50 (screensize(3)-150) (screensize(4)-150)*1/3]); %(x,y,width,height) screensize(3)-100
% for iRegion = 1:size(contrastAvg2.Spec,1)
%     
%     subplot(1,4,iRegion)
%     imagesc(tvec,1:numel(foi),squeeze(contrastAvg2.Spec(iRegion,:,:)));
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
%     ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%     caxis([-0.5 0.5]);xlim(xLim);
%     cl = colorbar('northoutside'); ylabel(cl,['Contra-ipsi power: ' regionNames{iRegion}],'FontSize',12)
% 
% end
% colormap(jet)
% 
% savefig(fig, [GroupAnalysisDir 'contrastAvg2_Normed_Spec.fig'],'compact');
% saveas(fig, [GroupAnalysisDir 'contrastAvg2_Normed_Spec.png']);
%     
% %-------- plot contrast of FC for each session    
% screensize = get( groot, 'Screensize' );
% fig = figure('Position',[10 50 (screensize(3)-150)*numel(regionPairs)/4 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
%     
%     
% for iRegionPair = 1:numel(regionPairs)
%     region_pair = [regionNames{regionPairs{iRegionPair}(1)} '_' regionNames{regionPairs{iRegionPair}(2)}];
%     region_pair2 = [regionNames{regionPairs{iRegionPair}(2)} '_' regionNames{regionPairs{iRegionPair}(1)}];
%     subplot(3,numel(regionPairs),iRegionPair)
%     imagesc(tvec,1:numel(foi),squeeze(contrastAvg2.PLV(iRegionPair,:,:)));
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
%     ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%     caxis([-0.5 0.5]);xlim(xLim);
%     cl = colorbar('northoutside'); ylabel(cl,['Contra-ipsi PLV: ' regionPairNames{iRegionPair}],'FontSize',12)
% 
%     subplot(3,numel(regionPairs),iRegionPair+numel(regionPairs))
%     imagesc(tvecGC,1:numel(foi),squeeze(contrastAvg2.GC(2*iRegionPair-1,:,:)));
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
%     ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%     caxis([-0.5 0.5]);xlim(xLim);
%     cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames1{iRegionPair}],'FontSize',12)
% 
%     subplot(3,numel(regionPairs),iRegionPair+2*numel(regionPairs))
%     imagesc(tvecGC,1:numel(foi),squeeze(contrastAvg2.GC(2*iRegionPair,:,:)));
%     xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
%     ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%     caxis([-0.5 0.5]);xlim(xLim);
%     cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNames2{iRegionPair}],'FontSize',12)
% 
%     colormap(jet)
% 
% end
% savefig(fig, [GroupAnalysisDir 'contrastAvg2_Normed_FC.fig'],'compact');
% saveas(fig, [GroupAnalysisDir 'contrastAvg2_Normed_FC.png']);
