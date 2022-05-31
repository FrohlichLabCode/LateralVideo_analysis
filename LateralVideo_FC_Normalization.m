%% This script is used after LateralVideo_FunConn.m to get baseline normalized
% PLV and GC (ideally on trial level, but that needs reprocess all session,
% so here I took a shortcut by normalizing PLV and GC for each channel pair
% and then average. 

%% prepare directory
clear
tic

cluster = 0;
skipRec = 1;
animalCode = '0169';
analysisType = 'FC';
folderSuffix = '_validChns_new';


if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];  
end

fileInfo   = dir([AnalysisDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto


% animal info
regionNames = {'lLPl','lPPC','rPPC','rLPl','lVCEEG','rVCEEG'};
numRegion   = numel(regionNames);
regionPairs = {[2,5],[1,5],[1,2],[4,6],[3,6],[4,3],[1,4],[2,3],[5,6]};
for i = 1:numel(regionPairs)
    regionPairNames{i}  = [regionNames{regionPairs{i}(1)} '-' regionNames{regionPairs{i}(2)}];
    regionPairNames1{i} = [regionNames{regionPairs{i}(1)} '->' regionNames{regionPairs{i}(2)}];
    regionPairNames2{i} = [regionNames{regionPairs{i}(2)} '->' regionNames{regionPairs{i}(1)}];
end % eg. lLPl-lPPC
fullVideoSessions = {'039','042','045','048','051','054','058','061','062','064','067','068','069','073','074',...
                     '077','078','083','084','086','087','091','092','097','100'};
numAllSessions     = numel(fileInfo);
numFullSessions    = numel(fullVideoSessions);
numLateralSessions = numAllSessions - numFullSessions;

% parameteres
twin = [-3 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
baseTwin = [-2.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
tvec = is_load('D:\FerretData\0169\Analyzed\0169_LateralVideo_001_20180713\FC_validChns_new\lLPl-lVCEEG\specAll_Left.mat','tvec');
tvecGC = is_load('D:\FerretData\0169\Analyzed\0169_LateralVideo_001_20180713\FC_validChns_new\lLPl-lVCEEG\GCAll_Left.mat','tvecGC');
foi  = is_load('D:\FerretData\0169\Analyzed\0169_LateralVideo_001_20180713\FC_validChns_new\lLPl-lVCEEG\specAll_Left.mat','foi');
basetvecMask = tvec>= baseTwin(1) & tvec<= baseTwin(2);
basetvecMaskGC = tvecGC>= baseTwin(1) & tvecGC<= baseTwin(2);
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

%% 
% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName = strsplit(recName,'_');
    sessionID = splitName{3}; 
    if find(contains(fullVideoSessions,sessionID))>0 ; continue;end % skip fullVideoSessions
    if cluster == 0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180820', 'InputFormat', 'yyyyMMdd'); continue;end

    for iRegionPair = 1:numel(regionPairs)
        regionPairName = regionPairNames{iRegionPair};
        rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/' regionPairName '/']; %<---------------
        % check if already processed, skip if wanted
        if length(dir([rootAnalysisDir '*.fig'])) >= 3 %each condition generate a fig file
        fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
        if skipRec == 1; continue; end  % to skip already analyzed records
    end
    fprintf('Normalizing %s \n',rootAnalysisDir');
        % first calculate normalized PLV and GC from plvAll_Left.mat
        % add real to avoid sometimes GC being complex number (real+0i),
        % have checked that real(plvAll) == plvAll 
        rawMat = real(is_load([rootAnalysisDir 'plvAll_Left.mat'],'plvAll')); %nChan x nChan x foi x tvec
        avgbaseAmp = repmat(nanmean(rawMat(:,:,:,basetvecMask),4),[1,1,1,size(rawMat,4)]); % checked, last dimension is the same values
        stdbaseAmp = repmat(nanstd(rawMat(:,:,:,basetvecMask),0,4),[1,1,1,size(rawMat,4)]);%nanstd(X,flag,dim)
        plvAllNormed_Left = (rawMat - avgbaseAmp)./stdbaseAmp;
        plvAvgNormed_Left = squeeze(nanmean(nanmean(plvAllNormed_Left ,1),2));

        rawMat = real(is_load([rootAnalysisDir 'plvAll_Right.mat'],'plvAll')); %nChan x nChan x foi x tvec
        avgbaseAmp = repmat(nanmean(rawMat(:,:,:,basetvecMask),4),[1,1,1,size(rawMat,4)]); % checked, last dimension is the same values
        stdbaseAmp = repmat(nanstd(rawMat(:,:,:,basetvecMask),0,4),[1,1,1,size(rawMat,4)]);%nanstd(X,flag,dim)
        plvAllNormed_Right = (rawMat - avgbaseAmp)./stdbaseAmp;
        plvAvgNormed_Right = squeeze(nanmean(nanmean(plvAllNormed_Right ,1),2));
        
        % GC Left video
        rawMat = real(is_load([rootAnalysisDir 'GCAll_Left.mat'],'GC_XtoY_All')); %nChan x nChan x foi x tvec
        avgbaseAmp = repmat(nanmean(rawMat(:,:,:,basetvecMaskGC),4),[1,1,1,size(rawMat,4)]); % checked, last dimension is the same values
        stdbaseAmp = repmat(nanstd(rawMat(:,:,:,basetvecMaskGC),0,4),[1,1,1,size(rawMat,4)]);%nanstd(X,flag,dim)
        normedGC_XtoY_All_Left = (rawMat - avgbaseAmp)./stdbaseAmp;
        normedGC_XtoY_Avg_Left = squeeze(nanmean(nanmean(normedGC_XtoY_All_Left ,1),2));

        rawMat = real(is_load([rootAnalysisDir 'GCAll_Left.mat'],'GC_YtoX_All')); %nChan x nChan x foi x tvec
        avgbaseAmp = repmat(nanmean(rawMat(:,:,:,basetvecMaskGC),4),[1,1,1,size(rawMat,4)]); % checked, last dimension is the same values
        stdbaseAmp = repmat(nanstd(rawMat(:,:,:,basetvecMaskGC),0,4),[1,1,1,size(rawMat,4)]);%nanstd(X,flag,dim)
        normedGC_YtoX_All_Left = (rawMat - avgbaseAmp)./stdbaseAmp;
        normedGC_YtoX_Avg_Left = squeeze(nanmean(nanmean(normedGC_YtoX_All_Left ,1),2));

        % GC Right video
        rawMat = real(is_load([rootAnalysisDir 'GCAll_Right.mat'],'GC_XtoY_All')); %nChan x nChan x foi x tvec
        avgbaseAmp = repmat(nanmean(rawMat(:,:,:,basetvecMaskGC),4),[1,1,1,size(rawMat,4)]); % checked, last dimension is the same values
        stdbaseAmp = repmat(nanstd(rawMat(:,:,:,basetvecMaskGC),0,4),[1,1,1,size(rawMat,4)]);%nanstd(X,flag,dim)
        normedGC_XtoY_All_Right = (rawMat - avgbaseAmp)./stdbaseAmp;
        normedGC_XtoY_Avg_Right = squeeze(nanmean(nanmean(normedGC_XtoY_All_Right ,1),2));

        rawMat = real(is_load([rootAnalysisDir 'GCAll_Right.mat'],'GC_YtoX_All')); %nChan x nChan x foi x tvec
        avgbaseAmp = repmat(nanmean(rawMat(:,:,:,basetvecMaskGC),4),[1,1,1,size(rawMat,4)]); % checked, last dimension is the same values
        stdbaseAmp = repmat(nanstd(rawMat(:,:,:,basetvecMaskGC),0,4),[1,1,1,size(rawMat,4)]);%nanstd(X,flag,dim)
        normedGC_YtoX_All_Right = (rawMat - avgbaseAmp)./stdbaseAmp;
        normedGC_YtoX_Avg_Right = squeeze(nanmean(nanmean(normedGC_YtoX_All_Right ,1),2));
        

        save([rootAnalysisDir 'plvAllNormed.mat'],'plvAllNormed_Left','plvAllNormed_Right','-v7.3');
        save([rootAnalysisDir 'plvAvgNormed.mat'],'plvAvgNormed_Left','plvAvgNormed_Right','-v7.3');
        save([rootAnalysisDir 'GCAllNormed.mat'],'normedGC_XtoY_All_Left','normedGC_YtoX_All_Left','normedGC_XtoY_All_Right','normedGC_YtoX_All_Right','-v7.3');
        save([rootAnalysisDir 'GCAvgNormed.mat'],'normedGC_XtoY_Avg_Left','normedGC_YtoX_Avg_Left','normedGC_XtoY_Avg_Right','normedGC_YtoX_Avg_Right','-v7.3');
    
        
        % plot normalized PLV and GC
        screensize = get( groot, 'Screensize' );
        fig = figure('Position',[10 50 (screensize(3)-150)*2/3 (screensize(4)-150)*2/3]); %(x,y,width,height) screensize(3)-100

        % plot phase locking value
        subplot(2,3,1)
        imagesc(tvec,1:numel(foi),plvAvgNormed_Left);
        xlabel('Time to event [s]'); ylabel('Left Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-3 3]);
        cl = colorbar('northoutside'); ylabel(cl,['%PLV change: ' regionPairName],'FontSize',12)

        subplot(2,3,2)
        imagesc(tvec,1:numel(foi),normedGC_XtoY_Avg_Left);
        xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);
        cl = colorbar('northoutside'); ylabel(cl,['%GC change: ' regionPairNames1{iRegionPair}],'FontSize',12)

        subplot(2,3,3)
        imagesc(tvec,1:numel(foi),normedGC_YtoX_Avg_Left); 
        xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);
        cl = colorbar('northoutside'); ylabel(cl,['%GC change: ' regionPairNames2{iRegionPair}],'FontSize',12)

        subplot(2,3,4)
        imagesc(tvec,1:numel(foi),plvAvgNormed_Right);
        xlabel('Time to event [s]'); ylabel('Right Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-3 3]);
        cl = colorbar('northoutside'); ylabel(cl,['%PLV change: ' regionPairName],'FontSize',12)

        subplot(2,3,5)
        imagesc(tvec,1:numel(foi),normedGC_XtoY_Avg_Right);
        xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);
        cl = colorbar('northoutside'); ylabel(cl,['%GC change: ' regionPairNames1{iRegionPair}],'FontSize',12)

        subplot(2,3,6)
        imagesc(tvec,1:numel(foi),normedGC_YtoX_Avg_Right);
        xlabel('Time to event [s]'); %ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-5 5]);
        cl = colorbar('northoutside'); ylabel(cl,['%GC change: ' regionPairNames2{iRegionPair}],'FontSize',12)
        
        colormap(jet);       
        
        savefig(fig, [rootAnalysisDir 'normedPLV_GC_bilateral.fig'],'compact');
        saveas(fig, [rootAnalysisDir 'normedPLV_GC_bilateral.png']);
        close all
        clear plvAvgNormed_Left plvAvgNormed_Right normedGC_XtoY_Avg_Left normedGC_YtoX_Avg_Left normedGC_XtoY_Avg_Right normedGC_YtoX_Avg_Right
        clear plvAllNormed_Left plvAllNormed_Right normedGC_XtoY_All_Left normedGC_YtoX_All_Left normedGC_XtoY_All_Right normedGC_YtoX_All_Right

   end
        
end
    

