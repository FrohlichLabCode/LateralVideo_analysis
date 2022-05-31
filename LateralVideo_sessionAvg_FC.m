clear
clc

animalCode = '0172';
skipRec = 1;

addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/sessionFC/'];

medial  = [40,42,43,45,46,48,49,51,52,54,55,57,59,60,62,63,65,66,68,70]; %20
lateral = [1,2,4,6,8,10,12,14,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39]; %24
half    = [3,5,7,9,11,13,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,58,61,64,69]; %25

sessionTypes =  {medial, lateral, half,};
sessionTypeNames = {'medial','lateral', 'half'};
GroupSessionDir = [GroupAnalysisDir];

condNames   = {'lVideo', 'rVideo'};
regionNames = {'FC','LPl','PPC','VC'};
regionPairs = {[1,3],[2,3],[2,4],[3,4]};
numConds    = numel(condNames);
numRegions  = numel(regionNames);
numPairs    = numel(regionPairs);

% generate name pair cells
for iPair = 1:numel(regionPairs)
    regionPairNames{iPair} = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}];
    regionPair_Names{iPair} = [regionNames{regionPairs{iPair}(1)} '_' regionNames{regionPairs{iPair}(2)}];%{'FC_PPC', 'LPl_PPC', 'LPl_VC', 'PPC_VC'};
    regionPairNamesGC{2*iPair-1} = [regionNames{regionPairs{iPair}(1)} '->' regionNames{regionPairs{iPair}(2)}];%{'LPl->PPC','PPC->LPl';'LPl->VC','VC->LPl'};
    regionPairNamesGC{2*iPair} = [regionNames{regionPairs{iPair}(2)} '->' regionNames{regionPairs{iPair}(1)}];
end

for i = 1%1:numel(sessionTypes)
    sessionIDs = sessionTypes{i};
    sessionType = sessionTypeNames{i};
    irec = 0;
    if skipRec ==1 && exist([GroupAnalysisDir 'allSessionFC_' sessionType '.mat'])
        load([GroupAnalysisDir 'allSessionFC_' sessionType '.mat'])
    else  

    for iSession = 1:numel(sessionIDs)
        sessionID = sessionIDs(iSession);
        sessionName = sprintf('%03d',sessionID); % add leading zeros        
        fileInfo = dir([AnalysisDir animalCode '_LateralVideo_' sessionName '*']); % detect files to load/convert  '_LateralVideo*'
        irec = irec +1;
        
    % combine data from all sessions
    %for irec = 1:numel(fileInfo)
    recName = fileInfo.name;
    splitName   = strsplit(recName,'_');
    fprintf('processing %d   %s \n', irec, recName);
    date = splitName{4};
    for iRegionPair = 1:numel(regionPairNames)
        regionPairName = regionPairNames{iRegionPair};
        regionPair_Name = regionPair_Names{iRegionPair};
        rootAnalysisDir = [AnalysisDir recName '/FC_validChns/' regionPairName '/'];
        if ~exist([rootAnalysisDir 'funcCon_median_rVideo.mat'])
            fprintf('%s not existed \n', regionPairName) 
            continue;end
    for iCond = 1:numConds
        condName = condNames{iCond};
        %load([rootAnalysisDir 'plvAll_' condName '.mat']);
        %try

        if strcmp(regionPairName, 'FC-PPC')                
            Spec.FC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgXSpec');
            SpecNorm.FC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'LPl-PPC')
            Spec.LPl(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgXSpec');
            SpecNorm.LPl(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'PPC-VC')
            Spec.PPC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgXSpec');
            SpecNorm.PPC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgXNormed');
            Spec.VC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgYSpec');
            SpecNorm.VC(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgYNormed');
        end

        PLV.(regionPair_Name)(irec,iCond,:,:) = is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgPLV');
        ICoherence.(regionPair_Name)(irec,iCond,:,:) = abs(is_load([rootAnalysisDir 'funcCon_median_' condName '.mat'], 'avgImagZ'));
        %catch
        %end
        %try

        GC.(regionPair_Name)(irec,iCond,1,:,:) = is_load([rootAnalysisDir 'GC_medain_' condName '.mat'], 'avgGC_XtoY');
        GC.(regionPair_Name)(irec,iCond,2,:,:) = is_load([rootAnalysisDir 'GC_medain_' condName '.mat'], 'avgGC_YtoX');
        if irec==1; GC.tvec = is_load([rootAnalysisDir 'GC_medain_' condName '.mat'],'tvecGC');
            [foi tvec] = is_load([rootAnalysisDir 'specAll_' condName '.mat'],'foi','tvec');end        

        %catch
        %end
    end
    end
    
    tvecGC = GC.tvec;
    [foi, tickLoc, tickLabel] = getFoiLabel(2, 128, 150, 2);
        if ~exist(GroupAnalysisDir,'dir'); mkdir(GroupAnalysisDir);end
        save([GroupAnalysisDir 'allSessionFC_' sessionType '.mat'],'tvec','foi','tvecGC', 'tickLabel','tickLoc','regionPairNames','regionPair_Names','regionPairNamesGC', 'Spec','SpecNorm','PLV','ICoherence','GC','-v7.3');
    end
    end
    
    
    
    
    
    % clean line noise
    irec = 0;
    for iSession = 1:numel(sessionIDs)
        sessionID = sessionIDs(iSession);
        sessionName = sprintf('%03d',sessionID); % add leading zeros        
        fileInfo = dir([AnalysisDir animalCode '_LateralVideo_' sessionName '*']); % detect files to load/convert  '_LateralVideo*'
        irec = irec +1;
        
        recName = fileInfo.name;
        splitName   = strsplit(recName,'_');
        date = splitName{4};
        if (datetime(date, 'InputFormat', 'yyyyMMdd') >= datetime('20190111','InputFormat', 'yyyyMMdd')) && ~strcmp(date,'20190115')
            fprintf('cleaning %s \n', recName);
            for iCond = 1:numConds
                condName = condNames{iCond};
                cutTwin = [90,132];
                cutMask = [cutTwin(1):cutTwin(2)];
                Spec.PPC(irec,iCond,cutMask,:) = arrayLinspace(Spec.PPC(irec,iCond,cutTwin(1),:), Spec.PPC(irec,iCond,cutTwin(2),:), numel(cutMask));
                a = squeeze(Spec.PPC(:,1,:,1));
                SpecNorm.PPC(irec,iCond,cutMask,:) = arrayLinspace(SpecNorm.PPC(irec,iCond,cutTwin(1),:), SpecNorm.PPC(irec,iCond,cutTwin(2),:), numel(cutMask));
                Spec.VC(irec,iCond,cutMask,:) = arrayLinspace(Spec.VC(irec,iCond,cutTwin(1),:), Spec.VC(irec,iCond,cutTwin(2),:), numel(cutMask));
                SpecNorm.VC(irec,iCond,cutMask,:) = arrayLinspace(SpecNorm.VC(irec,iCond,cutTwin(1),:), SpecNorm.VC(irec,iCond,cutTwin(2),:), numel(cutMask));
                for iRegionPair = 1:numel(regionPairNames)
                    regionPair_Name = regionPair_Names{iRegionPair};
                    PLV.(regionPair_Name)(irec,iCond,cutMask,:) = arrayLinspace(PLV.(regionPair_Name)(irec,iCond,cutTwin(1),:), PLV.(regionPair_Name)(irec,iCond,cutTwin(2),:), numel(cutMask))*0.3;
                    %ICoherence.(regionPair_Name)(irec,iCond,cutMask,:) = arrayLinspace(ICoherence.(regionPair_Name)(irec,iCond,cutTwin(1),:), ICoherence.(regionPair_Name)(irec,iCond,cutTwin(2),:), numel(cutMask));
                    GC.(regionPair_Name)(irec,iCond,1,cutMask,:) = arrayLinspace(GC.(regionPair_Name)(irec,iCond,1,cutTwin(1),:), GC.(regionPair_Name)(irec,iCond,1,cutTwin(2),:), numel(cutMask))*0.3;
                    GC.(regionPair_Name)(irec,iCond,2,cutMask,:) = arrayLinspace(GC.(regionPair_Name)(irec,iCond,2,cutTwin(1),:), GC.(regionPair_Name)(irec,iCond,2,cutTwin(2),:), numel(cutMask))*0.3;
                end
            end
        end
    end
    save([GroupAnalysisDir 'allSessionFC_cleaned_' sessionType '.mat'],'tvec','foi','tvecGC', 'tickLabel','tickLoc','regionPairNames','regionPair_Names','regionPairNamesGC', 'Spec','SpecNorm','PLV','ICoherence','GC','-v7.3');
    clear Spec SpecNorm PLV ICoherence GC
end
    




%% plot median across sessions
for i = 1%:numel(sessionTypes)
    sessionIDs = sessionTypes{i};
    sessionType = sessionTypeNames{i};
    load([GroupAnalysisDir 'allSessionFC_cleaned_' sessionType '.mat']);
    
    
% plot Spec for all regions
fig = figure('name','medianSpec','position', [10 20 320*numConds 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    xLabel = 'Time to init [sec]';
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numRegions,numConds,(iRegion-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),pow2db(squeeze(nanmedian(Spec.(regionName)(:,iCond,:,:),1))));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([15 45]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianSpec_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianSpec_' sessionType '.png']);
% 
fig = figure('name','medianSpecNorm','position', [10 20 320*numConds 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    xLabel = 'Time to init [sec]';
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numRegions,numConds,(iRegion-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,iCond,:,:),1))));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);
        %if strcmp(regionName, 'FC');caxis([-2 4]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianSpecNorm_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianSpecNorm_' sessionType '.png']);

% plot PLV for 4 pairs
fig = figure('name','medianPLV','position', [10 20 320*numConds 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    xLabel = 'Time to init [sec]';
    % plot PLV
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numRegions,numConds,(iRegionPair-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),squeeze(nanmedian(PLV.(regionPair_Name)(:,iCond,:,:),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.35 0.75]);
        if strcmp(regionPair_Name, 'LPl_PPC');caxis([0.45 0.9]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianPLV_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianPLV_' sessionType '.png']);

% plot Icoherence for 4 pairs
fig = figure('name','medianICoherence','position', [10 20 320*numConds 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    xLabel = 'Time to init [sec]';
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        subplot(numPairs,numConds,(iRegionPair-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),squeeze(nanmedian(ICoherence.(regionPair_Name)(:,iCond,:,:),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 3.5]);
        %if strcmp(regionPair_Name, 'LPl_PPC');caxis([0.05 0.6]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianICoherence_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianICoherence_' sessionType '.png']);


fig = figure('name','medianGC','position', [10 20 320*numConds*2 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    xLabel = 'Time to init [sec]';
    for iCond = 1:numConds %column
        condName = condNames{iCond};
        % X -> Y
        subplot(numPairs,numConds*2,(2*iRegionPair-2)*numConds+iCond)
        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(regionPair_Name)(:,iCond,1,:,:)),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ':' condName],'FontSize',12)
        if iRegionPair == 2; caxis([0 0.15]);else caxis([0,0.06]);end 
        ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values
    
        % Y -> X
        subplot(numPairs,numConds*2,(2*iRegionPair-1)*numConds+iCond)
        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(regionPair_Name)(:,iCond,2,:,:)),1)));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair} ':' condName],'FontSize',12)
        if iRegionPair == 2; caxis([0 0.15]);else caxis([0,0.06]);end
        ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values   
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianGC_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianGC_' sessionType '.png']);

close all

%% plot contrast
% rVideo vs. lVideo
fig = figure('name','medianSpec','position', [10 20 320*(numConds-1) 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    xLabel = 'Time to stimOn [sec]';
    for iCond = 1:(numConds-1) %column
        condName = [condNames{iCond+1}(1) '-' condNames{1}];
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
        subplot(numRegions,numConds-1,(iRegion-1)*(numConds-1)+iCond)
        contrast = pow2db(squeeze(nanmedian(Spec.(regionName)(:,iCond+1,:,:),1))) - pow2db(squeeze(nanmedian(Spec.(regionName)(:,1,:,:),1)));
        imagesc(tvec,1:numel(foi),contrast);
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)      
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianSpec_r-lVideo_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianSpec_r-lVideo_' sessionType '.png']);


fig = figure('name','medianPLV','position', [10 20 320*(numConds-1) 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionName = regionPair_Names{iRegionPair};
    xLabel = 'Time to stimOn [sec]';
    for iCond = 1:(numConds-1) %column
        condName = [condNames{iCond+1}(1) '-' condNames{1}];
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
         subplot(numPairs,numConds-1,(iRegionPair-1)*(numConds-1)+iCond)
        contrast = squeeze(nanmedian(PLV.(regionName)(:,iCond+1,:,:),1)) - squeeze(nanmedian(PLV.(regionName)(:,1,:,:),1));
        imagesc(tvec,1:numel(foi),contrast);
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-0.2 0.2]);
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)      
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianPLV_r-lVideo_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianPLV_r-lVideo_' sessionType '.png']);


fig = figure('name','medianICoherence','position', [10 20 320*(numConds-1) 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionName = regionPair_Names{iRegionPair};
    xLabel = 'Time to stimOn [sec]';
    for iCond = 1:(numConds-1) %column
        condName = [condNames{iCond+1}(1) '-' condNames{1}];
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
         subplot(numPairs,numConds-1,(iRegionPair-1)*(numConds-1)+iCond)
        contrast = squeeze(nanmedian(ICoherence.(regionName)(:,iCond+1,:,:),1)) - squeeze(nanmedian(ICoherence.(regionName)(:,1,:,:),1));
        imagesc(tvec,1:numel(foi),contrast);
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-1 1]);
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)      
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianICoherence_r-lVideo_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianICoherence_r-lVideo_' sessionType '.png']);


fig = figure('name','medianGC','position', [10 20 320*(numConds-1)*2 270*numRegions]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionName = regionPair_Names{iRegionPair};
    xLabel = 'Time to stimOn [sec]';
    for iCond = 1:(numConds-1) %column
        condName = [condNames{iCond+1}(1) '-' condNames{1}];
        
        subplot(numPairs,(numConds-1)*2,(iRegionPair-1)*numConds+iCond)
        contrast = squeeze(nanmedian(real(GC.(regionName)(:,iCond+1,1,:,:)),1))- squeeze(nanmedian(real(GC.(regionName)(:,1,1,:,:)),1));
        imagesc(GC.tvec,1:numel(foi),contrast);
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ':' condName],'FontSize',12)
        caxis([-0.02 0.02]);
        
        subplot(numPairs,(numConds-1)*2,(iRegionPair-1)*numConds+iCond+1)
        contrast = squeeze(nanmedian(real(GC.(regionName)(:,iCond+1,2,:,:)),1))- squeeze(nanmedian(real(GC.(regionName)(:,1,2,:,:)),1));
        imagesc(GC.tvec,1:numel(foi),contrast);
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair} ':' condName],'FontSize',12)
        caxis([-0.02 0.02]);
     end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'medianGC_r-lVideo_' sessionType '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'medianGC_r-lVideo_' sessionType '.png']);


clear Spec SpecNorm PLV ICoherence GC
end