function lateralVideo_PSTH(evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
    validChn,rootPreprocessDir, rootAnalysisDir, animalCode)

% AH 20190510: added combineCondition so that units were detected for both conditions and then plotted
% use smoothdata gaussian 30 to smooth the line
doVisualUnit = 1;

% load and process ttl data in ephys
load([rootPreprocessDir 'validChn']);
files = dir([rootPreprocessDir 'spikes\spk*.mat']);
totalNumChn  = length(files); % get total number of channels before exclusion
numConds     = numel(condID);
switch animalCode
    case '0172'
        numSpkRegion = numel(regionNames);
    case '0168'
        numSpkRegion = numel(regionNames);
    case '0169'
        numSpkRegion = numel(regionNames) - 2; % last two are EEG
end
binSize = 0.020;    % in seconds
twin = round([-3 8]); %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
baseTwin = [-2 -1]; %<<<--- baseline to subtract from spectrogram 
stimTwin = [0 2];
numDivX = twin(2)-twin(1)+1;
numTvec = numel(twin(1):binSize:twin(2))-1; %to match the lenth of time PSTH

% calculate PSTH
condCount = 1;
frZ = nan(numConds,totalNumChn,numTvec);
visualUnits = []; % to populate the channels we want to keep

for iChn = 1:totalNumChn
    load([rootPreprocessDir 'spikes\spk_' num2str(iChn)]);    
    spks  = spkTime;                
    for iCond = 1:numConds
        iCondID = condID(iCond);
        evtTime = evtTimes{iCondID};       
        % spike times in seconds % OLD: ./1000; % convert spike times (ms) to seconds
        [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,spks,twin,binSize); % CZ: PSTHrate's 1st timept is time of saccade
        % PSTH(condCount,iChn,:)  = PSTHrate; %if want to save raw PSTHrate
        stimtvecMask = timePSTH>= stimTwin(1) & timePSTH <= stimTwin(2);
        basetvecMask = timePSTH>= baseTwin(1) & timePSTH <= baseTwin(2);
        trialStimFR = nanmean(psthTrial(:,stimtvecMask),2);
        trialBaseFR = nanmean(psthTrial(:,basetvecMask),2);
        
        % normalise to pre saccade firing rate
        preBins = (timePSTH>baseTwin(1) & timePSTH<baseTwin(2)); % -2~0sec before stimulus
        frMean  = mean(PSTHrate(preBins));
        frSTD   = std(PSTHrate(preBins));
        frZ(iCondID,iChn,:) = (PSTHrate-frMean)/frSTD; % Spike z score
        %spkCell{iCondID,iChn} = spks; % save spike times in sec for later
        
        if doVisualUnit == 1
            if mean(trialStimFR) >= 1.1* mean(trialBaseFR)
                visualUnits = [visualUnits, iChn];
            end
        end
    end
    
    numBins = numel(timePSTH); % always the same, only need to save once
    
    %condCount = condCount + 1;
    
end
visualUnits = unique(visualUnits);
%visualValidUnits = intersect(visualUnits, [validChn{:}]); % combine visual Units with validChns
data2Analyze = frZ;
%% 
if doVisualUnit == 1 
    for iRegion = 1:numel(validChn)
        validChn{iRegion} = intersect(validChn{iRegion}, visualUnits);
    end
    %enddata2Analyze = frZ(:,visualValidUnits,:);
end


ipanel = 1;

screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)*2/3 (screensize(4)-150)/4]);

for iCond = 1:numConds
    iCondID = condID(iCond);
    for iRegion = 1:numSpkRegion
        
        toPlot = reshape(data2Analyze(iCondID,validChn{iRegion},:),[],size(data2Analyze,3));
        
        subplot(numConds,numSpkRegion,ipanel)
        
        imagesc(toPlot)
        title(['Z-score FR PSTH: ' regionNames{iRegion} '; ' condNames{iCondID}])
        xlabel('Time [s]');
        ylabel('Channel');
        set(gca,'XTick',linspace(1,numBins,numDivX))
        set(gca,'XTickLabel',linspace(twin(1),twin(2),numDivX))
        h = colorbar;
        ylabel(h, 'Z-score FR')
        caxis([-0.5 6])
        axis tight
        
        ipanel = ipanel + 1;
    end
    
end
clear toPlot
if ~exist(rootAnalysisDir,'dir'); mkdir(rootAnalysisDir); end

if doVisualUnit == 1
    savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_visual0-2.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_visual0-2.png']);
else
    savefig(fig, [rootAnalysisDir 'Z-score FR PSTH.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score FR PSTH.png']);
end

%% chn avged
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)*2/3 (screensize(4)-150)/4]);

for iRegion = 1:numSpkRegion % last 2 regions are EEG
    subplot(1,numSpkRegion,iRegion)
    hold on
    legendName = {};
    if length(condID) == 1
        iCondID = condID(1);
        % delete all channels with all zeros
        sliceData = reshape(data2Analyze(iCondID,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
        data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
        toPlot(iRegion,:) = squeeze(nanmean(data2Average,1));
        %toPlot(iRegion,:) = smoothts(toPlot(iRegion,:),'g',3,0.65); % smooth the line
        toPlot(iRegion,:) = smoothdata(toPlot(iRegion,:),'gaussian',30);
        plot(timePSTH, toPlot(iRegion,:), 'LineWidth', 2)        
        legend(condNames{condID(1)});
        xlim(twin);
    else
        for iCond = numConds:-1:1 % control in the back
            iCondID = condID(iCond);
            sliceData = reshape(data2Analyze(iCondID,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
            data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
            
            toPlot(iCondID,iRegion,:) = squeeze(nanmean(data2Average,1));
            %toPlot(iCondID,iRegion,:) = smoothdata(toPlot(iCondID,iRegion,:),'g',10,0.65); % smooth the line
            toPlot(iCondID,iRegion,:) = smoothdata(toPlot(iCondID,iRegion,:),'gaussian',30); % window, length
            plot(timePSTH, squeeze(toPlot(iCondID,iRegion,:)), 'LineWidth', 2)
            xlim(twin);
            legendName{end+1} = condNames{iCondID};
        end          
    end
    
    if iRegion == 1; legend(legendName); end % NOTE: match plotting order
    title(['Z-score FR PSTH: ' regionNames{iRegion} ])
    xlabel('Time [s]');
    ylabel('Z-score Firing Rate [Hz]');
    
    
    axis tight
    ylim([-1 2])
end

if doVisualUnit == 1
    save([rootAnalysisDir 'zPSTH_mean_visual.mat'],'timePSTH','toPlot','frZ','visualUnits','validChn', '-v7.3');
    savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg_visual0-2.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg_visual0-2.png']);
else
    save([rootAnalysisDir 'zPSTH_mean.mat'],'timePSTH','toPlot','frZ','validChn', '-v7.3');
    savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg.png']);
end


%% Median
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)*2/3 (screensize(4)-150)/4]);

for iRegion = 1:numSpkRegion % last 2 regions are EEG
    subplot(1,numSpkRegion,iRegion)
    hold on
    
    if length(condID) == 1
        iCondID = condID(1);
        % delete all channels with all zeros
        sliceData = reshape(data2Analyze(iCond,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
        data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
        toPlot(iRegion,:) = squeeze(nanmedian(data2Average,1));
        %toPlot(iRegion,:) = smoothts(toPlot(iRegion,:),'g',3,0.65); % smooth the line
        toPlot(iRegion,:) = smoothdata(toPlot(iRegion,:),'gaussian',30);
        plot(timePSTH, toPlot(iRegion,:), 'LineWidth', 2)        
        legend(condNames{condID(1)});
        xlim(twin);
    else
        for iCond = numConds:-1:1 % control in the back
            iCondID = condID(iCond);
            sliceData = reshape(data2Analyze(iCondID,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
            data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
            toPlot(iCondID,iRegion,:) = squeeze(nanmedian(data2Average,1));
            %toPlot(iCondID,iRegion,:) = smoothts(toPlot(iCondID,iRegion,:),'g',3,0.65); % smooth the line
            toPlot(iCondID,iRegion,:) = smoothdata(toPlot(iCondID,iRegion,:),'gaussian',30); % window, length
            plot(timePSTH, squeeze(toPlot(iCondID,iRegion,:)), 'LineWidth', 2)
            xlim(twin);
            legendName{end+1} = condNames{iCondID};
        end        
    end
    
    if iRegion == 1; legend(legendName); end
    title(['Z-score FR PSTH: ' regionNames{iRegion} ])
    xlabel('Time [s]');
    ylabel('Z-score Firing Rate [Hz]');    
    
    axis tight
    ylim([-1 2])
end

if doVisualUnit == 1
    save([rootAnalysisDir 'zPSTH_median_visual.mat'],'timePSTH','toPlot','frZ','visualUnits','validChn', '-v7.3');
    savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-median_visual0-2.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-median_visual0-2.png']);
else
    save([rootAnalysisDir 'zPSTH_median.mat'],'timePSTH','toPlot','frZ','validChn', '-v7.3');
    savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-median.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-median.png']);
end
close all % close figures
end
