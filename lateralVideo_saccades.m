function lateralVideo_saccades(evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
    validChn,rootPreprocessDir, rootAnalysisDir, animalCode)

doVisualUnit = 0;

% load and process ttl data in ephys
load([rootPreprocessDir 'saccades.mat']);
totalNumChn  = length(sacTime); % get total number of channels before exclusion
numConds     = numel(condID);

binSize = 0.5; %.020;    % in seconds
twin = round([-3 8]); %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
baseTwin = [-3 0]; %<<<--- baseline to subtract from spectrogram 
stimTwin = [0 5];
numDivX = twin(2)-twin(1)+1;
numTvec = numel(twin(1):binSize:twin(2))-1; %to match the lenth of time PSTH
regionNames = {'lEye','rEye'};
%condID = [1,3,4]; 
% calculate PSTH
condCount = 1;
frZ = nan(numConds,totalNumChn,numTvec);
for iCond = 1:numConds
    iCondID = condID(iCond);
    evtTime = evtTimes{iCondID};
    for iChn = 1:totalNumChn
        [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,sacTime{iChn},twin,binSize); % CZ: PSTHrate's 1st timept is time of saccade
        % PSTH(condCount,iChn,:)  = PSTHrate; %if want to save raw PSTHrate
        
        if doVisualUnit == 1
            stimtvecMask = timePSTH>= stimTwin(1) & timePSTH <= stimTwin(2);
            basetvecMask = timePSTH>= baseTwin(1) & timePSTH <= baseTwin(2);
            trialStimFR = nanmean(psthTrial(:,stimtvecMask),2);
            trialBaseFR = nanmean(psthTrial(:,basetvecMask),2);
            [h,p,ci,stats] = ttest2(trialStimFR,trialBaseFR); % test if stim duration has diff FR than baseline
            %if p > 0.1;continue;end %if not stim modulated, skip this
            %channal -- too strict, few channels left
            if mean(trialStimFR) < 1.2* mean(trialBaseFR); continue;end
        end
        
        % normalise to pre saccade firing rate
        preBins = (timePSTH>=baseTwin(1) & timePSTH<=baseTwin(2)); % -2~0sec before stimulus
        frMean  = mean(PSTHrate(preBins));
        frSTD   = std(PSTHrate(preBins));
        frZ(condCount,iChn,:) = (PSTHrate-frMean)/frSTD; % Spike z score
        %spkCell{condCount,iChn} = spks; % save spike times in sec for later
    end
    
    numBins = numel(timePSTH); % always the same, only need to save once
    
    condCount = condCount + 1;
    
end


%% 
data2Analyze = frZ;

%% chn avged
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)*2/3 (screensize(4)-150)/4]);

for iRegion = 1:totalNumChn % last 2 regions are EEG
    subplot(1,totalNumChn,iRegion)
    hold on
    legendName = {};
    if length(condID) == 1
        % delete all channels with all zeros
        sliceData = reshape(data2Analyze(iCond,iRegion,:),[1,size(data2Analyze,3)]);
        %data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
        toPlot(iRegion,:) = squeeze(nanmean(sliceData,1));
        toPlot(iRegion,:) = smoothts(toPlot(iRegion,:),'g',3,0.65); % smooth the line
        plot(timePSTH, toPlot(iRegion,:), 'LineWidth', 4)        
        legend(condNames{condID(1)});
        xlim(twin);
    else
        for iCond = numConds:-1:1 % control in the back
            sliceData = reshape(data2Analyze(iCond,iRegion,:),[1,size(data2Analyze,3)]);
            %data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
            toPlot(iCond,iRegion,:) = squeeze(nanmean(sliceData,1));
            toPlot(iCond,iRegion,:) = smoothts(toPlot(iCond,iRegion,:),'g',3,0.65); % smooth the line
            plot(timePSTH, squeeze(toPlot(iCond,iRegion,:)), 'LineWidth', 4);
            xlim(twin);
            legendName{end+1} = condNames{condID(iCond)};
    
        end          
    end
    
    if iRegion == 1; legend(legendName); end% NOTE: match plotting order
    title(['Z-score saccades: ' regionNames{iRegion} ])
    xlabel('Time [s]');
    ylabel('Z-score saccades');
    %axis tight
    ylim([-3 3]);
    vline(0,'k-');vline(5,'k--');
    hline(0,'k-');
    

end

if doVisualUnit == 1
    save([rootAnalysisDir 'saccade_mean_visual.mat'],'timePSTH','toPlot','frZ','regionNames', '-v7.3');
    savefig(fig, [rootAnalysisDir 'Z-score saccade_chn-avg_visual0-2.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score saccade_chn-avg_visual0-2.png']);
else
    save([rootAnalysisDir 'saccade_mean.mat'],'timePSTH','toPlot','frZ','regionNames', '-v7.3');
    savefig(fig, [rootAnalysisDir 'Z-score saccade_chn-avg.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Z-score saccade_chn-avg.png']);
end

close all % close figures
end
