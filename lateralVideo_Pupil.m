function lateralVideo_Pupil(skipRec, lfpFs,...
    evtTimes, condNames, condID,...
    rootPreprocessDir, rootAnalysisDir)

% if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
%     fprintf('Record %s already analyzed \n',recName'); 
%     if skipRec == 1; continue; end; end     

load([rootPreprocessDir ,'adc_data']);
% if new session, run is_detectSaccadesEphysFun_V3 first; if already ran,
% can just load pupil data

if exist(join([rootPreprocessDir, 'pupil.mat']),'file')
%    load([rootPreprocessDir ,'pupil.mat']); % directly load saved pupil data
else
    [sacTime,cazi,cele,cpup] = is_detectSaccadesEphysFun_V3(rootPreprocessDir);
end
eyeData = [cazi;cele;cpup]; % left eye is first, right eye is second

numConds = numel(condID);
numMaxEvt = 20;
% adc channel names
numPul = size(eyeData,1);
if numPul == 3
    eyeChn = {'rightX','rightY','rightD'};
elseif numPul == 6
    eyeChn = {'leftX','rightX','leftY','rightY','leftD','rightD'};
end

%% parameteres
twin = [-2 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
baseTwin = [-2 0]; %<<<--- baseline to subtract from 
twinSamp = twin*lfpFs;
baseSamp = baseTwin*lfpFs;
numWinsamp = numel(twinSamp(1):twinSamp(2)); 

%% extract snippits of pupil measure

evtDat = nan(numConds,numMaxEvt,numPul,numWinsamp);


for iCond = 1:numel(condID)
    iCondID = condID(iCond);
    evtTime = evtTimes{iCondID};
    numEvt  = numel(evtTime); 
    for iEvt = 1:numMaxEvt
        if iEvt<= length(evtTime) % prevent iEvt exceed evtTime length
            evtSamp = evtTime(iEvt)*lfpFs;
        end
        evtSampWin = round(evtSamp+twinSamp(1):evtSamp+twinSamp(2));
        % normalise to pre-stim pupil
        preBinsWin = round(evtSamp+baseSamp(1):evtSamp+baseSamp(2));
        
        for iChn = 1:numPul
            if evtSampWin(end) <= size(eyeData,2)
                pupMean = nanmedian(eyeData(iChn,preBinsWin),2); %mean
                pupSTD  = nanstd(eyeData(iChn,preBinsWin),[],2);
                evtDat(iCond,iEvt,iChn,:) = (eyeData(iChn, evtSampWin) - pupMean)/pupSTD; %evtSampWin must >0
            end
        end
        
    end  

end


evtDat_trialAvg = reshape(nanmean(evtDat(:,2:end,:,:),2),[numConds,numPul,numWinsamp]); % avg across events
% skip the first trial(evt) since pupil might change dramatically

% don't use squeeze since first dimension might also be 1, but don't want
% to squeeze that
%% plot
timeVec = twin(1):1/lfpFs:twin(2);
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2*numPul/3 (screensize(4)-150)/2]);
legendName = {};
for iChn = 1:size(evtDat_trialAvg,2)
    subplot(3,numPul/3,iChn)
    hold on
    
    for iCond = numConds:-1:1
        iCondID = condID(iCond);
        toPlot = squeeze(evtDat_trialAvg(iCond,iChn,:));
        %toPlot = smoothts(toPlot,'g',3,0.65); % not much difference
        plot(timeVec,toPlot, 'LineWidth', 1)
        legendName{end+1} = condNames{iCondID};
    end
    if iChn ==1; legend(legendName);end
    title([eyeChn{iChn}])
    xlabel('Time [s]');
    ylabel('Normed amplitude');
%     nXtick = round((twin(2)-twin(1))/2)+1;
%     set(gca,'XTick',linspace(1,numWinsamp,nXtick))
%     set(gca,'XTickLabel',linspace(twin(1),twin(2),nXtick))
    vline(0,'k-');vline(5,'k--');
    hline(0,'k-');
    axis tight
    %ylim([-30 30])
end

if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
savefig(fig, [rootAnalysisDir 'Pupil_baselineNormed_mean_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Pupil_baselineNormed_mean_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.png']);

% plot median
evtDat_trialAvg = reshape(nanmedian(evtDat(:,2:end,:,:),2),[numConds,numPul,numWinsamp]); % avg across events
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2*numPul/3 (screensize(4)-150)/2]);
legendName = {};
for iChn = 1:size(evtDat_trialAvg,2)
    subplot(3,numPul/3,iChn)
    hold on
    
    for iCond = numConds:-1:1
        iCondID = condID(iCond);
        toPlot = squeeze(evtDat_trialAvg(iCond,iChn,:));
        %toPlot = smoothts(toPlot,'g',3,0.65); % not much difference
        plot(timeVec,toPlot, 'LineWidth', 1)
        legendName{end+1} = condNames{iCondID};
    end
    if iChn ==1; legend(legendName);end
    title([eyeChn{iChn}])
    xlabel('Time [s]');
    ylabel('Normed amplitude');
%     nXtick = round((twin(2)-twin(1))/2)+1;
%     set(gca,'XTick',linspace(1,numWinsamp,nXtick))
%     set(gca,'XTickLabel',linspace(twin(1),twin(2),nXtick))
    vline(0,'k-');vline(5,'k--');
    hline(0,'k-');
    axis tight
    %ylim([-30 30])
end
savefig(fig, [rootAnalysisDir 'Pupil_baselineNormed_median_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Pupil_baselineNormed_median_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.png']);


% plot contrast
if numConds == 2
    evtDat_trialAvg_contrast = evtDat_trialAvg(1,:,:) - evtDat_trialAvg(2,:,:); %left - right
    evtDat_trialAvg_contrast = reshape(evtDat_trialAvg_contrast,[numPul,numWinsamp]);
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-100)/2*numPul/3 (screensize(4)-150)/2]);
    for iChn = 1:size(evtDat_trialAvg,2)
        subplot(3,numPul/3,iChn)
        if numPul == 3 % right eye
            toPlot = squeeze(-evtDat_trialAvg_contrast(iChn,:)); %right-left
            %toPlot = smoothts(toPlot,'g',3,0.65); % not much difference
            plot(timeVec,toPlot,'k', 'LineWidth', 1) % color has to follow x,y value
        elseif numPul == 6 % left and right eye
            if ismember(iChn,[1,3,5]) %left eye
                toPlot = squeeze(evtDat_trialAvg_contrast(iChn,:)); %left-right
                %toPlot = smoothts(toPlot,'g',3,0.65); % not much difference
                plot(timeVec,toPlot,'k', 'LineWidth', 1)
            elseif ismember(iChn,[2,4,6]) %right eye
                toPlot = squeeze(-evtDat_trialAvg_contrast(iChn,:)); %right-left
                %toPlot = smoothts(toPlot,'g',3,0.65); % not much difference
                plot(timeVec,toPlot,'k', 'LineWidth', 1)
            end
        end

        legend('ipsi-contra video')
        title([eyeChn{iChn}])
        xlabel('Time [s]');
        ylabel('Normed amplitude');
    %     nXtick = round((twin(2)-twin(1))/2)+1;
    %     set(gca,'XTick',linspace(1,numWinsamp,nXtick))
    %     set(gca,'XTickLabel',linspace(twin(1),twin(2),nXtick))
        vline(0,'k-');vline(5,'k--');
        axis tight
        %ylim([-30 30])
    end
    savefig(fig, [rootAnalysisDir 'Pupil_contrast_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Pupil_contrast_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.png']);
    save([rootAnalysisDir 'Pupil_contrast'], 'timeVec','evtDat_trialAvg_contrast');
end
save([rootAnalysisDir 'Pupil_baselineNormed_median'],'timeVec','evtDat_trialAvg');
close all
end