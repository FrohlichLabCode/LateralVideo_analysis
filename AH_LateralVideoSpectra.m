function AH_LateralVideoSpectra(recPath,recName,animalCode)
% modified from Iain's code is_imagesVideoSpectra.m
% AH 2018.6.29

linORlog = 2; %freqs of interest: 1=linear 2=log

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if linORlog == 1
    numFreqs = 100;
    lowFreq  = 1;
    highFreq = 80;
    foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
    numFreqs = 100;
    lowFreq  = 2;
    highFreq = 128;
    foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end

% region info
switch animalCode
    case '0168'
    regionNames = {'LPl','PPC','VC'};
    numRegion   = numel(regionNames);
    %regionPairs = {[1,2],[1,3],[2,3]};
    noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd'))...
                  +(datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') == datetime('20180718', 'InputFormat', 'yyyyMMdd'));
    
    case '0169' % 0169 region info
    regionNames = {'lLPl','lPPC','rPPC','rLPl','lVCEEG','rVCEEG'};
    numRegion   = numel(regionNames);
    %regionPairs = {[1,2],[1,5],[2,5],[4,3],[4,6],[3,6],[1,4],[2,3],[5,6]};
    noPhotoData = datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd');
end

%% select LFP to load
if MedianorPCA == 0
    lfpMat = lfp.reorderedSig;
    for i = 1:numel(lfp.validChn) 
        regionChn{i} = lfp.validChn{i}; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(lfp.reorderedChn{i},:); % reordered channel correspond to reordered lfp
    end
    
elseif MedianorPCA == 1
    lfpMat = lfp.median;
    for i = 1:size(lfp.median,1) %lfp.median is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end

elseif MedianorPCA == 2
    lfpMat = lfp.PCA;
    for i = 1:size(lfp.PCA,1) %lfp.PCA is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end
end

%% parameteres
lfpFs  = lfp.Fs;
twin = [-3 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
baseTwin = [-2.5 -1]; %<<<--- baseline to subtract from spectrogram 
% define window
win = round([-3 8]*lfpFs);% 5s movie + 3s gray screen
tvec = (win(1):win(2))/1000; % time vector in sec


%% preprocess session log file
% %to deal with more complicated file names:
% logPath = 'C:/Users/angel/Dropbox (Frohlich Lab)/Codebase/CodeAngel/PresentationScripts/Logs/';
% F = sprintf('%s*%s*.log',animalCode, recName(end-7:end)); %eg.0168*20180627*.log
% S = dir(fullfile(logPath,F));
% fileName = S.name;
% fileID   = fopen([logPath fileName]);

fileID  = fopen([rootBehavDatDir '.log']);
formatSpec = '%f %s %s %f %f %f %f %f %f %f %s %f';
LogFile = textscan(fileID,formatSpec,'HeaderLines',5,'Delimiter', '\t');
fclose(fileID);

% detect location
Col_TrialNumber = 1;
Col_TrialType = 2; % video
Col_Code = 3;      % trial side and video name

cellLeft = cellfun(@regexp,LogFile(Col_Code),{'1;'},'UniformOutput',false);
leftInd = cellfun(@isempty,cellLeft{:}) == 0;
cellRight = cellfun(@regexp,LogFile(Col_Code),{'2;'},'UniformOutput',false);
rightInd = cellfun(@isempty,cellRight{:}) == 0;

% detect stim onsets
ttlInd = 1;
rawFs = 30000;
TTLonset  = find(diff(triggerData(ttlInd,:))==1)./rawFs; % in sec
TTLoffset = find(diff(triggerData(ttlInd,:))==-1)./rawFs;

if usePhotodiode == 1 && noPhotoData == 0 % not a trial without photo data
    load([rootPreprocessDir 'adc_data']);
    photoInd = 4;
    photoRaw = adc_data(photoInd,:);
    temp = zeros(1,length(photoRaw));
    temp(photoRaw>3.2) = 1;
    photoOnset = [];
    for i = TTLonset*rawFs % for each trigger time point    
        photoOnset = [photoOnset, (i+find(diff(temp(i:i+rawFs*0.1))==1,1,'last'))./rawFs]; % the last drop->0 is when gray->black in photodiode
    end
    onset = photoOnset;
else
    onset = TTLonset + 0.05; % on average, trigger data is 0.0529 sec earlier than visual stimuli onset
end

% % reject noise in recording (high amplitude noise can destroy coherence estimates in particular)
% for i = 1:numel(size(lfpMat,1)); lfpDenoise(i,:) = sub_rejectNoise(lfpMat(i,:),fs,2,1); 
% xsig = sub_rejectNoise(xser,fs,2,1); 
% ysig = sub_rejectNoise(yser,fs,2,1); 

condNames = {'Left','Right'};
condID = [1 2];
numConds = numel(condID);

for iCond = condID

    if iCond == 1
        evtTime = onset(leftInd); % i.e. leftOnset
    elseif iCond == 2
        evtTime = onset(rightInd); % i.e. rightOnset
    end
    evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window
    
    evtTimes{iCond} = evtTime;
end

% initialise matrices
powLvideo_pul = nan(length(foi),numel(regionChn{1}),diff(win)+1); 
powLvideo_ppc = nan(length(foi),numel(regionChn{2}),diff(win)+1); 
powLvideo_vc  = nan(length(foi),numel(regionChn{3}),diff(win)+1);
powRvideo_pul = nan(length(foi),numel(regionChn{1}),diff(win)+1);
powRvideo_ppc = nan(length(foi),numel(regionChn{2}),diff(win)+1);
powRvideo_vc  = nan(length(foi),numel(regionChn{3}),diff(win)+1);

% get wavelets
morWav = is_makeWavelet(foi,lfpFs);

for f = 1:numel(foi)
    display([recName ' ' num2str(f) '/' num2str(numFreqs)])
    % Pulvinar first
    C_pul = conv2(lfpMat(regionChn{1},:),morWav{f},'same'); % convole LFP with wavelet, matrix of complex numbers
    zpow  = zscore(abs(C_pul).^2,[],2); % z-score power over channels (If dim = 2, then zscore uses the means and standard deviations along the rows of X.)
    angPul = angle(C_pul); clear C_pul
    
    % LeftVideo Pul
    smat = nan(length(evtTimes{1}),size(zpow,1),diff(win)+1);
    for iev = 1:numel(evtTimes{1})
        evSamp = round(evtTimes{1}(iev)*lfpFs);
        smat(iev,:,:) = zpow(:,evSamp+win(1):evSamp+win(2)); end
    powLvideo_pul(f,:,:) = squeeze(nanmean(smat,1));
    
    % RightVideo Pul
    smat = nan(length(rightOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(rightOnset); smat(iev,:,:) = zpow(:,rightOnset(iev)+win(1):rightOnset(iev)+win(2)); end
    powRvideo_pul(f,:,:) = squeeze(nanmean(smat,1));
    
    % PPC second
    C_ppc = conv2(lfpMat(ppcChans,:),morWav{f},'same'); % convole PPC LFP with wavelet
    zpow  = zscore(abs(C_ppc).^2,[],2); % z-score power
    angPPC = angle(C_ppc); clear C_ppc
    
    % LeftVideo PPC
    smat = nan(length(leftOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(leftOnset); smat(iev,:,:) = zpow(:,leftOnset(iev)+win(1):leftOnset(iev)+win(2)); end
    powLvideo_ppc(f,:,:) = squeeze(nanmean(smat,1));
    
    % RightVideo PPC
    smat = nan(length(rightOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(rightOnset); smat(iev,:,:) = zpow(:,rightOnset(iev)+win(1):rightOnset(iev)+win(2)); end
    powRvideo_ppc(f,:,:) = squeeze(nanmean(smat,1));

    % VC third
    C_vc = conv2(lfpMat(vcChans,:),morWav{f},'same'); % convole PPC LFP with wavelet
    zpow  = zscore(abs(C_vc).^2,[],2); % z-score power
    angVC = angle(C_vc); clear C_vc
    
    % LeftVideo VC
    smat = nan(length(leftOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(leftOnset); smat(iev,:,:) = zpow(:,leftOnset(iev)+win(1):leftOnset(iev)+win(2)); end
    powLvideo_vc(f,:,:) = squeeze(nanmean(smat,1));
    
    % RightVideo VC
    smat = nan(length(rightOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(rightOnset); smat(iev,:,:) = zpow(:,rightOnset(iev)+win(1):rightOnset(iev)+win(2)); end
    powRvideo_vc(f,:,:) = squeeze(nanmean(smat,1));

    
    clear zpow
    
%     % compute PLV
%     plvWin = round([-0.5 0.5]*lfpFs);
%     steps  = -3:0.01:8; %<-------------------
%     % define leftVideo samples to use
%     leftSamps = cell(1,numel(steps));
%     for k = 1:numel(steps)
%         ts = [];
%         for iev = 1:numel(leftOnset)
%             ts = [ts leftOnset(iev)+plvWin(1)+steps(k)*lfpFs:leftOnset(iev)+plvWin(2)+steps(k)*lfpFs];
%         end
%         leftSamps{k} = uint8(ts);
%     end
%     % define rightVideo samples to use
%     rightSamps = cell(1,numel(steps));
%     for k = 1:numel(steps)
%         ts = [];
%         for iev = 1:numel(rightOnset)
%             ts = [ts rightOnset(iev)+plvWin(1)+steps(k)*lfpFs:rightOnset(iev)+plvWin(2)+steps(k)*lfpFs];
%         end
%         rightSamps{k} = uint8(ts); % make sure they are integers so that we can use as index later
%     end
% 
%     pmat   = nan(numel(ppcChans)*numel(pulChans),numel(steps));
%     vmat   = nan(numel(ppcChans)*numel(pulChans),numel(steps));
%     count = 0;
%     for ichan = 1:numel(ppcChans)
%         display(num2str(ichan))
%         for jchan = 1:numel(pulChans)
%             count = count + 1;
%             for istep = 1:numel(steps)
%                 display(num2str(istep));
%                 phaseDiff = angPPC(ichan,leftSamps{istep}) - angPul(jchan,leftSamps{istep});
%                 pmat(count,istep) = abs(nanmean(exp(1i*phaseDiff)));
%                 phaseDiff = angPPC(ichan,rightSamps{istep}) - angPul(jchan,rightSamps{istep});
%                 vmat(count,istep) = abs(nanmean(exp(1i*phaseDiff)));                
%             end
%         end
%     end
%     avLvideoPLV(f,:) = nanmean(pmat);
%     avRvideoPLV(f,:) = nanmean(vmat);
    
    cd(recPath);
    savePath = ['..\..\Analyzed\' recName '\Spectra\'];
    if ~exist(savePath,'dir') 
    mkdir(savePath); end
    save([savePath 'LateralVideoSpectra.mat'],'foi','tvec','powLvideo_pul','powLvideo_ppc','powLvideo_vc','powRvideo_pul','powRvideo_ppc','powRvideo_vc','-v7.3');
end

doPlot = 1;
if doPlot == 1
    load('LateralVideoSpectra.mat');
    avgXSpec = squeeze(nanmean(powLvideo_pul ,2));
    avgYSpec = squeeze(nanmean(powLvideo_ppc ,2));
    avgZSpec = squeeze(nanmean(powLvideo_vc ,2));
    avgXSpecR = squeeze(nanmean(powRvideo_pul ,2));
    avgYSpecR = squeeze(nanmean(powRvideo_ppc ,2));
    avgZSpecR = squeeze(nanmean(powRvideo_vc ,2));
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    % Compute ticks for plotting
    fois = [1 2 4 8 16 32 64];
    %fois = [5 10 15 20 25 30];
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(foi-fois(fi)));
        tickLoc(fi) = bb(1);
    end
    tickLabel = string(fois);
    
    subplot(2,3,1)
        imagesc(tvec, foi, avgXSpec)
        xlabel('Time to video [s]'); ylabel('Frequency [Hz]');  title('Left side video');
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-0.7 2.4]);
        cl = colorbar('northoutside'); ylabel(cl,'LPl z-scored power','FontSize',12)
    subplot(2,3,2)
        imagesc(tvec, foi, avgYSpec)
        xlabel('Time to video [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-0.7 2.4]);
        cl = colorbar('northoutside'); ylabel(cl,'PPC z-scored power','FontSize',12)   
    subplot(2,3,3)
        imagesc(tvec, foi, avgZSpec)
        xlabel('Time to video [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-0.7 2.4]);
        cl = colorbar('northoutside'); ylabel(cl,'VC z-scored power','FontSize',12)  
    subplot(2,3,4)
        imagesc(tvec, foi, avgXSpecR)
        xlabel('Time to video [s]'); ylabel('Frequency [Hz]');  title('Right side video');
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-0.7 2.4]);
        cl = colorbar('northoutside'); ylabel(cl,'LPl z-scored power','FontSize',12)
    subplot(2,3,5)
        imagesc(tvec, foi, avgYSpecR)
        xlabel('Time to video [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-0.7 2.4]);
        cl = colorbar('northoutside'); ylabel(cl,'PPC z-scored power','FontSize',12)   
    subplot(2,3,6)
        imagesc(tvec, foi, avgZSpecR)
        xlabel('Time to video [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-0.7 2.4]);
        cl = colorbar('northoutside'); ylabel(cl,'VC z-scored power','FontSize',12)
        
    colormap(jet)
    savefig(fig, 'LateralVideoSpectrogram.fig','compact');
    saveas(fig, 'LateralVideoSpectrogram.png');
    %savefig(fig, [savePath 'avgFuncCon.fig'],'compact');
    %saveas(fig, [savePath 'avgFuncCon.png']);

end
