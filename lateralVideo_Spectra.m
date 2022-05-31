function lateralVideo_Spectra(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
    evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
    regionLFP, regionChn, rootAnalysisDir)

% modified from Iain's code is_imagesVideoSpectra.m
% AH 2018.6.29
spectra = struct;
doPlot = 1;
dispstat('','init'); % One time only initialization


win = round(twin*lfpFs);
tvec = (win(1):win(2))/1000; % time vector in sec
leftOnset = round(evtTimes{1}*lfpFs);
rightOnset = round(evtTimes{2}*lfpFs);

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if     linORlog == 1
        foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
        foi      = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end

saveAnalysisDir    = rootAnalysisDir;

    if ~exist(join(saveAnalysisDir),'dir'); mkdir(join(saveAnalysisDir));
    elseif length(dir([saveAnalysisDir '*.fig'])) >= length(condID) %each condition generate a fig file
        fprintf('Record %s already analyzed \n',saveAnalysisDir'); 
        %if skipRec == 1; continue; end  % to skip already analyzed records
    end
    fprintf('\nWorking on %s \n',saveAnalysisDir');  
    
% initialise matrices
powLvideo_pul = nan(length(foi),numel(regionChn{1}),diff(win)+1); 
powLvideo_ppc = nan(length(foi),numel(regionChn{2}),diff(win)+1); 
powLvideo_vc  = nan(length(foi),numel(regionChn{3}),diff(win)+1);
powRvideo_pul = nan(length(foi),numel(regionChn{1}),diff(win)+1);
powRvideo_ppc = nan(length(foi),numel(regionChn{2}),diff(win)+1);
powRvideo_vc  = nan(length(foi),numel(regionChn{3}),diff(win)+1);

% get wavelets
morWav = is_makeWavelet(foi,lfpFs);

% for iCond = condID
%     %if iCond == 1 && length(dir([saveAnalysisDir '*.fig']))>=1;continue;end
%     evtTime = evtTimes{iCond};
    
for f = 1:numFreqs
    dispstat(sprintf([num2str(f) '/' num2str(numFreqs)])); %dispstat(sprintf('Progress %d%%',i),'timestamp');
    % Convolve nFreq x nChannel x nTime
    xspec(f,:)  = conv2(regionLFP{1},morWav{f},'same'); % convolve LFP with wavelet, matrix of complex numbers
    yspec(f,:)  = conv2(regionLFP{2},morWav{f},'same');
    zspec(f,:)  = conv2(regionLFP{3},morWav{f},'same');
end

% z-score power over channels (If dim = 2, then zscore uses the means and standard deviations along the rows of X.)
xZpow  = zscore(abs(xspec).^2,[],2);
yZpow  = zscore(abs(yspec).^2,[],2);
zZpow  = zscore(abs(zspec).^2,[],2);

% get angle
xang = angle(xspec); clear xspec
yang = angle(xspec); clear xspec
zang = angle(xspec); clear xspec 

for iCond = condID
    if iCond == 1 && length(dir([saveAnalysisDir '*.fig']))>=1;continue;end
    event = evtTimes{iCond};     
    
    % Cut out spectral data around events
    tsamps = round(twin*fs);
    baseSamps = round(baseTwin*fs);
    fprintf('numel(event)=%d, numFreqs=%d, diff(tsamps)+1 =%d\n', ...
        numel(event),numFreqs,diff(tsamps)+1);
    xmat = nan(numel(event),numFreqs,size(xZpow,2),diff(tsamps)+1);
    ymat = nan(numel(event),numFreqs,size(yZpow,2),diff(tsamps)+1);
    zmat = nan(numel(event),numFreqs,size(zZpow,2),diff(tsamps)+1);
    xBasemat = nan(numel(event),numFreqs,size(xZpow,2),diff(baseSamps)+1);
    yBasemat = nan(numel(event),numFreqs,size(yZpow,2),diff(baseSamps)+1);
    zBasemat = nan(numel(event),numFreqs,size(zZpow,2),diff(baseSamps)+1);
    
    for iev = 1:numel(event)
        evSamp = round(event(iev)*fs);
        xmat(iev,:,:,:) = xZpow(:,:,evSamp+tsamps(1):evSamp+tsamps(2)
    spectra
end

    % LeftVideo Pul
    smat = nan(length(leftOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(leftOnset); smat(iev,:,:) = zpow(:,leftOnset(iev)+win(1):leftOnset(iev)+win(2)); end
    powLvideo_pul(f,:,:) = squeeze(nanmean(smat,1));
    
    % RightVideo Pul
    smat = nan(length(rightOnset),size(zpow,1),diff(win)+1);
    for iev = 1:numel(rightOnset); smat(iev,:,:) = zpow(:,rightOnset(iev)+win(1):rightOnset(iev)+win(2)); end
    powRvideo_pul(f,:,:) = squeeze(nanmean(smat,1));
    
    % PPC second
    C_ppc = conv2(regionLFP{2},morWav{f},'same'); % convole PPC LFP with wavelet
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
    C_vc = conv2(regionLFP{3},morWav{f},'same'); % convole PPC LFP with wavelet
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
    
    
    save([saveAnalysisDir 'LateralVideoSpectra.mat'],'foi','tvec','powLvideo_pul','powLvideo_ppc','powLvideo_vc','powRvideo_pul','powRvideo_ppc','powRvideo_vc','-v7.3');
end


if doPlot == 1
    load([saveAnalysisDir 'LateralVideoSpectra.mat']);
    avgXSpec = squeeze(nanmean(powLvideo_pul ,2));
    avgYSpec = squeeze(nanmean(powLvideo_ppc ,2));
    avgZSpec = squeeze(nanmean(powLvideo_vc ,2));
    avgXSpecR = squeeze(nanmean(powRvideo_pul ,2));
    avgYSpecR = squeeze(nanmean(powRvideo_ppc ,2));
    avgZSpecR = squeeze(nanmean(powRvideo_vc ,2));
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 screensize(4)-150]);
    % Compute ticks for plotting
    fois = [2 4 8 16 32 64 128];
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

end
    

