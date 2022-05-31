clear

addpath(genpath('C:/Users/angel/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
skipRec = 1;
usePhotodiode = 0;

animalCode = '0172';

twin = [-3 8];
baseTwin = [-2.5 -0.5];
twins = {twin};
baseTwins = {baseTwin};
    
PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
fileInfo   = dir([PreprocessDir animalCode '_MozartAudio*']); % detect files to load/convert  '_LateralVideo*'

linORlog = 2;
if linORlog == 1
    numFreqs = 200;
    lowFreq  = 1;
    highFreq = 80;
    foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
    fois = [2, 5:5:highFreq];
    tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
elseif linORlog == 2
    numFreqs = 80;
    lowFreq  = 2;
    highFreq = 128;
    foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
    fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
    tickLabel = string(fois);
end

for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end



% loop through each recording
for irec = 7:numel(fileInfo)
    recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
    splitName = strsplit(recName,'_');
    if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180823', 'InputFormat', 'yyyyMMdd'); continue;end

    rootPreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/' recName '/'];
    rootAnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/' recName '/LFP/'];
    rootBehavDatDir   = ['D:/FerretData/' animalCode '/behav/' recName];
    if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
        fprintf('Record %s already analyzed \n',recName'); 
        if skipRec == 1; continue; end; end

% region info
switch animalCode
    case '0172'
    regionNames = {'FC','LPl','PPC','VC'};
    numRegion   = numel(regionNames);
    noPhotoData = 0;
    %noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd'))...
    %              +(datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') == datetime('20180718', 'InputFormat', 'yyyyMMdd'));
    allChn = {[1:16],[17:32],[33:48],[49:64]};
end

% load and process ttl data in ephys
display(['Processing TTL and LFP for rec ' recName]);
load([rootPreprocessDir 'triggerData']);
load([rootPreprocessDir 'lfp/lfpMat']); % including lfpFs=1000

validNumChn = numel(allChn);
fileID  = fopen([rootBehavDatDir '.log']);
formatSpec = '%f %s %s %f %f %f %f %f %f %f %s %f';
LogFile = textscan(fileID,formatSpec,'HeaderLines',5,'Delimiter', '\t');
fclose(fileID);

% detect location
Col_TrialNumber = 1;
Col_TrialType = 2; % video
Col_Code = 3;      % trial side and video name

if str2num(splitName{3}) <= 2 % first 2 sessions only Mozart condition
    cellInd = [1:length(LogFile{Col_Code})];
else
    for condID = 1:4
    cellCond = cellfun(@regexp,LogFile(Col_Code),{[num2str(condID) +';']},'UniformOutput',false); %match {1;...} output {[],1,[],[],1...}
    cellInd{condID} = cellfun(@isempty,cellCond{:}) == 0; %output [0,1,0,0,1...] 
    end
end

% detect stim onsets
ttlInd = 1;
rawFs = 30000;
ttlOnset  = find(diff(triggerData(ttlInd,:))==1)./rawFs; % in sec
ttlOffset = find(diff(triggerData(ttlInd,:))==-1)./rawFs;
onset = ttlOnset;

% if usePhotodiode == 1 && noPhotoData == 0 % not a trial without photo data
%     load([rootPreprocessDir 'adc_data']);
%     photoInd = 4; % one eye-tracking
%     photoThreshold = 3.2; % used 3.3V ground
%     if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') >= datetime('20180807', 'InputFormat', 'yyyyMMdd')
%         photoInd = 7; % bilateral eye-tracking
%         photoThreshold = 0.12; % used 0 ground
%     end
%     photoRaw = adc_data(photoInd,:);
% %     temp = zeros(1,length(photoRaw)); % this was used with find(diff(temp(i:i+rawFs*0.1))==1,1,'last');
% %     temp(photoRaw>photoThreshold) = 1; % 1=gray(no stim), 0 = black(stimOn)
%     photoOnset = [];
%     for i = ttlOnset*rawFs % for each trigger time point
%         rawAvg = movmean(photoRaw(i:i+rawFs*0.1), 200);% about every <200 sample, photodiode has a dip, this will smooth it out
%         offset = find(rawAvg<photoThreshold+0.04,1,'first'); %moving average, use 0.16 as threshold %0169 use:find(diff(temp(i:i+rawFs*0.1))==1,1,'last');
%         if isempty(offset); offset = 0;end
%         iphotoOnset = (i+offset)/rawFs; 
%         photoOnset = [photoOnset, iphotoOnset]; % the last drop->0 is when gray->black in photodiode
%     end
%     onset = photoOnset;
% else
%     onset = ttlOnset + 0.05; % on average, trigger data is 0.0529 sec earlier than visual stimuli onset
% end

%Visualize the diff in timing
%figure();plot(onset,1,'r.');hold on;plot(photoOnset,1,'b.'); 

condNames = {'Mozart','Bach','pureC','whitenoise'};
if str2num(splitName{3}) <= 2
    condID = [1];
elseif str2num(splitName{3}) <= 4
    condID = [1,3,4];
else
    condID = [1,2,3,4];
end
numConds = numel(condID);

%% extract snippits of lfp
condCount = 1;
evtTime = onset;
evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window
evtTimes{1} = evtTime;
for iCond = 1:numConds
    iCondID = condID(iCond);
    evtTime = onset(cellInd{iCondID}); % i.e. Mozart Onset
    evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window
    evtTimes{iCondID} = evtTime;
end

save([rootPreprocessDir 'eventTimes'], 'condNames','condID', 'cellInd', 'evtTimes','twins','baseTwins', 'lfpFs');

xRange = [evtTimes{1}(2)+twin(1),evtTimes{1}(3)+twin(2)]; % in sec, window around 2nd left video
tRange = round(xRange(1)*lfpFs):round(xRange(end)*lfpFs); %column index

for iRegion = 1:numRegion
    nChn = numel(allChn{iRegion});
    nCol = 4;
    nRow = ceil(nChn/nCol);
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
    for iChn = 1:nChn
        subplot(nRow,nCol,iChn)
        plot(tRange/lfpFs, lfpMat(allChn{iRegion}(iChn),tRange));
        vline(evtTimes{1},'k-'); %left video
        vline(evtTimes{1}+5,'k--'); %left video          
%             vline(evtTimes{2},'r-'); %right video
%             vline(evtTimes{2}+5,'r--'); %left video
        xlim(xRange);
        title(['chn ' num2str(iChn)]);
        ylim([-200,200]);
        if iChn>(nRow-1)*nCol % last row
            xlabel('Time [s]');
        end

        if mod(iChn,nCol)==0 
            set(gca,'yticklabel',{[]})
            ylabel('')
        else
            ylabel('uV') %left most column has unit
        end
    end
    if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
    savefig(fig, [rootAnalysisDir regionNames{iRegion} '_LFP.fig'],'compact');
    saveas(fig, [rootAnalysisDir regionNames{iRegion} '_LFP.png']);

    % Compute spectrogram of each channel
    window   = 0.5*1024;
    noverlap = round(3*window/4);
    specEpoch     = nan(size(lfpMat(allChn{iRegion}(iChn),tRange),1),numFreqs);

    fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
    for iChn = 1:nChn
        subplot(nRow,nCol,iChn)
        [s,f,t] = spectrogram(lfpMat(allChn{iRegion}(iChn),tRange),window,noverlap,foi,lfpFs);
        tvec = t+tRange(1)/lfpFs;
        imagesc(tvec,1:numel(foi),pow2db(abs(s).^2));
        vline(evtTimes{1},'k-'); %left video
        vline(evtTimes{1}+5,'k--'); %left video
%             vline(evtTimes{2},'r-'); %right video
%             vline(evtTimes{2}+5,'r--'); %left video
        title(['chn ' num2str(iChn)]);
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        %ylim([1,30]); 
        if iChn>(nRow-1)*nCol % last row
            xlabel('Time [s]');
        end

        if mod(iChn,nCol)==1
            ylabel('uV') %left most column has unit
        else
            set(gca,'yticklabel',{[]})
            ylabel('')                
        end

        if mod(iChn,nCol)==0
            cl = colorbar('eastoutside'); ylabel(cl,'Power [db]','FontSize',12)        
        end
    end
    colormap(jet)
    if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
    savefig(fig, [rootAnalysisDir regionNames{iRegion} '_spectrogram.fig'],'compact');
    saveas(fig, [rootAnalysisDir regionNames{iRegion} '_spectrogram.png']);

    % plot spectra for ~half of the recording
    %tRange2 = round(xRange(1)*lfpFs):round(size(lfpMat,2)*1/2);  % get the 2nd quarter of the recording
    [pxx,f] = pwelch(lfpMat(allChn{iRegion},tRange)',window,noverlap,foi,lfpFs); % power density spectra; PSD is computed independently for each column and stored in the corresponding column of pxx
    fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
    for iChn = 1:nChn
        subplot(nRow,nCol,iChn)
        loglog(f, pxx, 'color', 0.6*[1 1 1]); hold on
        loglog(f, nanmedian(pxx, 2), 'color', 'r')
        loglog(f, pxx(:, iChn), 'linewidth', 2,'color', 'b');
        xlim([f(1) f(end)]); %ylim([1e2 1e6]); 
        title(['chn ' num2str(iChn)],'color', 'b');
        if iChn>(nRow-1)*nCol % last row
            xlabel('Frequency [Hz]');
        end
        if mod(iChn,nCol)==1
            ylabel('Power') %left most column has unit
        else
            set(gca,'yticklabel',{[]})
            ylabel('')                
        end            
    end       

    if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
    savefig(fig, [rootAnalysisDir regionNames{iRegion} '_spectra.fig'],'compact');
    saveas(fig, [rootAnalysisDir regionNames{iRegion} '_spectra.png']);     

    close all % close all figure
end    
    
    
    %condCount = condCount + 1;    

end