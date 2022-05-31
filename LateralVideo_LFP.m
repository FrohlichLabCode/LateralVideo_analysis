% This script will prepross events and save eventTimes in preprocess
% directory, if choose to plot, it will also plot LFP for inspection
% purpose

% AH

clear

addpath('C:/Users/angel/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/');
skipRec = 0;
usePhotodiode = 1;
doPlot = 0;

animalCode = '0172';

PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
fileInfo   = dir([PreprocessDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*'

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
for irec = 1:numel(fileInfo)
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
    noPhotoData = strcmp(splitName{4}, '20190111') + strcmp(splitName{4}, '20190112');%photodiode to INTAN connection loose
    allChn = {[1:16],[17:32],[33:48],[49:64]};
    
    case '0168'
    regionNames = {'LPl','PPC','VC'};
    numRegion   = numel(regionNames);
    noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd'))...
                  +(datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') == datetime('20180718', 'InputFormat', 'yyyyMMdd'));
    allChn = {[1:16],[17:48],[49:64]};
    case '0169' % 0169 region info
    regionNames = {'lLPl','lPPC','rPPC','rLPl','lrVCEEG'};
    numRegion   = numel(regionNames);
    noPhotoData = datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd')...
                  + strcmp(splitName{3},'039')+ strcmp(splitName{3},'042'); % full size video
    allChn = {[1:16],[17:32],[33:48],[49:64],[65:66]};
end

twin = [-3 8];
baseTwin = [-3,-0.5];

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

cellLeft = cellfun(@regexp,LogFile(Col_Code),{'1;'},'UniformOutput',false);
leftInd = cellfun(@isempty,cellLeft{:}) == 0;
cellRight = cellfun(@regexp,LogFile(Col_Code),{'2;'},'UniformOutput',false);
rightInd = cellfun(@isempty,cellRight{:}) == 0;

% detect stim onsets
ttlInd = 1;
rawFs = 30000;

if usePhotodiode == 1 && noPhotoData == 0 % not a trial without photo data
    load([rootPreprocessDir 'adc_data']);
    photoInd = 4; % one eye-tracking
    photoBaseline = 3.2; % used 3.3V ground
    if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') >= datetime('20180807', 'InputFormat', 'yyyyMMdd')
        photoInd = 7; % bilateral eye-tracking
        photoBaseline = 0.12; % used 0 ground
    end
    photoRaw = adc_data(photoInd,:);
%     temp = zeros(1,length(photoRaw)); % this was used with find(diff(temp(i:i+rawFs*0.1))==1,1,'last');
%     temp(photoRaw>photoThreshold) = 1; % 1=gray(no stim), 0 = black(stimOn)
    onset = getStimOnset(triggerData(ttlInd,:), photoRaw, photoBaseline);
end

%Visualize the diff in timing
%figure();plot(onset,1,'r.');hold on;plot(photoOnset,1,'b.'); 

condNames = {'lVideo','rVideo'};
condID = [1 2];
numConds = numel(condID);

%% extract snippits of lfp
condCount = 1;
for iCond = 1:numConds
    
    if iCond == 1
        evtTime = onset(leftInd); % i.e. leftOnset
    elseif iCond == 2
        evtTime = onset(rightInd); % i.e. rightOnset
    end
    evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window    
    evtTimes{iCond} = evtTime;
end

save([rootPreprocessDir 'eventTimes'], 'condNames', 'evtTimes','twin','baseTwin','condID', 'lfpFs','usePhotodiode'); 

if doPlot == 1
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
            vline(evtTimes{2},'r-'); %right video
            vline(evtTimes{2}+5,'r--'); %left video
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
            vline(evtTimes{2},'r-'); %right video
            vline(evtTimes{2}+5,'r--'); %left video
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
        savefig(fig, [rootAnalysisDir regionNames{iRegion} '_spec.fig'],'compact');
        saveas(fig, [rootAnalysisDir regionNames{iRegion} '_spec.png']);
        close all % close all figure
    end
end

end