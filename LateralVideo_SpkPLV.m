clear

cluster = 0;
skipRec = 0;
linORlog = 2; %freqs of interest: 1=linear 2=log
usePhotodiode = 1;
MedianorPCA = 1; 
animals = {'0168','0169'};

wavHil = 1; % 0 for hilbert mean lfp, 1 for wavelet each chn, 2 for hilbert each chn
transformLab = {'Wave', 'Hil'};
newFs = 400; % to downsample the lfp for faster computing
lfpLab = {'Evoked', 'Indu'};

%% Define frequencies of interest.
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
if wavHil == 1
    wavs = is_makeWavelet(foi,newFs); % make wavelets from FOI frequencies
end
dat.foi = foi;

for iAnimal = 1:numel(animals)
    animalCode = animals{iAnimal};

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    fileInfo      = dir([PreprocessDir '0*']);
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];

    %code for initialising parallel computing
    numCore = 24; % USR DEFINE
    parpool('local',numCore);   
end

% Different folder name if use median or PCA LFP
if MedianorPCA == 0;     folderSuffix = '_validChns'; %use all valid channels
elseif MedianorPCA == 1; folderSuffix = '_median'; %use median of all valid channels
elseif MedianorPCA == 2; folderSuffix = '_PCA';
end

fileInfo   = dir([PreprocessDir animalCode '_LateralVideo*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto


%% extract snippits of lfp

for irec = 1%:numel(fileInfo)     
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    %if cluster==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180724', 'InputFormat', 'yyyyMMdd'); continue;end
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/SpkPLV' folderSuffix '/'];
    rootBehavDatDir   = [BehavDatDir recName];
        
if ~exist(join(rootAnalysisDir),'dir') 
    mkdir(join(rootAnalysisDir)); fprintf('\nWorking on record %s =============== \n',recName'); end


%% load and process ttl data in ephys
load([rootPreprocessDir 'lfp/lfpValid']);
load([rootPreprocessDir 'triggerData']);
% don't need to load behavior here 

%% load spike data
numSpkFile = numel(dir([rootPreprocessDir 'spikes/' 'spk*'])); % should match number of channels in recording
spkCell = cell(1,numSpkFile);
for iChn = 1:numSpkFile
    load([rootPreprocessDir 'spikes/' 'spk_' num2str(iChn)])
    spkCell{iChn} = spkTime;
end

% region info
switch animalCode
    case '0168'
    regionNames = {'LPl','PPC','VC'};
    numRegion   = numel(regionNames);
    regionPairs = {[1,2],[1,3],[2,3]};
    noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd'))...
                  +(datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') == datetime('20180718', 'InputFormat', 'yyyyMMdd'));

    case '0169' % 0169 region info
    regionNames = {'lLPl','lPPC','rPPC','rLPl','lVCEEG','rVCEEG'};
    numRegion   = numel(regionNames);
    regionPairs = {[1,2],[4,3],[1,4],[2,3],[5,6],[1,5],[4,6],[2,5],[3,6]};
    noPhotoData = datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd');
end

%% select LFP type to load
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
onset  = find(diff(triggerData(ttlInd,:))==1)./rawFs; % in sec
offset = find(diff(triggerData(ttlInd,:))==-1)./rawFs;

% use photodiode data if available since it's more accurate than ttl trigger
if usePhotodiode == 1 && noPhotoData == 0 % not a trial without photo data
    load([rootPreprocessDir 'adc_data']);
    photoInd = 4;
    photoRaw = adc_data(photoInd,:);
    temp = zeros(1,length(photoRaw));
    temp(photoRaw>3.2) = 1;
    photoOnset = [];
    for i = onset*rawFs % for each trigger time point    
        photoOnset = [photoOnset, (i+find(diff(temp(i:i+rawFs*0.1))==1,1,'last'))./rawFs]; % the last drop->0 is when gray->black in photodiode
    end
    onset = photoOnset;
else
    onset = onset + 0.05; % on average, trigger data is 0.0529 sec earlier than visual stimuli onset (1 sample session)
end

%% get conds from behav dat
condNames = {'Left','Right'};
condID = [1 2];
numConds = numel(condID);

%% get event times for each condition
for iCond = condID
    if iCond == 1
        evtTime = round(onset(leftInd)); % i.e. leftOnset
    elseif iCond == 2
        evtTime = round(onset(rightInd)); % i.e. rightOnset
    end
    evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window
    evtTimes{iCond} = evtTime;
end
    
    
    numTotChn = size(lfpMat,1);
    numTotSamp = size(lfpMat,2);
    
    lfpDownsampled = resample(lfpMat',newFs,lfpFs)';
    lfpInput = lfpDownsampled;
    numLFPSamp = size(lfpInput,2);    
    
    avgPowCond = nan(numConds,numFreqs,numTotChn);
    stdPowCond = nan(numConds,numFreqs,numTotChn);
    for iCond = 1:numConds
        
        display(['Loading data for ' recName ' cond: ' num2str(iCond)]);        
        evtTime = evtTimes{iCond}; % trialOnset in seconds
        numEvt = numel(evtTime);
        
        for f = 1:numFreqs
            
            % define frequency parameters and wavelets; CHANGE
            Freq = foi(f);
            if wavHil == 0
                C = [];
            elseif wavHil == 1
                C = conv2(lfpInput,wavs{f},'same'); % dims are channels by time
            elseif wavHil == 2
                % Begin: filter and perform hilbert transform
                C = cz_computeHilbert(lfpInput, Freq, newFs);
            end
            %% compute spike phase locking
            
            C_mean = nan(numRegion,numLFPSamp);
            for iRegion = 1:numRegion
                meanLFP = nanmedian(lfpInput(regionChn{iRegion},:),1); % take average lfp
                C_mean(iRegion,:) = cz_computeHilbert(meanLFP, Freq, newFs);
            
            [evtSpkPLV,evtSpkAngle,numBins] = is_spkPLV_MeanLFP_V3_Carelli(spkCell,evtTime,analyzeTheseTrials,C_mean,newLFPFs,f,regionChn,twin);

            %sponSpkPLVAll(iCond,f,:)      = sponSpkPLV;
            evtSpkPLVAll(iCond,iRegion,f,:,:) = evtSpkPLV; % dims are condition, frequency, channel, bins
            evtSpkAngleAll(iCond,iRegion,f,:,:) = evtSpkAngle;

            [sponSpkPLV,evtSpkPLV,evtSpkAngle,evtPhase,numBins] = is_spkPLV_MeanLFP_V3(spkCell,evtTime,C,C_mean,newFs,f,regionChn,twin);
            %sponSpkPLVAll(iCond,f,:)      = sponSpkPLV;
            evtSpkPLVAll(iCond,iRegion,f,:,:) = evtSpkPLV; % dims are condition,region, frequency, channel, bins 
            evtSpkAngleAll(iCond,iRegion,f,:,:) = evtSpkAngle;
            evtPhaseAll{iCond,iRegion,f} = evtPhase;
            
            end
            
            
            
        end
        
    end
 
    dat.regionChn = regionChn;
    dat.twin = twin;
    dat.numBins = numBins;
    %dat.sponSpkPLVAll = sponSpkPLVAll;
    dat.evtSpkPLVAll = evtSpkPLVAll;
    dat.evtSpkAngleAll = evtSpkAngleAll;
    dat.evtPhaseAll = evtPhaseAll;
    dat.foi = foi;
    dat.condNames = session_output_data.condNames;
    
    save([rootPreprocessDir recName '_SpkPLV'],'dat','-v7.3');
    
end
end



