function lateralVideo_rec_cluster(animalCode,irec, fileInfo,analysisType,folderSuffix, PreprocessDir, AnalysisDir, BehavDatDir, GroupAnalysisDir,...
        cluster, skipRec, linORlog, MedianorPCA,usePhotodiode)

    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionID = splitName{3}; 

    %if cluster==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180822', 'InputFormat', 'yyyyMMdd'); continue;end
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/']; %<---------------
    rootBehavDatDir   = [BehavDatDir recName];
        
    if ~exist(join(rootAnalysisDir),'dir') 
        mkdir(join(rootAnalysisDir)); fprintf('\nWorking on record %s =============== \n',recName'); end

%%
% Define frequencies of interest. Linear spacing for Phase slope index, and
% logarithmic spacing for all other methods.

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if linORlog == 1
    numFreqs = 100;
    lowFreq  = 2;
    highFreq = 30;
    %foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
    numFreqs = 150;
    lowFreq  = 2;
    highFreq = 128;
    %foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end

% region info
switch animalCode
    case '0172'
        regionNames = {'FC','LPl','PPC','VC'};
        numRegion   = numel(regionNames);
        regionPairs = {[1,3],[2,3],[2,4],[3,4]};
        noPhotoData = 0;
        noPhotoData = strcmp(splitName{4}, '20190111') + strcmp(splitName{4}, '20190112');%photodiode to INTAN connection loose
        allChn = {[1:16],[17:32],[33:48],[49:64]};       
        fullVideo   = 0; % there is no fullVideo session
    case '0168'
        regionNames = {'LPl','PPC','VC'};
        numRegion   = numel(regionNames);
        regionPairs = {[1,2],[1,3],[2,3]};
        noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180806', 'InputFormat', 'yyyyMMdd'))...
                      +strcmp(splitName{4}, '20180718');
        fullVideo   = 0; % there is no fullVideo session
    
    case '0169' % 0169 region info
    regionNames = {'lLPl','lPPC','rPPC','rLPl','lVCEEG','rVCEEG'};
    numRegion   = numel(regionNames);
    regionPairs = {[2,5],[1,5],[1,2],[4,6],[3,6],[4,3],[1,4],[2,3],[5,6]};
    fullVideoSessions = {'039','042','045','048','051','054','058','062','064','069','073','077','083','086','092','100','136'};
    if find(contains(fullVideoSessions,sessionID))>0; fullVideo = 1; else fullVideo = 0; end
    noPhotoSessions = {'039','042','045','048','051','054','058','061','062','064','067','068','069','073','074',...
                       '077','078','083','084','086','087','091','092','097','100','103','109','113','117','123','130','131','132','133','134','135','136'}; % full size video blocked photodiode
    noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd'));
    if find(contains(noPhotoSessions,sessionID))>0; noPhotoData = noPhotoData + 1;end
end

%% select LFP to load
% load and process ttl data in ephys
[lfpMat,lfpFs] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs'); % don't feed in denoised data with NaN values
[lfp.validChn, lfp.reorderedChn] = keepChn(recName);
if MedianorPCA == 0
    %lfpMat = lfpMat; % use the raw lfp signal since GC can't deal with NaN 
    for i = 1:numel(lfp.validChn) 
        regionChn{i} = lfp.validChn{i}; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(lfp.validChn{i},:); % use original channel number!
    end
    
elseif MedianorPCA == 1
    load([rootPreprocessDir 'lfp/lfpValid']);
    lfpMat = lfp.median;
    for i = 1:size(lfp.median,1) %lfp.median is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end

elseif MedianorPCA == 2
    load([rootPreprocessDir 'lfp/lfpValid']);
    lfpMat = lfp.PCA;
    for i = 1:size(lfp.PCA,1) %lfp.PCA is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end
    
elseif MedianorPCA == 3
    for i = 1:numel(lfp.validChn)
        regionLFP{i} = lfpMat(lfp.validChn{i}(1),:); %pick the first valid channel from each region
        regionChn{i} = i; % Pulvinar, PPC, VC
    end
end


%% parameteres
if exist([rootPreprocessDir 'eventTimes.mat']) 
    load([rootPreprocessDir 'eventTimes']);
    %condNames = {'lVideo','rVideo'};
    condID = [1 2];
    twin = [-3 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
    baseTwin = [-2.5 -0.5]; %<<<--- baseline to subtract from spectrogram 

else
    twin = [-3 8]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
    baseTwin = [-2.5 -0.5]; %<<<--- baseline to subtract from spectrogram 

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
    ttlOnset = find(diff(triggerData(ttlInd,:))==1)./rawFs; % in sec
    ttlOnset = find(diff(triggerData(ttlInd,:))==-1)./rawFs;

    if usePhotodiode == 1 && noPhotoData == 0 % not a trial without photo data
        load([rootPreprocessDir 'adc_data']);
        photoInd = 4; % one eye-tracking
        photoThreshold = 3.2; % used 3.3V ground
        if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') >= datetime('20180807', 'InputFormat', 'yyyyMMdd')
            photoInd = 7; % bilateral eye-tracking
            photoThreshold = 0.12; % used 0 ground
        end
        photoRaw = adc_data(photoInd,:);
    %     temp = zeros(1,length(photoRaw)); % this was used with find(diff(temp(i:i+rawFs*0.1))==1,1,'last');
    %     temp(photoRaw>photoThreshold) = 1; % 1=gray(no stim), 0 = black(stimOn)
        photoOnset = [];
        for i = ttlOnset*rawFs % for each trigger time point
            rawAvg = movmean(photoRaw(i:i+rawFs*0.1), 200);% about every <200 sample, photodiode has a dip, this will smooth it out
            offset = find(rawAvg<photoThreshold+0.04,1,'first'); %moving average, use 0.16 as threshold %0169 use:find(diff(temp(i:i+rawFs*0.1))==1,1,'last');
            if isempty(offset); offset = 0;end
            iphotoOnset = (i+offset)/rawFs; 
            photoOnset = [photoOnset, iphotoOnset]; % the last drop->0 is when gray->black in photodiode
        end
        onset = photoOnset;
    else
        onset = ttlOnset + 0.05; % on average, trigger data is 0.0529 sec earlier than visual stimuli onset
    end

    %Visualize the diff in timing
    %figure();plot(onset,1,'r.');hold on;plot(photoOnset,1,'b.'); 

    condNames = {'lVideo','rVideo'};
    condID = [1 2];
    numConds = numel(condID);

    %% get event times for each condition
    for iCond = condID

        if iCond == 1
            evtTime = onset(leftInd); % i.e. leftOnset
        elseif iCond == 2
            evtTime = onset(rightInd); % i.e. rightOnset
        end
        evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window

        evtTimes{iCond} = evtTime;
        twins{iCond}    = twin;
        baseTwins{iCond}= baseTwin;
    end
    % full screen sessions: merge left and right stim event time
    if fullVideo
        clear evtTimes;
        evtTimes{1} = sort([onset(leftInd) onset(rightInd)]);
        condID = 1;
        condNames = {'Full'};
    end

    % sanity check for event times
    for i = 1:numel(evtTimes)
        display(['evtTimes for ' condNames{i} ':' num2str(round(evtTimes{i}))]); 
    end

    % only need to save this for the first time running
    save([rootPreprocessDir 'eventTimes'], 'evtTimes', 'condNames', 'condID', 'lfpFs', 'usePhotodiode');
end



%% for all region pairs
% lateralVideo_Spectra(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%     evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
%     regionLFP, regionChn, rootAnalysisDir);

% regionPair_FunConn(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs, sessionID,...
%    evtTimes,twin,baseTwin, condNames, condID, regionPairs, regionNames, ...
%    regionLFP, regionChn, rootAnalysisDir, GroupAnalysisDir);


% if length(dir([rootAnalysisDir '*.fig'])) >= 2 %each condition generate a fig file
%     fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
%     if skipRec == 1; return; end  % to skip already analyzed records
% end
% fprintf('\nWorking on %s \n',rootAnalysisDir');
% lateralVideo_Pupil(skipRec, lfpFs,...
%     evtTimes,condNames, condID,...
%     rootPreprocessDir, rootAnalysisDir);

if length(dir([rootAnalysisDir '*.fig'])) >= 3 %each condition generate a fig file
    fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
    if skipRec == 1; return; end  % to skip already analyzed records
end
display(['computing PSTH ' recName]);
lateralVideo_PSTH(evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
    regionChn, rootPreprocessDir, rootAnalysisDir, animalCode);

% region_spec_by_trial(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%     evtTimes,twins,baseTwins, condNames, condID, regionNames, ...
%     regionLFP, regionChn, sessionName, GroupAnalysisDir);

