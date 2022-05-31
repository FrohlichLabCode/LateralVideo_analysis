function MozartAudio_rec_cluster(animalCode,irec, fileInfo,analysisType,folderSuffix, PreprocessDir, AnalysisDir, BehavDatDir, GroupAnalysisDir,...
        cluster, skipRec, linORlog, MedianorPCA,usePhotodiode)

    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionName = splitName{3}; 
    %if cluster==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180822', 'InputFormat', 'yyyyMMdd'); continue;end
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/']; %<---------------
    rootBehavDatDir   = [BehavDatDir recName];
        
    if ~exist(join(rootAnalysisDir),'dir') 
        mkdir(join(rootAnalysisDir)); fprintf('\nWorking on record %s =============== \n',recName'); end


    % load and process ttl data in ephys
    lfpMat = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat'); % don't feed in denoised data with NaN values
    load([rootPreprocessDir 'triggerData']);
    % don't need to load behavior here
        

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
        noPhotoData = 1;
        fullVideo = 1; % there is no fullVideo session
%     case '0168'
%         regionNames = {'LPl','PPC','VC'};
%         numRegion   = numel(regionNames);
%         regionPairs = {[1,2],[1,3],[2,3]};
%         noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180806', 'InputFormat', 'yyyyMMdd'))...
%                       +strcmp(splitName{4}, '20180718');
%         fullVideo   = 0; % there is no fullVideo session
%     
%     case '0169' % 0169 region info
%     regionNames = {'lLPl','lPPC','rPPC','rLPl','lVCEEG','rVCEEG'};
%     numRegion   = numel(regionNames);
%     regionPairs = {[2,5],[1,5],[1,2],[4,6],[3,6],[4,3],[1,4],[2,3],[5,6]};
%     fullVideoSessions = {'039','042','045','048','051','054','058','062','064','069','073','077','083','086','092','100','136'};
%     if find(contains(fullVideoSessions,sessionName))>0; fullVideo = 1; else fullVideo = 0; end
%     noPhotoSessions = {'039','042','045','048','051','054','058','061','062','064','067','068','069','073','074',...
%                        '077','078','083','084','086','087','091','092','097','100','103','109','113','117','123','130','131','132','133','134','135','136'}; % full size video blocked photodiode
%     noPhotoData = (datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180628', 'InputFormat', 'yyyyMMdd'));
%     if find(contains(noPhotoSessions,sessionName))>0; noPhotoData = noPhotoData + 1;end
end

%% select LFP to load
[lfp.validChn, lfp.reorderedChn] = keepChn(recName);
if MedianorPCA == 0
    %lfpMat = lfpMat; % use the raw lfp signal since GC can't deal with NaN
    lfpFs = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpFs');  
    for i = 1:numel(lfp.validChn) 
        regionChn{i} = lfp.validChn{i}; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(lfp.validChn{i},:); % use original channel number!
    end
    
elseif MedianorPCA == 1
    load([rootPreprocessDir 'lfp/lfpValid']);
    lfpMat = lfp.median;
    lfpFs  = lfp.Fs;
    for i = 1:size(lfp.median,1) %lfp.median is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end

elseif MedianorPCA == 2
    load([rootPreprocessDir 'lfp/lfpValid']);
    lfpMat = lfp.PCA;
    lfpFs  = lfp.Fs;
    for i = 1:size(lfp.PCA,1) %lfp.PCA is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end
    
elseif MedianorPCA == 3
    lfpFs = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpFs');
    for i = 1:numel(lfp.validChn)
        regionLFP{i} = lfpMat(lfp.validChn{i}(1),:); %pick the first valid channel from each region
        regionChn{i} = i; % Pulvinar, PPC, VC
    end
end

load([rootPreprocessDir 'eventTimes.mat']);

% sanity check for event times
for i = 1:numel(evtTimes)
    evtTime = evtTimes{i};
    display(['evtTimes for ' condNames{i} ':' num2str(round(evtTime))]);
end

% twin = twins{1};
% baseTwin = baseTwins{1};
twin = [-3,8];
baseTwin = [-2.5,-0.5];
%% for each region pairs
% lateralVideo_Spectra(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%     evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
%     regionLFP, regionChn, rootAnalysisDir);
% 
% regionPair_FunConn(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%    evtTimes,twin,baseTwin, condNames, condID, regionPairs, regionNames, ...
%    regionLFP, regionChn, rootAnalysisDir);


% if length(dir([rootAnalysisDir '*.fig'])) >= 1 %each condition generate a fig file
%     fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
%     if skipRec == 1; return; end  % to skip already analyzed records
% end
% fprintf('\nWorking on %s \n',rootAnalysisDir');
% lateralVideo_Pupil(skipRec, lfpFs,...
%     evtTimes, condNames, condID,...
%     rootPreprocessDir, rootAnalysisDir);

% if length(dir([rootAnalysisDir '*.fig'])) >= 3 %each condition generate a fig file
%     fprintf('Record %s already analyzed \n',rootAnalysisDir'); 
%     if skipRec == 1; return; end  % to skip already analyzed records
% end
% display(['computing PSTH ' recName]);
% lateralVideo_PSTH(evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
%     regionChn, rootPreprocessDir, rootAnalysisDir, animalCode);
lateralVideo_saccades(evtTimes,twin,baseTwin, condNames, condID, regionNames, ...
    regionChn, rootPreprocessDir, rootAnalysisDir, animalCode);


% region_spec_by_trial(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%     evtTimes,twins,baseTwins, condNames, condID, regionNames, ...
%     regionLFP, regionChn, sessionName, GroupAnalysisDir);

