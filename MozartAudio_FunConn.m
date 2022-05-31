clear
tic

cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
usePhotodiode = 1;
MedianorPCA = 3; 
animals = {'0168','0169','0172'};
experiment = 'LateralVideo'; %'LateralVideo';
analysisType = 'FC';  %'Spectra','FC','Pupil','PSTH','Saccades'


for iAnimal = 3%1:numel(animals)
    animalCode = animals{iAnimal};

    if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
        addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
        PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
        GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/FC/'];
    elseif cluster == 1
        addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
        PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
        AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];
        GroupAnalysisDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/GroupAnalysis/FC/'];
       
        %code for initialising parallel computing
        numCore = 24; % USR DEFINE
        myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
    
    % Different folder name if use median or PCA LFP
    if strcmp(analysisType, 'FC') % if it is functional connectivity analysis
        if MedianorPCA == 0;     folderSuffix = '_validChns'; %use all valid channels
        elseif MedianorPCA == 1; folderSuffix = '_median'; %use median of all valid channels
        elseif MedianorPCA == 2; folderSuffix = '_PCA';
        elseif MedianorPCA == 3; folderSuffix = '_firstChn'; % randomly pick 1 channel            
        end
    else folderSuffix = '';
    end
    
    fileInfo   = dir([PreprocessDir animalCode '_' experiment '*']); % detect files to load/convert  '_LateralVideo*' can't process opto

% loop through each recording
for irec = 1:numel(fileInfo)%0169 -- 54:numel(fileInfo)
    % so that if one record doesn't work, others can still run
    MozartAudio_rec_cluster(animalCode,irec,fileInfo,analysisType,folderSuffix, PreprocessDir, AnalysisDir, BehavDatDir, GroupAnalysisDir,...
        cluster, skipRec, linORlog, MedianorPCA, usePhotodiode);

end % end of all records for an animal
end % end of all animals

if cluster == 1; delete(myPool);end

x = toc;
fprintf('time required =%f sec\n', x);


