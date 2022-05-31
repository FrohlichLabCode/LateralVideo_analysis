% This master file calls AH_LateralVideoSpectra.m
clear
tic

cluster = 0;
skipRec = 0;
linORlog = 2; %freqs of interest: 1=linear 2=log
usePhotodiode = 1;
MedianorPCA = 0; 
animals = {'0168','0169'};


for iAnimal = 1:numel(animals)
    animalCode = animals{iAnimal};

    if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
        addpath(genpath( 'D:/Dropbox (Frohlich Lab)/Codebase/CodeAngel/Ephys'));
        PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
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
    
    fileInfo   = dir([PreprocessDir animalCode '_LateralVideo*']); % detect files to load/convert  '_LateralVideo*' can't process opto


% loop through each recording
for irec = 1:numel(recNames)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    if cluster==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180724', 'InputFormat', 'yyyyMMdd'); continue;end
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/'];
    rootBehavDatDir   = [BehavDatDir recName];

    AH_LateralVideoSpectra(rootPreprocessDir,recName,animalCode)
    
end
end