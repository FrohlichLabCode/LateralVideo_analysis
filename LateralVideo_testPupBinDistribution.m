clear
tic

cluster = 0;
animalCode = '0169';
excludeWin = [-0.1 0.6]; % in seconds

if cluster == 0
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    fileInfo     = dir([PreprocessDir animalCode '_LateralVideo_*']);
elseif cluster == 1
    addpath(genpath('/nas/longleaf/home/zhouz/Code/is_pupilState/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir       = [ '/pine/scr/z/h/zhouz/OptoStim/toAnalyze/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    saveRootPath = [ '/pine/scr/z/h/zhouz/OptoStim/analyzed/'];
    mkdir(saveRootPath)
    %code for initialising parallel computing
    
    c = parcluster('local'); % build the 'local' cluster object
    numCore = c.NumWorkers;        % get the number of workers
    saveProfile(c);
    
    display(num2str(numCore))
    
%     parpool('local',numCore,'SpmdEnabled',false);
%     fileInfo     = dir([PreprocessDir animalNum '*']);
end

theseRecs = [29:136]; %29:136 all trials with bilateral eye-tracking
numRecs = numel(theseRecs);

for irec =  1:numRecs%1:numel(fileInfo) %12:15%
    
    rec2Analyze = theseRecs(irec);
    
    recName = fileInfo(rec2Analyze).name;
    
    recPath = [PreprocessDir recName '\'];
    rootAnalysisDir = [AnalysisDir recName '\Pupil\']; 
    %is_detectSaccadesEphys(recPath)
    is_pupilState_2eyes(recPath)

    load([recPath 'pupilStates']);
    [pupDiv,eyeNames,binLims] = is_load([recPath 'pupilStates'],'pupDiv','eyeNames','binLims');
    
    fig = figure();
    
    for iEye = 1:size(pupDiv,1)
        subplot(1,2,iEye)
        cpup = pupDiv(iEye,:);
        minMaxPup = [min(cpup,[],2) max(cpup,[],2)]; %nEye x 2
        [h,c] = hist(cpup,minMaxPup(1):0.006:minMaxPup(2));
        bm       = blueMap; %64 rows
        cls      = flip(bm(1:10:end,:),1); % from light to dark
    
        %binLims(7,2) = 2.3305;
        for ibin = 1:5
            ind = (c>=binLims(iEye,ibin,1) & c<= binLims(iEye,ibin,2));
            bar(c(ind),h(ind),0.9,'FaceColor',cls(ibin,:),'EdgeColor','none'); hold on
        end
        title([eyeNames{iEye} ' dia distribution' ]);
        if iEye == 1; legend('Low','2','3','4','High','Location','northwest'); ylabel('Binned count');end
        xlabel('Pupil diameter [a.u]')        
        %xlim([minMaxPup(1) minMaxPup(2)])
        %ylim([0 20000])
        set(gca,'TickDir','out')
    end
    
    savefig(fig, [rootAnalysisDir 'Pupil_DiaDistribution.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'Pupil_DiaDistribution.png']); 
    close all
end