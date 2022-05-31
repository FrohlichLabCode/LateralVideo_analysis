load('D:\FerretData\0172\Analyzed\0172_MozartAudio_007_20190212\Pupil\Pupil_baselineNormed_median.mat')
rootAnalysisDir = 'D:\FerretData\0172\Analyzed\0172_MozartAudio_007_20190212\Pupil\';
numPul = 6;
eyeChn = {'leftX','rightX','leftY','rightY','leftD','rightD'};
creensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2*numPul/3 (screensize(4)-150)/2]);
condNames = {'Mozart','Bach','pureC','whitenoise'};
condID = [1,3,4];
numConds = numel(condID);
legendName = {};
iPlot = 1;
for iChn = [5,6] %1:size(evtDat_trialAvg,2)
    
    subplot(1,numPul/3,iPlot)
    hold on
    
    for iCond = numConds:-1:1
        iCondID = condID(iCond);
        toPlot = squeeze(evtDat_trialAvg(iCondID,iChn,:));
        toPlot = smooth(squeeze(evtDat_trialAvg(iCondID,iChn,:)),1000); %lowpass filter
        %toPlot = smoothdata(toPlot,'gaussian',200); % not much difference
        plot(timeVec,toPlot, 'LineWidth', 4)
        legendName{end+1} = condNames{iCondID};
    end
    if iPlot ==1; legend(legendName);end
    title([eyeChn{iChn}])
    xlabel('Time [s]');
    ylabel('Normed amplitude');
%     nXtick = round((twin(2)-twin(1))/2)+1;
%     set(gca,'XTick',linspace(1,numWinsamp,nXtick))
%     set(gca,'XTickLabel',linspace(twin(1),twin(2),nXtick))    
    hline(0,'k-');
    axis tight
    iPlot = iPlot +1;
    xlim([-1.9,7.9]);
    ylim([-1 2]);
    vline(0,'k-');vline(5,'k--');
end
if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
savefig(fig, [rootAnalysisDir 'Pupil_baselineNormed_median_smooth_new.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Pupil_baselineNormed_median_smooth_new.png']);
