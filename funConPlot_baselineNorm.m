baseTwin = [-2.5,-0.5];
doMedian = 1;

% Mozart, 002 session
tvec = is_load('D:\FerretData\0172\Analyzed\0172_MozartAudio_002_20190204\FC_validChns\LPl-PPC\specAll_lrAudio.mat','tvec');
avgPLV = is_load('D:\FerretData\0172\Analyzed\0172_MozartAudio_002_20190204\FC_validChns\LPl-PPC\funcCon_avg_lrAudio.mat', 'avgPLV');
[tvecGC,avgGC_XtoY,avgGC_YtoX] = is_load('D:\FerretData\0172\Analyzed\0172_MozartAudio_002_20190204\FC_validChns\LPl-PPC\GCAll_lrAudio.mat','tvecGC','avgGC_XtoY','avgGC_YtoX');
[normedPLV_Mozart,~,~] = baselineNorm(avgPLV, tvec, baseTwin, doMedian);
[normedGC_Mozart_LPl_PPC,~,~] = baselineNorm(real(avgGC_XtoY), tvecGC, baseTwin, doMedian);
[normedGC_Mozart_PPC_LPl,~,~] = baselineNorm(real(avgGC_YtoX), tvecGC, baseTwin, doMedian);

% Whitenoise, 007 session
avgPLV = is_load('D:\FerretData\0172\Analyzed\0172_MozartAudio_007_20190212\FC_firstChn\LPl-PPC\funcCon_avg_whitenoise.mat', 'avgPLV');
[tvecGC,avgGC_XtoY,avgGC_YtoX] = is_load('D:\FerretData\0172\Analyzed\0172_MozartAudio_007_20190212\FC_firstChn\LPl-PPC\GCAll_whitenoise.mat','tvecGC','avgGC_XtoY','avgGC_YtoX');
[normedPLV_WN,~,~] = baselineNorm(avgPLV, tvec, baseTwin, doMedian);
[normedGC_WN_LPl_PPC,~,~] = baselineNorm(real(avgGC_XtoY), tvecGC, baseTwin, doMedian);
[normedGC_WN_PPC_LPl,~,~] = baselineNorm(real(avgGC_YtoX), tvecGC, baseTwin, doMedian);



regionXname = 'LPl';
regionYname = 'PPC';

%% plot
lowFreq = 2; highFreq = 128; numFreqs = 150;
foi     = logspace(log10(lowFreq),log10(highFreq),numFreqs);
fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
tickLabel = string(fois);

for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end

% Plot
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100
Ylim = [1,100]; %1-30Hz
% PLV
subplot(2,3,1)
imagesc(tvec,1:numel(foi),normedPLV_Mozart);
%imagesc(tvec,foi,avgPLV);
xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
ylim(Ylim);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
caxis([-3,6]);
cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)

% plot GC
subplot(2,3,2)
imagesc(tvecGC,1:numel(foi),normedGC_Mozart_LPl_PPC);
xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12)
%caxis([0 0.3]); 
ylim(Ylim); % 90+ Hz has saturated values

subplot(2,3,3)
imagesc(tvecGC,1:numel(foi),normedGC_Mozart_PPC_LPl);
xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
%caxis([0 0.3]);
ylim(Ylim); % 90+ Hz has saturated values

% whitenoise
subplot(2,3,4)
imagesc(tvec,1:numel(foi),normedPLV_WN);
%imagesc(tvec,foi,avgPLV);
xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
ylim(Ylim);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%caxis([0.1 0.7]);
cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)

% plot GC
subplot(2,3,5)
imagesc(tvecGC,1:numel(foi),normedGC_WN_LPl_PPC);
xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12)
%caxis([0 0.3]);
ylim(Ylim); % 90+ Hz has saturated values


subplot(2,3,6)
imagesc(tvecGC,1:numel(foi),normedGC_WN_PPC_LPl);
xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
%caxis([0 0.3]);
ylim(Ylim); % 90+ Hz has saturated values

colormap(jet)
GroupAnalysisDir = 'D:\FerretData\0172\GroupAnalysis\forGrant\'; 
savefig(fig, [GroupAnalysisDir 'normedFC_Mozart+WN_2.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'normedFC_Mozart+WN_2.png']);
        