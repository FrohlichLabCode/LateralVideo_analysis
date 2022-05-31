function regionPair_FunConn(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs, sessionID,...
    evtTimes,twin,baseTwin, condNames, condID, regionPairs, regionNames, ...
    regionLFP, regionChn, rootAnalysisDir, GroupAnalysisDir)

cluster = 0;
doPlot = 1;
doMean = 0;
doMedian = 1;
numConds  = numel(condID);

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if     linORlog == 1
        foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
        foi      = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end    

for iPair = 1:numel(regionPairs)
    regionXind = regionPairs{iPair}(1);
    regionYind = regionPairs{iPair}(2);
    regionXname = regionNames{regionPairs{iPair}(1)};
    regionYname = regionNames{regionPairs{iPair}(2)};

    
    saveAnalysisDir    = [rootAnalysisDir regionXname '-' regionYname '/']; %FC/Pul-PPC/
    if ~exist(join(GroupAnalysisDir),'dir'); mkdir(join(GroupAnalysisDir));end
    if ~exist(join(saveAnalysisDir),'dir'); mkdir(join(saveAnalysisDir));
%    elseif length(dir([saveAnalysisDir '*.fig'])) >= length(condID) %each condition generate a fig file
%        fprintf('Record %s already analyzed \n',saveAnalysisDir'); 
%        if skipRec == 1; continue; end  % to skip already analyzed records
    end
    fprintf('\nWorking on %s \n',saveAnalysisDir');     
    
    regionxLFP = regionLFP{regionXind};
    regionyLFP = regionLFP{regionYind};
    numChnReg1 = numel(regionChn{regionXind});
    numChnReg2 = numel(regionChn{regionYind});

    for iCond = 1:numConds
        condName = condNames{condID(iCond)};
        if skipRec==1 && (length(dir([saveAnalysisDir '*FC_' condName '*.fig'])) >= (doMean+doMedian))
            continue; end % skip already processed condition             
        evtTime = evtTimes{condID(iCond)};     
        %twin    = twins{condID(iCond)};
        %baseTwin = baseTwins{condID(iCond)};
        
        % pre-allocate size for speed
        numT = (twin(2)-twin(1))*lfpFs/10+1; % depends on the suppression ratio in is_functionalConnectivity
        numTGC = (twin(2)-twin(1))*10+1; % depends on the suppression ratio in is_functionalConnectivity
        numFOI = length(foi);
        tvecAll = nan(numChnReg1,numChnReg2,numT);
        xSpecAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        ySpecAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        xSpecNormedAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        ySpecNormedAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        plvAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        coherencyAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        coherenceAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        iCoherenceAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        imagZAll = nan(numChnReg1,numChnReg2,numFOI,numT);
        %psiAll = nan(numChnReg1,numChnReg2,numFOI,numT); %size is not numFOI
        %psiNormAll = nan(numChnReg1,numChnReg2,numFOI,numT); %size is not numFOI
        
        tvecGCAll = nan(numChnReg1,numChnReg2,numTGC);
        GC_XtoY_All = nan(numChnReg1,numChnReg2,numFOI,numTGC);
        GC_YtoX_All = nan(numChnReg1,numChnReg2,numFOI,numTGC);

%         if cluster == 0; parforArg = 0;else parforArg = Inf; end %flag for whether use parfor or for
%         parfor (iReg1Chn = 1:numChnReg1, parforArg)
        for iReg1Chn = 1:numChnReg1
            xser = regionxLFP(iReg1Chn,:);

            for iReg2Chn = 1:numChnReg2
                fprintf('\nWorking on pair %d/%d, OuterCounter=%d, InnerCounter=%d ================================ \n',...
                    iReg2Chn+(iReg1Chn-1)*numChnReg2 , numChnReg1*numChnReg2, iReg1Chn, iReg2Chn);

                yser = regionyLFP(iReg2Chn,:);
                
                %clear funcCon % to clear out memory
                funcCon = is_functionalConnectivity_V2(xser,yser,regionXname,regionYname,condNames{condID(iCond)},...
                    lfpFs,evtTime,twin,baseTwin,regionChn{regionXind}(iReg1Chn),regionChn{regionYind}(iReg2Chn),...
                    saveAnalysisDir, linORlog, lowFreq, highFreq, numFreqs);
                close all
            
                % for each struct matrix, it's frequency by time
                tvecAll(iReg1Chn,iReg2Chn,:) = funcCon.tvec;
                psiFreqAll(iReg1Chn,iReg2Chn,:) = funcCon.psiFreq;
                try
                tvecGCAll(iReg1Chn,iReg2Chn,:) = funcCon.grangerCausality.tvec;
                catch
                end
                xSpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.xspec;
                ySpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.yspec;
                xSpecNormedAll(iReg1Chn,iReg2Chn,:,:) = funcCon.xspecNormed;
                ySpecNormedAll(iReg1Chn,iReg2Chn,:,:) = funcCon.yspecNormed;
                %only need to change iReg2Chn
                %but Subscripted assignment dimension mismatch.for
                %if iReg1Chn == 1;  ySpecAll(iReg2Chn,:,:) = funcCon.yspec;
                plvAll(iReg1Chn,iReg2Chn,:,:) = funcCon.plv;
                coherencyAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherency; % might want to keep the complex part to look at phase lags
                coherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherence;
                iCoherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imaginaryCoherence;
                imagZAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imagZ;
                psiAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psi;
                psiNormAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psiNorm;
                try
                GC_XtoY_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.X_to_Y;
                GC_YtoX_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.Y_to_X;
                catch
                end                
            end % 
        end
        % only need one instance, all rows should be the same
        tvec = squeeze(tvecAll(1,1,:))';
        psiFreq = squeeze(psiFreqAll(1,1,:))';
        tvecGC = squeeze(tvecGCAll(1,1,:))';
        
        
    %% Eliminate nan sessions
    % The results from some sessions have nan elements. Those sessions need to
    % be eliminated. Note that simple use of "nanmean" does not work at this stage, but will be used later (EN).

    for i=1:size(xSpecAll, 1)
        for j=1:size(xSpecAll, 2)
            temp1 = squeeze(xSpecAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                xSpecAll(i, j, :, :) = nan;
            end
            clear temp1 temp2


            temp1 = squeeze(ySpecAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in ySpec   ', i, j);
                ySpecAll(i, j, :, :) = nan;
            end
            clear temp1 temp2

    %         temp1 = squeeze(GC_XtoY_All(i, j, :, :));
    %         temp2 = isnan(temp1);
    %         if any(temp2(:))
    %             fprintf('\n[i, j]= [%d, %d], NaN in GC_XtoY   ', i, j);
    %             GC_XtoY_All(i, j, :, :) = nan;
    %         end
    %         clear temp1 temp2
    %         
    %         temp1 = squeeze(GC_YtoX_All(i, j, :, :));
    %         temp2 = isnan(temp1);
    %         if any(temp2(:))
    %             fprintf('\n[i, j]= [%d, %d], NaN in GC_YtoX   ', i, j);
    %             GC_YtoX_All(i, j, :, :) = nan;
    %         end
    %         clear temp1 temp2

            temp1 = squeeze(plvAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in plv   ', i, j);
                plvAll(i, j, :, :) = nan;
            end
            clear temp1 temp2

            temp1 = squeeze(coherencyAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in coherency   ', i, j);
                coherencyAll(i, j, :, :) = nan;
            end
            clear temp1 temp2

            temp1 = squeeze(imagZAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in coherency   ', i, j);
                imagZAll(i, j, :, :) = nan;
            end
            clear temp1 temp2


        end
    end

    %% Compute the mean values across channel pairs

    % Compute ticks for plotting
    if linORlog == 1
        fois = [2, 5:5:highFreq];
        tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
        psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]); % generate a string array for phase slope index
    elseif linORlog == 2
        fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
        tickLabel = string(fois);
        psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]);
    end
    for fi = 1:numel(fois)
        [bi,bb] = sort(abs(foi-fois(fi)));
        tickLoc(fi) = bb(1);
        [bi,bb] = sort(abs(psiFreq-fois(fi)));
        psitickLoc(fi) = bb(1);
    end

    if doMean == 1
    % calculating mean
    avgXSpec     = squeeze(nanmean(nanmean( xSpecAll ,1),2)); % change mean to nanmean
    avgYSpec     = squeeze(nanmean(nanmean( ySpecAll ,1),2));
    avgXNormed   = squeeze(nanmean(nanmean( xSpecNormedAll ,1),2));
    avgYNormed   = squeeze(nanmean(nanmean( ySpecNormedAll ,1),2));
    avgPLV       = squeeze(nanmean(nanmean( plvAll ,1),2));
    avgCoherency = squeeze(nanmean(nanmean( coherencyAll ,1),2));
    avgImagZ     = squeeze(nanmean(nanmean( imagZAll ,1),2));
    avgpsiNorm   = squeeze(nanmean(nanmean( psiNormAll ,1),2));
    try
    avgGC_XtoY   = squeeze(nanmean(nanmean( GC_XtoY_All ,1),2));
    avgGC_YtoX   = squeeze(nanmean(nanmean( GC_YtoX_All ,1),2));
    catch
    end

    
    cd([saveAnalysisDir]);
    % save means
    save(['funcCon_avg_' condName '.mat'],'tvec','fois','avgXSpec','avgYSpec','avgXNormed', 'avgYNormed','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');
    fprintf('\nDone saving mean');
    end
    
    % save all channel data
    cd([saveAnalysisDir]);
    
    save(['specAll_' condName '.mat'],'tvec','foi','xSpecAll','ySpecAll','-v7.3');
    save(['plvAll_' condName '.mat'],'plvAll','-v7.3');
    save(['coherencyAll_' condName '.mat'],'coherencyAll','-v7.3');
    save(['psiAll_' condName '.mat'],'psiAll','psiFreq','psiNormAll','-v7.3');
    try
    save(['GCAll_' condName '.mat'],'tvecGC','GC_XtoY_All','GC_YtoX_All','avgGC_XtoY','avgGC_YtoX','-v7.3');
    catch
    end
    fprintf('\nDone saving ============================================\n')
    
    if doMedian == 1
    % calculate median
    avgXSpec     = squeeze(nanmedian(nanmedian( xSpecAll ,1),2)); % change mean to nanmean
    avgYSpec     = squeeze(nanmedian(nanmedian( ySpecAll ,1),2));
    avgXNormed   = squeeze(nanmedian(nanmedian( xSpecNormedAll ,1),2));
    avgYNormed   = squeeze(nanmedian(nanmedian( ySpecNormedAll ,1),2));
    avgPLV       = squeeze(nanmedian(nanmedian( plvAll ,1),2));
        tempReal = nanmedian(nanmedian(real(coherencyAll),1),2);
        tempImag = 1i*nanmedian(nanmedian(imag(coherencyAll),1),2);
    avgCoherency = squeeze(tempReal + tempImag); 
        tempReal = nanmedian(nanmedian(real(imagZAll),1),2);
        tempImag = 1i*nanmedian(nanmedian(imag(imagZAll),1),2); 
    avgImagZ     = squeeze(tempReal + tempImag);         
    avgpsiNorm   = squeeze(nanmedian(nanmedian( psiNormAll ,1),2));
    try
        tempReal = nanmedian(nanmedian(real(GC_XtoY_All),1),2);
        tempImag = 1i*nanmedian(nanmedian(imag(GC_XtoY_All),1),2);     
    avgGC_XtoY   = squeeze(tempReal + tempImag); 
        tempReal = nanmedian(nanmedian(real(GC_YtoX_All),1),2);
        tempImag = 1i*nanmedian(nanmedian(imag(GC_YtoX_All),1),2);    
    avgGC_YtoX   = squeeze(tempReal + tempImag); 
    catch
    end    
    
    cd([saveAnalysisDir]);
    save(['funcCon_median_' condName '.mat'],'tvec','fois','avgXSpec','avgYSpec','avgXNormed', 'avgYNormed','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');
    try
    save(['GC_medain_' condName '.mat'],'tvecGC','avgGC_XtoY','avgGC_YtoX','-v7.3');
    catch
    end
    fprintf('\nDone saving median ============================================\n')
    end
    
    % to save memory space and also avoid using old data
    clear xSpecAll ySpecAll xSpecNormedAll ySpecNormedAll 
    clear plvAll coherencyAll coherenceAll iCoherenceAll imagZAll psiAll psiNormAll GC_XtoY_All GC_YtoX_All 


    %% plotting
    if doPlot == 1
        %% first plot based on median 
        
        % Compute ticks for plotting
        if linORlog == 1
            fois = [2, 5:5:highFreq];
            tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
            psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]); % generate a string array for phase slope index
        elseif linORlog == 2
            fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
            tickLabel = string(fois);
            psitickLabel = string([round(psiFreq(1)) fois(2:end-1) round(psiFreq(end))]);
        end
        for fi = 1:numel(fois)
            [bi,bb] = sort(abs(foi-fois(fi)));
            tickLoc(fi) = bb(1);
            [bi,bb] = sort(abs(psiFreq-fois(fi)));
            psitickLoc(fi) = bb(1);
        end

        % Plot
        screensize = get( groot, 'Screensize' );
        fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100
        if doMedian == 1
        % plot power spectrum for signal x
        subplot(3,4,1)
        try
        imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
    %    imagesc(tvec,foi,pow2db(avgXSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12)
        catch
        end
        % plot power spectrum for signal y
        try
        subplot(3,4,2)
        imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12)
        catch
        end
        try  
        subplot(3,4,3)
        imagesc(tvec,1:numel(foi),avgXNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
        catch
        end
        % plot power spectrum for signal y
        try
        % plot power spectrum for signal y
        subplot(3,4,4)
        imagesc(tvec,1:numel(foi),avgYNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
        catch
        end
        try
        % plot phase locking value
        subplot(3,4,5)
        imagesc(tvec,1:numel(foi),avgPLV);
        %imagesc(tvec,foi,avgPLV);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.1 0.7]);
        cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
        catch
        end

        try
        % plot coherence
        subplot(3,4,6)
        imagesc(tvec,1:numel(foi),abs(avgCoherency));
        %imagesc(tvec,foi,abs(avgCoherencey));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([0 0.6]);
        cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
        catch
        end

        try
        % plot imaginary coherence
        subplot(3,4,7)
        imagesc(tvec,1:numel(foi),abs(avgImagZ));
        %imagesc(tvec,foi,imag(avgCoherency));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Imaginary coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 4]);
        cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
        % plot phase slope index
        catch
        end

        try
            subplot(3,4,8)
            imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
            %imagesc(tvec,psiFreq,avgpsiNorm);
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Phase slope index')
            ylim([psitickLoc(1) psitickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
            % plot granger causality X to Y
            caxis([-4 4])
        catch
        end
        try
            subplot(3,4,9)
            imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
            %imagesc(tvecGC,foi,real(avgGC_XtoY));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
            % plot granger causality Y to X
            subplot(3,4,10)
            imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
            %imagesc(tvecGC,foi,real(avgGC_YtoX));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: Y to X')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
        catch
        end
        colormap(jet)

        savefig(fig, ['medianFC_' condName '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz.fig'],'compact');
        saveas(fig, ['medianFC_' condName '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz.png']);
        saveas(fig, [GroupAnalysisDir 'medianFC_' sessionID '_' regionXname '-' regionYname '_' condName '.png']);        

        close all
        clear avgXSpec avgYSpec avgXNormed avgYNormed
        clear avgPLV avgCoherency avgImagZ avgpsiNorm avgGC_XtoY avgGC_YtoX
        end
        
        %% Then plot based on mean
        if doMean == 1
        load(['funcCon_avg_' condName '.mat']);
        try
        load(['GCAll_' condName '.mat'],'avgGC_XtoY','avgGC_YtoX');
        catch
        end
        % Plot
        screensize = get( groot, 'Screensize' );
        fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100

        % plot power spectrum for signal x
        try
        subplot(3,4,1)
        imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
    %    imagesc(tvec,foi,pow2db(avgXSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12)
        catch
        end
        
        % plot power spectrum for signal y
        try
        subplot(3,4,2)
        imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12)
        catch
        end
        
        try        
        subplot(3,4,3)
        imagesc(tvec,1:numel(foi),avgXNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
        catch
        end
        
        % plot power spectrum for signal y
        try
        subplot(3,4,4)
        imagesc(tvec,1:numel(foi),avgYNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
        catch
        end
        
        try
        % plot phase locking value
        subplot(3,4,5)
        imagesc(tvec,1:numel(foi),avgPLV);
        %imagesc(tvec,foi,avgPLV);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.1 0.7]);
        cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
        catch
        end

        try
        % plot coherence
        subplot(3,4,6)
        imagesc(tvec,1:numel(foi),abs(avgCoherency));
        %imagesc(tvec,foi,abs(avgCoherencey));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([0 0.6]);
        cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
        catch
        end

        try
        % plot imaginary coherence
        subplot(3,4,7)
        imagesc(tvec,1:numel(foi),abs(avgImagZ));
        %imagesc(tvec,foi,imag(avgCoherency));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Imaginary coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 4]);
        cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
        % plot phase slope index
        catch
        end

        try
            subplot(3,4,8)
            imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
            %imagesc(tvec,psiFreq,avgpsiNorm);
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Phase slope index')
            ylim([psitickLoc(1) psitickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
            % plot granger causality X to Y
            caxis([-4 4])
        catch
        end
        try
            subplot(3,4,9)
            imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
            %imagesc(tvecGC,foi,real(avgGC_XtoY));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
            % plot granger causality Y to X
            subplot(3,4,10)
            imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
            %imagesc(tvecGC,foi,real(avgGC_YtoX));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: Y to X')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
        catch
        end
        colormap(jet)
        cd([saveAnalysisDir]);
        savefig(fig, ['meanFC_' condName '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz.fig'],'compact');
        saveas(fig, ['meanFC_' condName '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz.png']);
        saveas(fig, [GroupAnalysisDir 'meanFC_' sessionID '_' regionXname '-' regionYname '_' condName '.png']);        

        close all
        clear avgXSpec avgYSpec avgXNormed avgYNormed
        clear avgPLV avgCoherency avgImagZ avgpsiNorm avgGC_XtoY avgGC_YtoX        
        end
        
        end
    end
end