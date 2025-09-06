function stimActivityAnalysis(spikeData, Params, Info, figFolder, oneFigureHandle)
% Performs analysis of spike data relative to stimulation times 
% INPUT
% -------
% spikeData : struct 
% Params : struct 
% Info : struct 
% figFolder : path 
% oneFigureHandle : matlab figure object

    
    %% Gather stimulation times
    allStimTimes = [];
    
    for channelIdx = 1:length(spikeData.stimInfo)
        
        allStimTimes = [allStimTimes, spikeData.stimInfo{channelIdx}.elecStimTimes];
    
    end

    %% Firing rate before and after stimulation 
    figName = '9_FR_before_after_stimulation';
    plotPrePostStimFR(spikeData, allStimTimes, Params, figFolder, figName, Info);

    % Do it for each pattern 
    for patternIdx = 1:length(spikeData.stimPatterns)
        figName = sprintf('9_FR_before_after_stimulation_pattern_%.f', patternIdx);
        plotPrePostStimFR(spikeData, spikeData.stimPatterns{patternIdx}, Params, figFolder, figName, Info);
    end


    %% Plot population raster align to stim times 
    % rasterWindow = [Params.preStimWindow(1), Params.postStimWindow(2)];
    Params.rasterBinWidth = 0.01; %  TODO: move this to the app

    [frAlignedToStim, rasterBins] = getFrAlignedToStim(spikeData, allStimTimes, Params);
    numChannels = size(frAlignedToStim, 1);

    figureHandle = figure;
    figName = '10_stimulation_raster_and_psth';
    ylabel_txt = 'Mean firing rate (spikes/s)';
    plotMetricAlignedToStim(frAlignedToStim, rasterBins, Info, Params, ...
    ylabel_txt, figFolder, figName, figureHandle)
    
     % Do it for each pattern 
    for patternIdx = 1:length(spikeData.stimPatterns)
        figureHandle = figure;
        figName = sprintf('10_stimulation_raster_and_psth_pattern_%.f', patternIdx);
        [frAlignedToStim, rasterBins] = getFrAlignedToStim( ... 
            spikeData, spikeData.stimPatterns{patternIdx}, Params);
        plotMetricAlignedToStim(frAlignedToStim, rasterBins, Info, Params, ...
                ylabel_txt, figFolder, figName, figureHandle)
    end


    %% Look at spike amplitude aligned to stimulus 
    % TODO: Loop throgh patterns
    numStimEvent = length(allStimTimes);
    spikeAmps = getSpikeAmp(spikeData.spikeWaveforms); 
    spikeData.spikeAmps = spikeAmps;
    
    rasterWindow = [Params.preStimWindow(1), Params.postStimWindow(2)];
    rasterBinWidth = Params.rasterBinWidth;   % originally 0.025 
    
    rasterBins = rasterWindow(1):rasterBinWidth:rasterWindow(2);
    numBins = length(rasterBins) - 1; 
    
    ampAlignedToStim = zeros(numChannels, numStimEvent, numBins) + nan;

    for channelIdx= 1:numChannels
        channelSpikeTimes = spikeData.spikeTimes{channelIdx}.(Params.SpikesMethod);
        channelSpikeAmps = spikeData.spikeAmps{channelIdx}.(Params.SpikesMethod);
        
        % Process spike times to remove spikes near stimulus time 
        
        for stimTimeIdx = 1:length(allStimTimes)
            stimTime = allStimTimes(stimTimeIdx);
            removeIndex = find((channelSpikeTimes >=  stimTime + Params.stimRemoveSpikesWindow(1)) & ...
                               (channelSpikeTimes <=  stimTime + Params.stimRemoveSpikesWindow(2)));
            channelSpikeTimes(removeIndex) = [];
            channelSpikeAmps(removeIndex) = [];
        end  
         
    
         for stimEventIdx = 1:numStimEvent 
            stimTime = allStimTimes(stimEventIdx);
            
            for binIdx = 1:length(rasterBins)-1
                binStart = stimTime + rasterBins(binIdx);
                binEnd = stimTime + rasterBins(binIdx+1);
                spikeIdx = find((channelSpikeTimes >= binStart)  & (channelSpikeTimes < binEnd));
                meanAmps = mean(abs(channelSpikeAmps(spikeIdx)));
                ampAlignedToStim(channelIdx, stimEventIdx, binIdx) = meanAmps;
            end
    
            % ampAlignedToStim(channelIdx, stimEventIdx, :) = histcounts(channelSpikeTimes - stimTime, rasterBins) / rasterBinWidth;
    
         end 
    
        
    end
    
    figureHandle = figure;
    set(figureHandle, 'Position', [100, 100, 600, 500]);
    subplot(2, 1, 1)
    meanAmpalignedToStim = squeeze(nanmean(ampAlignedToStim, [1, 2]));
    plot(rasterBins(2:end), meanAmpalignedToStim)
    hold on 
    fill([Params.stimRemoveSpikesWindow(1), Params.stimRemoveSpikesWindow(2), ...
          Params.stimRemoveSpikesWindow(2), Params.stimRemoveSpikesWindow(1)], ...
         [0, 0, max(meanAmpalignedToStim), max(meanAmpalignedToStim)], [0.5, 0.5, 0.5], 'FaceAlpha', 0.3,'LineStyle','none')
    box off 
    set(gca, 'TickDir', 'out');
    ylabel('Mean absolute spike amplitude')
    title(Info.FN{1}, 'Interpreter', 'none');
    
    subplot(2, 1, 2)
    imagesc(rasterBins(2:end), 1:numChannels, squeeze(nanmean(ampAlignedToStim, 2)))
    box off
    ylabel('Channel')
    xlabel('Time from stimulation (s)')
    set(gca, 'TickDir', 'out');
    set(gcf, 'color', 'w');
    

    % save figure
    figName = '11_stimulation_amplitude_raster_and_psth';
    pipelineSaveFig(fullfile(figFolder, figName), Params.figExt, Params.fullSVG, gcf);

    %% Spike Latency / Time-to-first-spike 
    
    numPatterns = length(spikeData.stimPatterns);
    channelMeanSpikeLatency = zeros(numChannels, numPatterns) + nan;
    % NOTE: Here we remove spikes around each stimulus time, regardless
    % of which pattern, hence the use of "allStimTimes(:)" to flatten it
    spikeTimesStimRemoved = removeStimSpikes(allStimTimes(:), spikeData.spikeTimes, Params);
    for patternIdx = 1:numPatterns
        stimTimesToAlign = spikeData.stimPatterns{patternIdx};
        % for each electrode
        for channelIdx= 1:numChannels
            channelSpikeTimes = spikeTimesStimRemoved{channelIdx}.(Params.SpikesMethod);
            spikeLatencies = getSpikeLatencyRelStim(stimTimesToAlign, channelSpikeTimes);
            channelMeanSpikeLatency(channelIdx, patternIdx) = nanmean(spikeLatencies);
        end
    end 
    % TODO: Make null distribution of spike latency

    % Make plots of spike latency
    vrange = [min(channelMeanSpikeLatency(:)) max(channelMeanSpikeLatency(:))];
    cmap = flip(viridis);
    for patternIdx = 1:length(spikeData.stimPatterns)
        oneFigureHandle = figure();
        nodeMetric = channelMeanSpikeLatency(:, patternIdx);
        cmapLabel = 'Mean spike latency (s)';
        oneFigureHandle = plotStimHeatmapWmetric(nodeMetric, vrange, cmap, cmapLabel, ...
            spikeData.stimInfo, patternIdx, oneFigureHandle);
        title(sprintf('Pattern %.f', patternIdx));
        figName = sprintf('spikeLatency_pattern_%.f_heatmap', patternIdx);
        figPath = fullfile(figFolder, figName);
        pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, oneFigureHandle)
    end

    %% Plot heatmap of activity after stimulation 
    % stim_activity_window = [0.1, 0.3];  % seconds
    patternSpikeMatrixStore = {};
    stimActivityStore = {};

    % some parameter adjustments to look into smaller windows 
    % Params.stimRemoveSpikesWindow = [0, 0]
    % Params.rasterBinWidth = 0.001; % originally 0.01
    
    Params.stimDecodingTimeWindows = linspace(0, 0.01, 10);
    % Params.stimDecodingTimeWindows = [0, 0.002, 0.004, 0.006, 0.008, 0.01];

    for patternIdx = 1:length(spikeData.stimPatterns)
        stimTimesToAlign = spikeData.stimPatterns{patternIdx};
        [frAlignedToStim, rasterBins] = getFrAlignedToStim(spikeData, stimTimesToAlign, Params);
        
        numDecodingTimeBins = length(Params.stimDecodingTimeWindows)-1;
        numTrials = size(frAlignedToStim, 2);
        patternSpikeMatrix = zeros(numTrials, numChannels*numDecodingTimeBins);
        for window_idx = 1:numDecodingTimeBins
            subsetTimeIdx = find((rasterBins >= Params.stimDecodingTimeWindows(window_idx)) & ...
            (rasterBins <= Params.stimDecodingTimeWindows(window_idx+1))) - 1;
            patternSpikeMatrixAtWindow = mean(frAlignedToStim(:, :, subsetTimeIdx), 3); % mean across time [numChannels, numTrials]
            start_idx = 1 + (window_idx - 1) * numChannels; 
            end_idx = start_idx + numChannels - 1;
            patternSpikeMatrix(:, start_idx:end_idx) = patternSpikeMatrixAtWindow'; % Transpose to match [numTrials, numChannels]
        end

        % subsetTimeIdx = find((rasterBins >= Params.postStimWindow(1)) & ...
        %     (rasterBins <= Params.postStimWindow(2))) - 1;
        % patternSpikeMatrix = mean(frAlignedToStim(:, :, subsetTimeIdx), 3); % mean across time
        stimActivityStore{patternIdx} = patternSpikeMatrix; % for decoding 
        patternSpikeMatrix = squeeze(mean(patternSpikeMatrix, 1)); % mean across trials
       
        patternSpikeMatrixStore{patternIdx} = patternSpikeMatrix;
        
        % set stimulating electrode value to 0??? 

        % get spike matrix from alignment
        % electrodeHeatMaps(FN, spikeMatrix, channels, spikeFreqMax, Params, coords, figFolder, oneFigureHandle);
    end
    
    % get the max frequency to scale 
    spikeFreqMax = max(cellfun(@max, patternSpikeMatrixStore));
    spikeFreqMax = max([spikeFreqMax, 0.0001]);
    vrange = [0, spikeFreqMax];
    cmap = 'viridis';
    
    for patternIdx = 1:length(spikeData.stimPatterns)
        oneFigureHandle = figure();
        nodeMetric = patternSpikeMatrixStore{patternIdx};
        cmapLabel = 'Firing rate (spikes/s)';
        oneFigureHandle = plotStimHeatmapWmetric(nodeMetric, vrange, cmap, cmapLabel, ...
            spikeData.stimInfo, patternIdx, oneFigureHandle);
        title(sprintf('Pattern %.f', patternIdx));
        figName = sprintf('stimPattern_%.f_heatmap', patternIdx);
        figPath = fullfile(figFolder, figName);
        pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, oneFigureHandle)
    end


    %% Plot trial by trial stimulation response 
    oneFigureHandle = figure;
    t = tiledlayout(1, length(stimActivityStore), ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    all_activity = vertcat(stimActivityStore{:});
    min_activity = min(all_activity(:));
    max_activity = max(all_activity(:));
    
    if max_activity == 0 
        max_activity = 0.0001;
    end

    clim = [min_activity max_activity];
    for stimId = 1:length(stimActivityStore)
        
        ax = nexttile;
        imagesc(stimActivityStore{stimId}, clim)
        patternNumTrials = size(stimActivityStore{stimId}, 1);
        
        xlabel('Electrodes and time bins')
        ylabel('Trials')
        hTitle = title(['Stim pattern ' num2str(stimId)], '   ');  % add an empty line to create padding
        
        
        hold on
        for window_idx = 1:numDecodingTimeBins
            start_idx = 1 + (window_idx-1) * numChannels; 
            end_idx = start_idx + numChannels - 1;
            plot([start_idx, end_idx], [-0.5, -0.5], 'LineWidth', 3, 'Color', 'black');
            decodingStartMs = Params.stimDecodingTimeWindows(window_idx) * 1000;
            text(start_idx, -0.5, sprintf('%.f ms', decodingStartMs), 'VerticalAlignment', 'bottom');
            
            if window_idx == numDecodingTimeBins
                decodingEndMs = Params.stimDecodingTimeWindows(window_idx+1) * 1000;
                text(end_idx, -0.5, sprintf('%.f ms', decodingEndMs), ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment','left');
            end
            
            hold on 
        end

        ylim([-1, patternNumTrials]);
        yticks([1, patternNumTrials])
        box(ax, 'off');
        set(ax, 'TickDir', 'out');  

    end 
    cb = colorbar(ax, 'eastoutside');
    cb.Layout.Tile = 'east';  % Put the colorbar outside the tile layout
    ylabel(cb, 'spikes/s', 'FontSize',12)
    set(gcf, 'color', 'w')
    
    figName = 'stimPattern_activity_per_trial';
    figPath = fullfile(figFolder, figName);
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, oneFigureHandle)


    %% Do decoding of which pattern was stimulated
    numNodesToTry = 1:5:numChannels;
    numRepeatsPerNodeNumber = 10;
    decodingAccuracy = zeros(length(numNodesToTry), numRepeatsPerNodeNumber);

    % construct our X and y for decoding 
    % TODO: adjust which features to use (ie. include spike latency)
    X = vertcat(stimActivityStore{:});  % numTrials x numNodes
    
    % TEMP: Do some z-scoring 
    X = (X - mean(X, 1)) ./ std(X, 1);

    % Do some NaN imputation... 
    X(isnan(X)) = 0;

    y = []; 
    for patternIdx = 1:length(stimActivityStore)
        y = [y; repmat(patternIdx, size(stimActivityStore{patternIdx}, 1), 1)];
    end 
    
    numTrials = size(X, 1);
    % trialIdx = 1:numTrials;
    % propTrainTrials = 0.5;  % proportion of trials for training set
    % numTrainTrials = floor(propTrainTrials * numTrials);
    clf_num_kfold = 5;

    for numNodeIdx = 1:length(numNodesToTry)
        numNodesToUse = numNodesToTry(numNodeIdx);
        
        for repeatIdx = 1:numRepeatsPerNodeNumber
            
            % randomly get a subset of nodes (without replacement)
            nodeToUse = randsample(numChannels, numNodesToUse, false);
            featureIndices = [];
            for windowIdx = 1:numDecodingTimeBins
                startIdx = (windowIdx-1) * numChannels;
                featureIndices = [featureIndices; startIdx + nodeToUse];
            end

            % trialIdxShuffled = trialIdx(randperm(length(trialIdx));
            % train_idx = trialIdx(1:numTrainTrials);
            %test_idx = trialIdx((numTrainTrials+1):end);
            X_subset = X(:, featureIndices);

            if all(isnan(X_subset(:)))
                mean_model_loss = nan;
            else
                clf_model = fitcecoc(X_subset,y);  %multi-class model
                % clf_model = fitcsvm(X_subset,y);
                cross_val_model = crossval(clf_model, 'KFold', clf_num_kfold);
                mean_model_loss = mean(cross_val_model.kfoldLoss);
            end 
            decodingAccuracy(numNodeIdx, repeatIdx) = 1 - mean_model_loss;
        end

    end

    oneFigureHandle = figure;
    plot(numNodesToTry, mean(decodingAccuracy, 2));
    hold on 
    scatter(numNodesToTry, mean(decodingAccuracy, 2));
    ylim([0, 1])
    xlabel('Number of nodes used in classification')
    ylabel('Classification accuracy')
    set(gcf, 'color', 'w')
    figName = 'stimPattern_decoding';
    figPath = fullfile(figFolder, figName);
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, oneFigureHandle)

    %% Plot individual PSTHs
    % Integrated PSTH analysis for individual electrodes per stimulation pattern
    % Based on batch_psth_baseline_analysis.m functionality
    
    % PSTH Analysis Parameters
    psth_window_s = [0, 0.02];           % 20ms analysis window post-stimulus
    psth_bin_width_s = 0.001;            % 1ms bin width for PSTH
    num_baseline_psths = 100;            % Number of baseline PSTHs 
    baseline_duration_s = psth_window_s(2) - psth_window_s(1);  % Match analysis window duration
    % NEW APPROACH: Calculate artifact window based on blank durations from detectStimTimesTemplate
    % plus a hardcoded postBlankIgnore value
    postBlankIgnore = 0.5; % Hardcoded to 0.5 ms additional time after blank end
    
    % Get blank durations from any channel that has it (should be same across channels)
    allBlankDurations = [];
    for ch_idx = 1:length(spikeData.stimInfo)
        if isfield(spikeData.stimInfo{ch_idx}, 'allBlankDurations') && ...
           ~isempty(spikeData.stimInfo{ch_idx}.allBlankDurations)
            allBlankDurations = spikeData.stimInfo{ch_idx}.allBlankDurations;
            break;
        end
    end
    
    if isempty(allBlankDurations)
        error('allBlankDurations is required but not found in stimInfo. Ensure detectStimTimesTemplate has been run.');
    end
    
    % Calculate artifact window as mode of blank durations plus postBlankIgnore
    mode_blank_duration_ms = mode(allBlankDurations * 1000); % Convert to ms
    artifact_window_end_ms = mode_blank_duration_ms + postBlankIgnore;
    artifact_window_ms = [0, artifact_window_end_ms]; % Window from stimTime to mode+postBlankIgnore
    
    fprintf('Using blank duration-based artifact window: [0, %.2f] ms (mode duration: %.2f ms + %.2f ms post-blank ignore)\n', ...
            artifact_window_end_ms, mode_blank_duration_ms, postBlankIgnore);
    
    % Loop through each stimulation pattern
    for patternIdx = 1:length(spikeData.stimPatterns)
        stimTimes = spikeData.stimPatterns{patternIdx};  % Get stim times for this pattern
        
        if isempty(stimTimes)
            continue;  % Skip if no stimulation times for this pattern
        end
        
        % Create pattern-specific output folder
        patternFolderName = sprintf('Individual_PSTH_and_Raster_Pattern_%d', patternIdx);
        patternFigFolder = fullfile(figFolder, patternFolderName);
        if ~exist(patternFigFolder, 'dir')
            mkdir(patternFigFolder);
        end
        
        % Initialize storage for this pattern's results
        networkResponse = [];
        valid_channel_count = 0;
        
        % Loop through each channel
        for channelIdx = 1:numChannels
            % Exclude stimulated channels for the current pattern
            if isfield(spikeData, 'stimInfo') && ...
               channelIdx <= length(spikeData.stimInfo) && ...
               isfield(spikeData.stimInfo{channelIdx}, 'pattern') && ...
               ~isempty(spikeData.stimInfo{channelIdx}.pattern) && ...
               spikeData.stimInfo{channelIdx}.pattern == patternIdx
                continue; % Skip analysis for stimulated channels
            end

            % Extract spike times for current channel
            if channelIdx > length(spikeData.spikeTimes) || ...
               isempty(spikeData.spikeTimes{channelIdx}) || ...
               ~isfield(spikeData.spikeTimes{channelIdx}, Params.SpikesMethod)
                continue;  % Skip if no spike data
            end
            
            all_spike_times_s = spikeData.spikeTimes{channelIdx}.(Params.SpikesMethod);
            if isempty(all_spike_times_s)
                continue;  % Skip if no spikes
            end
            
            % Remove spikes within artifact window for each stimulus
            spikeTimes_cleaned_s = all_spike_times_s;
            for stimIdx = 1:length(stimTimes)
                stimTime = stimTimes(stimIdx);
                spikeTimes_cleaned_s = spikeTimes_cleaned_s(...
                    spikeTimes_cleaned_s < (stimTime + artifact_window_ms(1)/1000) | ...
                    spikeTimes_cleaned_s >= (stimTime + artifact_window_ms(2)/1000));
            end
            spikeTimes_cleaned_s = sort(spikeTimes_cleaned_s(:));
            
            % Calculate response PSTH
            [response, resp_metrics] = calculate_psth_metrics(...
                spikeTimes_cleaned_s, stimTimes, psth_window_s, psth_bin_width_s);
            
            if isempty(response.psth_samples)
                continue;  % Skip if no spikes in response window
            end
            
            % Calculate multiple baseline PSTHs
            baseline_aucs = zeros(num_baseline_psths, 1);
            all_baseline_psth_smooth = [];
            
            for i = 1:num_baseline_psths
                % Define baseline window moving backwards from stimulus
                start_s = -(i * baseline_duration_s);
                end_s = -((i-1) * baseline_duration_s);
                current_baseline_window_s = [start_s, end_s];
                
                % Apply artifact blanking to baseline window
                spikeTimes_cleaned_baseline_s = spikeTimes_cleaned_s;
                blank_start_offset_s = start_s + artifact_window_ms(1)/1000;
                blank_end_offset_s = start_s + artifact_window_ms(2)/1000;
                
                for stimIdx = 1:length(stimTimes)
                    stimTime = stimTimes(stimIdx);
                    spikeTimes_cleaned_baseline_s = spikeTimes_cleaned_baseline_s(...
                        spikeTimes_cleaned_baseline_s < (stimTime + blank_start_offset_s) | ...
                        spikeTimes_cleaned_baseline_s >= (stimTime + blank_end_offset_s));
                end
                
                [~, base_metrics] = calculate_psth_metrics(...
                    spikeTimes_cleaned_baseline_s, stimTimes, current_baseline_window_s, psth_bin_width_s);
                baseline_aucs(i) = base_metrics.auc;
                
                if isempty(all_baseline_psth_smooth)
                    all_baseline_psth_smooth = zeros(num_baseline_psths, length(base_metrics.psth_smooth));
                end
                all_baseline_psth_smooth(i, :) = base_metrics.psth_smooth;
            end
            
            % Calculate baseline-corrected metrics
            mean_baseline_auc = mean(baseline_aucs);
            auc_corrected = resp_metrics.auc - mean_baseline_auc;
            mean_baseline_psth = mean(all_baseline_psth_smooth, 1);
            
            % Decay analysis using half-max from peak
            [Rmax, Rmax_idx] = max(resp_metrics.psth_smooth);
            halfRmax = Rmax / 2;
            halfRmax_idx = find(resp_metrics.psth_smooth(Rmax_idx:end) <= halfRmax, 1, 'first');
            
            if ~isempty(halfRmax_idx)
                halfRmax_idx = halfRmax_idx + Rmax_idx - 1;
                halfRmax_time_s = resp_metrics.time_vector_s(halfRmax_idx);
                halfRmax_val = resp_metrics.psth_smooth(halfRmax_idx);
            else
                halfRmax_time_s = NaN;
                halfRmax_val = NaN;
            end
            
            % Generate PSTH plot for this channel
            fig = figure('Position', [100 100 1200 900], 'Visible', 'off');
            psth_window_ms = psth_window_s * 1000;
            
            % Get channel ID (use electrode info if available, otherwise use index)
            if isfield(spikeData, 'stimInfo') && channelIdx <= length(spikeData.stimInfo) && ...
               isfield(spikeData.stimInfo{channelIdx}, 'channelName')
                channel_id = spikeData.stimInfo{channelIdx}.channelName;
            else
                channel_id = channelIdx;  % Fallback to channel index
            end
            
            sgtitle(sprintf('Pattern %d | Channel %d | Peak Rate: %.1f Hz | Corrected AUC: %.3f', ...
                patternIdx, channel_id, resp_metrics.peak_firing_rate, auc_corrected), 'FontWeight', 'bold');
            
            % Subplot 1 (top right): Spike Raster Plot
            ax1 = subplot(2, 2, 2); hold on;
            for trial_idx = 1:length(response.spikeTimes_byEvent)
                trial_spikes_s = response.spikeTimes_byEvent{trial_idx};
                if ~isempty(trial_spikes_s)
                    plot(trial_spikes_s * 1000, trial_idx * ones(size(trial_spikes_s)), ...
                        'r.', 'MarkerSize', 5);
                end
            end
            hold off;
            set(gca, 'YDir', 'reverse');
            xlim(psth_window_ms);
            ylim([0 length(stimTimes)+1]);
            ylabel('Trial Number');
            xlabel('Time from stimulus (ms)');
            title('Spike Raster (Response)');
            grid on;
            
            % Subplot 2 (bottom right): Response vs Baselines comparison
            ax2 = subplot(2, 2, 4); hold on;
            baseline_time_ms = (base_metrics.time_vector_s - current_baseline_window_s(1)) * 1000;
            for i = 1:num_baseline_psths
                plot(baseline_time_ms, all_baseline_psth_smooth(i, :), ...
                    'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
            end
            p1_diag = plot(resp_metrics.time_vector_s*1000, resp_metrics.psth_smooth, ...
                'r-', 'LineWidth', 2);
            p2_diag = plot(baseline_time_ms, mean_baseline_psth, 'k-', 'LineWidth', 2);
            hold off;
            title('Diagnostic: Response vs. Baselines');
            ylabel('Firing Rate (spikes/s)');
            xlabel('Time from stimulus (ms)');
            legend([p1_diag, p2_diag], 'Response', 'Mean Baseline', 'Location', 'Best');
            grid on;
            
            % Subplot 3 (left): Smoothed response PSTH with metrics
            ax3 = subplot(2, 2, [1, 3]); hold on;
            edges_s = psth_window_s(1):psth_bin_width_s:psth_window_s(2);
            bar(edges_s * 1000, response.psth_histogram, 1, ...
                'FaceColor', 0.7*[1 1 1], 'EdgeColor', 0.8*[1 1 1], 'HandleVisibility', 'off');
            plot(resp_metrics.time_vector_s * 1000, resp_metrics.psth_smooth, ...
                'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed PSTH');
            plot(resp_metrics.peak_time_s*1000, resp_metrics.peak_firing_rate, ...
                'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'R_{max}');
            text(resp_metrics.peak_time_s*1000, resp_metrics.peak_firing_rate, ...
                sprintf(' R_{max}: %.1f Hz', resp_metrics.peak_firing_rate), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
                'Color', 'k', 'FontWeight', 'bold');
            
            if ~isnan(halfRmax_time_s)
                plot(halfRmax_time_s*1000, halfRmax_val, ...
                    'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'Half R_{max}');
                text(halfRmax_time_s*1000, halfRmax_val, ...
                    sprintf(' Half R_{max} @ %.1f ms', halfRmax_time_s*1000), ...
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
                    'Color', 'k', 'FontWeight', 'bold');
            end
            hold off;
            ylabel('Firing Rate (spikes/s)');
            xlabel('Time from stimulus (ms)');
            title('Smoothed Response PSTH & Metrics');
            legend('Location', 'northeast');
            grid on;
            
            linkaxes([ax1, ax2, ax3], 'x');
            xlim(psth_window_ms);
            
            % Save plot
            plot_filename = fullfile(patternFigFolder, sprintf('Individual_PSTH_and_Raster_channel_%d.png', channel_id));
            pipelineSaveFig(plot_filename, Params.figExt, Params.fullSVG, fig);
            close(fig);
            
            % Store results for this channel
            valid_channel_count = valid_channel_count + 1;
            networkResponse(valid_channel_count).channel_id = channel_id;
            networkResponse(valid_channel_count).file_index = channelIdx;
            networkResponse(valid_channel_count).pattern_id = patternIdx;
            networkResponse(valid_channel_count).auc_response = resp_metrics.auc;
            networkResponse(valid_channel_count).auc_baseline_mean = mean_baseline_auc;
            networkResponse(valid_channel_count).auc_corrected = auc_corrected;
            networkResponse(valid_channel_count).peak_firing_rate_hz = resp_metrics.peak_firing_rate;
            networkResponse(valid_channel_count).peak_time_ms = resp_metrics.peak_time_s * 1000;
            networkResponse(valid_channel_count).halfRmax_time_ms = halfRmax_time_s * 1000;
        end
        
        % Save networkResponse data for this pattern
        if ~isempty(networkResponse)
            timestamp = datestr(now, 'ddmmmyyyy_HHMMSS');
            output_filename = fullfile(patternFigFolder, ...
                sprintf('networkResponse_pattern_%d_%s.mat', patternIdx, timestamp));
            save(output_filename, 'networkResponse');
        end
    end

end
