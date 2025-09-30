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

            % trialIdxShuffled = trialIdx(randperm(length(trialIdx)));
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
    num_baseline_psths = 50;            % Number of baseline PSTHs 
    baseline_duration_s = psth_window_s(2) - psth_window_s(1);  % Match analysis window duration
    
    % PSTH Smoothing Parameters - MANUALLY SPECIFY HERE
    psth_smoothing_method = 'gaussian';  % Options: 'ssvkernel' or 'gaussian'
    psth_gaussian_width_ms = 1;          % Gaussian kernel width in ms (only used if method is 'gaussian')
    
    fprintf('PSTH Smoothing: %s', psth_smoothing_method);
    if ~strcmp(psth_smoothing_method, 'ssvkernel')
        fprintf(' (Gaussian width: %.1f ms)', psth_gaussian_width_ms);
    end
    fprintf('\n');
    % NEW APPROACH: Calculate artifact window based on blank durations from detectStimTimesTemplate
    % plus a hardcoded postBlankIgnore value
    postBlankIgnore = 0.5; % Hardcoded to 0.5 ms additional time after blank end
    
    % Get blank durations from any channel that has it (should be same across channels)
    blankDurations = [];
    for ch_idx = 1:length(spikeData.stimInfo)
        if isfield(spikeData.stimInfo{ch_idx}, 'blankDurations') && ...
           ~isempty(spikeData.stimInfo{ch_idx}.blankDurations)
            blankDurations = spikeData.stimInfo{ch_idx}.blankDurations;
            break;
        end
    end
    
    if isempty(blankDurations)
        error('blankDurations is required but not found in stimInfo. Ensure detectStimTimesTemplate has been run with longblank method.');
    end
    
    % Calculate artifact window as mode of blank durations plus postBlankIgnore
    mode_blank_duration_ms = mode(blankDurations * 1000); % Convert to ms
    artifact_window_end_ms = mode_blank_duration_ms + postBlankIgnore;
    artifact_window_ms = [0, artifact_window_end_ms]; % Window from stimTime to mode+postBlankIgnore
    
    fprintf('Using blank duration-based artifact window: [0, %.2f] ms (mode duration: %.2f ms + %.2f ms post-blank ignore)\n', ...
            artifact_window_end_ms, mode_blank_duration_ms, postBlankIgnore);
    
    % Create global list of stimulated channels to exclude from ALL pattern analyses
    % Only exclude channels with pattern > 0 (pattern 0 = no stimulation)
    stimulatedChannels = [];
    if isfield(spikeData, 'stimInfo')
        for channelIdx = 1:length(spikeData.stimInfo)
            if isfield(spikeData.stimInfo{channelIdx}, 'pattern') && ...
               ~isempty(spikeData.stimInfo{channelIdx}.pattern) && ...
               spikeData.stimInfo{channelIdx}.pattern > 0
                stimulatedChannels = [stimulatedChannels, channelIdx];
            end
        end
    end
    stimulatedChannels = unique(stimulatedChannels);  % Remove duplicates
    fprintf('Excluding %d stimulated channels from analysis across all patterns: [%s]\n', ...
            length(stimulatedChannels), num2str(stimulatedChannels));
    
    % Create consolidated stimulation times and pattern labels for all patterns
    allStimTimesConsolidated = [];
    stimPatternLabels = [];
    
    for patternIdx = 1:length(spikeData.stimPatterns)
        stimTimes = spikeData.stimPatterns{patternIdx};
        if ~isempty(stimTimes)
            allStimTimesConsolidated = [allStimTimesConsolidated, stimTimes];
            stimPatternLabels = [stimPatternLabels, repmat(patternIdx, 1, length(stimTimes))];
        end
    end
    
    % Sort stimulation times chronologically and keep track of pattern labels
    [allStimTimesConsolidated, sortIdx] = sort(allStimTimesConsolidated);
    stimPatternLabels = stimPatternLabels(sortIdx);
    
    % Initialize matrix for average firing rates in post-stim window
    % Dimensions: [numTrialsTotal x numChannels] where numTrialsTotal = all stim times across patterns
    numTrialsTotal = length(allStimTimesConsolidated);
    avgFiringRateMatrix = NaN(numTrialsTotal, numChannels);
    
    fprintf('Creating consolidated firing rate matrix: %d trials x %d channels\n', numTrialsTotal, numChannels);
    
    % Loop through each stimulation pattern for individual PSTH analysis
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
        
        % First pass: collect all channel metrics to identify top channels for plotting
        channel_metrics = [];
        temp_data = {};
        
        % Loop through each channel
        for channelIdx = 1:numChannels
            % Exclude any channels that are stimulated in ANY pattern
            if ismember(channelIdx, stimulatedChannels)
                continue; % Skip analysis for all stimulated channels
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
            
            % Calculate response PSTH with smoothing parameters
            if strcmp(psth_smoothing_method, 'ssvkernel')
                [response, resp_metrics] = calculate_psth_metrics(...
                    spikeTimes_cleaned_s, stimTimes, psth_window_s, psth_bin_width_s, ...
                    'smoothing_method', 'ssvkernel');
            else
                [response, resp_metrics] = calculate_psth_metrics(...
                    spikeTimes_cleaned_s, stimTimes, psth_window_s, psth_bin_width_s, ...
                    'smoothing_method', 'gaussian', 'gaussian_width_ms', psth_gaussian_width_ms);
            end
            
            if isempty(response.psth_samples)
                continue;  % Skip if no spikes in response window
            end
            
            % Calculate firing rates for each trial in post-stim window (for consolidated matrix)
            % Only calculate once per channel across all trials
            if patternIdx == 1  % Calculate consolidated matrix only once per channel
                % Duration of analysis window in seconds
                analysis_window_duration_s = psth_window_s(2) - psth_window_s(1);
                
                for trialIdx = 1:numTrialsTotal
                    stimTime = allStimTimesConsolidated(trialIdx);
                    
                    % Count spikes in post-stim window for this trial (excluding artifact window)
                    trial_window_start = stimTime + psth_window_s(1);
                    trial_window_end = stimTime + psth_window_s(2);
                    
                    % Remove artifact window from spike counting
                    artifact_start = stimTime + artifact_window_ms(1)/1000;
                    artifact_end = stimTime + artifact_window_ms(2)/1000;
                    
                    % Find spikes in analysis window but outside artifact window
                    spikes_in_window = all_spike_times_s(...
                        (all_spike_times_s >= trial_window_start & all_spike_times_s <= trial_window_end) & ...
                        ~(all_spike_times_s >= artifact_start & all_spike_times_s < artifact_end));
                    
                    % Calculate firing rate for this trial and electrode
                    num_spikes = length(spikes_in_window);
                    effective_window_duration = analysis_window_duration_s - (artifact_window_ms(2) - artifact_window_ms(1))/1000;
                    firing_rate_hz = num_spikes / effective_window_duration;
                    
                    avgFiringRateMatrix(trialIdx, channelIdx) = firing_rate_hz;
                end
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
                
                % Calculate baseline PSTH with same smoothing parameters
                if strcmp(psth_smoothing_method, 'ssvkernel')
                    [~, base_metrics] = calculate_psth_metrics(...
                        spikeTimes_cleaned_baseline_s, stimTimes, current_baseline_window_s, psth_bin_width_s, ...
                        'smoothing_method', 'ssvkernel');
                else
                    [~, base_metrics] = calculate_psth_metrics(...
                        spikeTimes_cleaned_baseline_s, stimTimes, current_baseline_window_s, psth_bin_width_s, ...
                        'smoothing_method', 'gaussian', 'gaussian_width_ms', psth_gaussian_width_ms);
                end
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
            
            % Calculate d-prime for stimulus response significance
            % First calculate trial-by-trial firing rates for baseline and post-stim periods
            baseline_window_s_dprime = [-psth_window_s(2), -psth_window_s(1)]; % Same duration as post-stim, but before
            baseline_firing_rates_dprime = [];
            poststim_firing_rates_dprime = [];
            
            for stimIdx = 1:length(stimTimes)
                stimTime = stimTimes(stimIdx);
                
                % Baseline firing rate for this trial (with artifact exclusion for consistency)
                baseline_start = stimTime + baseline_window_s_dprime(1);
                baseline_end = stimTime + baseline_window_s_dprime(2);
                % Apply artifact blanking to baseline window (for consistency with post-stim)
                artifact_start_baseline = stimTime + baseline_window_s_dprime(1) + artifact_window_ms(1)/1000;
                artifact_end_baseline = stimTime + baseline_window_s_dprime(1) + artifact_window_ms(2)/1000;
                baseline_spikes = all_spike_times_s(...
                    (all_spike_times_s >= baseline_start & all_spike_times_s < baseline_end) & ...
                    ~(all_spike_times_s >= artifact_start_baseline & all_spike_times_s < artifact_end_baseline));
                baseline_duration = baseline_window_s_dprime(2) - baseline_window_s_dprime(1) - (artifact_window_ms(2) - artifact_window_ms(1))/1000;
                baseline_firing_rates_dprime(stimIdx) = length(baseline_spikes) / baseline_duration;
                
                % Post-stim firing rate for this trial (with artifact exclusion)
                poststim_start = stimTime + psth_window_s(1);
                poststim_end = stimTime + psth_window_s(2);
                poststim_spikes = all_spike_times_s(...
                    (all_spike_times_s >= poststim_start & all_spike_times_s < poststim_end) & ...
                    ~(all_spike_times_s >= (stimTime + artifact_window_ms(1)/1000) & ...
                      all_spike_times_s < (stimTime + artifact_window_ms(2)/1000)));
                poststim_duration = psth_window_s(2) - psth_window_s(1) - (artifact_window_ms(2) - artifact_window_ms(1))/1000;
                poststim_firing_rates_dprime(stimIdx) = length(poststim_spikes) / poststim_duration;
            end
            
            % Calculate statistics for d-prime
            baseline_mean_hz_dprime = mean(baseline_firing_rates_dprime);
            baseline_std_hz_dprime = std(baseline_firing_rates_dprime);
            poststim_mean_hz_dprime = mean(poststim_firing_rates_dprime);
            poststim_std_hz_dprime = std(poststim_firing_rates_dprime);
            
            % Calculate d-prime (d') - signal detection theory measure
            % d' = (μ_post - μ_baseline) / √((σ²_post + σ²_baseline) / 2)
            if baseline_std_hz_dprime == 0 && poststim_std_hz_dprime == 0
                d_prime = abs(poststim_mean_hz_dprime - baseline_mean_hz_dprime);
            else
                pooled_variance = (baseline_std_hz_dprime^2 + poststim_std_hz_dprime^2) / 2;
                d_prime = (poststim_mean_hz_dprime - baseline_mean_hz_dprime) / sqrt(pooled_variance);
            end
            
            % Store channel data for later processing
            valid_channel_count = valid_channel_count + 1;
            
            % Get channel ID (use electrode info if available, otherwise use index)
            if isfield(spikeData, 'stimInfo') && channelIdx <= length(spikeData.stimInfo) && ...
               isfield(spikeData.stimInfo{channelIdx}, 'channelName')
                channel_id = spikeData.stimInfo{channelIdx}.channelName;
            else
                channel_id = channelIdx;  % Fallback to channel index
            end
            
            % Store all necessary data for this channel
            temp_data{valid_channel_count}.channelIdx = channelIdx;
            temp_data{valid_channel_count}.channel_id = channel_id;
            temp_data{valid_channel_count}.auc_corrected = auc_corrected;
            temp_data{valid_channel_count}.response = response;
            temp_data{valid_channel_count}.resp_metrics = resp_metrics;
            temp_data{valid_channel_count}.base_metrics = base_metrics;
            temp_data{valid_channel_count}.all_baseline_psth_smooth = all_baseline_psth_smooth;
            temp_data{valid_channel_count}.mean_baseline_auc = mean_baseline_auc;
            temp_data{valid_channel_count}.mean_baseline_psth = mean_baseline_psth;
            temp_data{valid_channel_count}.halfRmax_time_s = halfRmax_time_s;
            temp_data{valid_channel_count}.halfRmax_val = halfRmax_val;
            temp_data{valid_channel_count}.current_baseline_window_s = current_baseline_window_s;
            temp_data{valid_channel_count}.d_prime = d_prime;
            temp_data{valid_channel_count}.baseline_firing_rates_dprime = baseline_firing_rates_dprime;
            temp_data{valid_channel_count}.poststim_firing_rates_dprime = poststim_firing_rates_dprime;
            
            % Store results for networkResponse
            networkResponse(valid_channel_count).channel_id = channel_id;
            networkResponse(valid_channel_count).file_index = channelIdx;
            networkResponse(valid_channel_count).pattern_id = patternIdx;
            networkResponse(valid_channel_count).auc_response = resp_metrics.auc;
            networkResponse(valid_channel_count).auc_baseline_mean = mean_baseline_auc;
            networkResponse(valid_channel_count).auc_corrected = auc_corrected;
            networkResponse(valid_channel_count).peak_firing_rate_hz = resp_metrics.peak_firing_rate;
            networkResponse(valid_channel_count).peak_time_ms = resp_metrics.peak_time_s * 1000;
            networkResponse(valid_channel_count).halfRmax_time_ms = halfRmax_time_s * 1000;
            networkResponse(valid_channel_count).d_prime = d_prime;
            networkResponse(valid_channel_count).baseline_mean_hz = baseline_mean_hz_dprime;
            networkResponse(valid_channel_count).baseline_std_hz = baseline_std_hz_dprime;
            networkResponse(valid_channel_count).poststim_mean_hz = poststim_mean_hz_dprime;
            networkResponse(valid_channel_count).poststim_std_hz = poststim_std_hz_dprime;
        end
        
        % Second pass: identify top 5 channels with corrected AUC > 0.5 for plotting
        if ~isempty(temp_data)
            % Extract corrected AUC values
            auc_values = [networkResponse.auc_corrected];
            
            % Find channels with corrected AUC > 0.5
            high_auc_indices = find(auc_values > 0.5);
            
            % Sort by corrected AUC (descending) and take top 5
            [~, sort_indices] = sort(auc_values(high_auc_indices), 'descend');
            top_channels_for_plotting = high_auc_indices(sort_indices(1:min(5, length(sort_indices))));
            
            fprintf('Pattern %d: Found %d channels with corrected AUC > 0.5. Plotting top %d channels.\n', ...
                patternIdx, length(high_auc_indices), length(top_channels_for_plotting));
            
            % Generate plots only for selected channels
            for plot_idx = 1:length(top_channels_for_plotting)
                channel_data_idx = top_channels_for_plotting(plot_idx);
                data = temp_data{channel_data_idx};
                
                % Generate PSTH plot for this channel
                fig = figure('Position', [100 100 1200 900], 'Visible', 'off');
                psth_window_ms = psth_window_s * 1000;
                
                sgtitle(sprintf('Pattern %d | Channel %d | Peak Rate: %.1f Hz | Corrected AUC: %.3f', ...
                    patternIdx, data.channel_id, data.resp_metrics.peak_firing_rate, data.auc_corrected), 'FontWeight', 'bold');
                
                % Subplot 1 (top right): Spike Raster Plot
                ax1 = subplot(2, 2, 2); hold on;
                for trial_idx = 1:length(data.response.spikeTimes_byEvent)
                    trial_spikes_s = data.response.spikeTimes_byEvent{trial_idx};
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
                baseline_time_ms = (data.base_metrics.time_vector_s - data.current_baseline_window_s(1)) * 1000;
                for i = 1:num_baseline_psths
                    plot(baseline_time_ms, data.all_baseline_psth_smooth(i, :), ...
                        'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
                end
                p1_diag = plot(data.resp_metrics.time_vector_s*1000, data.resp_metrics.psth_smooth, ...
                    'r-', 'LineWidth', 2);
                p2_diag = plot(baseline_time_ms, data.mean_baseline_psth, 'k-', 'LineWidth', 2);
                hold off;
                title('Diagnostic: Response vs. Baselines');
                ylabel('Firing Rate (spikes/s)');
                xlabel('Time from stimulus (ms)');
                legend([p1_diag, p2_diag], 'Response', 'Mean Baseline', 'Location', 'Best');
                grid on;
                
                % Subplot 3: Smoothed response PSTH with metrics and dual y-axes
                ax3 = subplot(2, 2, [1, 3]); hold on;
                
                % Use already calculated d-prime statistics for z-score PSTH analysis
                % Handle zero std case (if all baseline trials have identical firing rates)
                baseline_std_hz_safe = baseline_std_hz_dprime;
                if baseline_std_hz_safe == 0
                    baseline_std_hz_safe = eps; % Use machine epsilon to avoid division by zero
                end
                
                % Calculate z-score for the smoothed PSTH using d-prime baseline statistics
                zscore_psth = (data.resp_metrics.psth_smooth - baseline_mean_hz_dprime) ./ baseline_std_hz_safe;
                
                % Calculate trial-based z-score for reference using d-prime statistics
                trial_zscore = (poststim_mean_hz_dprime - baseline_mean_hz_dprime) / baseline_std_hz_safe;
                
                fprintf('Ch %d: Baseline=%.1f±%.1f Hz (n=%d trials), PostStim=%.1f±%.1f Hz, Z=%.1f, d''=%.2f, Max PSTH Z=%.1f\n', ...
                    data.channel_id, baseline_mean_hz_dprime, baseline_std_hz_dprime, length(baseline_firing_rates_dprime), ...
                    poststim_mean_hz_dprime, poststim_std_hz_dprime, trial_zscore, d_prime, max(zscore_psth));
                
                % DEBUG: Check if spike data is being shared between channels
                fprintf('  DEBUG Ch %d: Total spikes=%d, Unique spike times (first 5): [%.3f, %.3f, %.3f, %.3f, %.3f]\n', ...
                    data.channel_id, length(all_spike_times_s), ...
                    all_spike_times_s(min(1,end)), all_spike_times_s(min(2,end)), all_spike_times_s(min(3,end)), ...
                    all_spike_times_s(min(4,end)), all_spike_times_s(min(5,end)));
                fprintf('  Individual trial baseline rates: [%.1f, %.1f, %.1f, %.1f, %.1f] Hz (showing first 5)\n', ...
                    baseline_firing_rates_dprime(min(1,end)), baseline_firing_rates_dprime(min(2,end)), baseline_firing_rates_dprime(min(3,end)), ...
                    baseline_firing_rates_dprime(min(4,end)), baseline_firing_rates_dprime(min(5,end)));
                fprintf('  Individual trial post-stim rates: [%.1f, %.1f, %.1f, %.1f, %.1f] Hz (showing first 5)\n', ...
                    poststim_firing_rates_dprime(min(1,end)), poststim_firing_rates_dprime(min(2,end)), poststim_firing_rates_dprime(min(3,end)), ...
                    poststim_firing_rates_dprime(min(4,end)), poststim_firing_rates_dprime(min(5,end)));
                
                % Debug: Print baseline window details
                fprintf('  Baseline window: [%.1f, %.1f] ms, Post-stim: [%.1f, %.1f] ms, Artifact: [%.1f, %.1f] ms\n', ...
                    baseline_window_s_dprime(1)*1000, baseline_window_s_dprime(2)*1000, ...
                    psth_window_s(1)*1000, psth_window_s(2)*1000, ...
                    artifact_window_ms(1), artifact_window_ms(2));
                
                % Plot Z-score PSTH only
                p_psth_zscore = plot(data.resp_metrics.time_vector_s * 1000, zscore_psth, ...
                    'r-', 'LineWidth', 2, 'DisplayName', 'Z-score PSTH');
                
                % Calculate z-score equivalent of peak and half-max markers
                [zscore_peak_val, zscore_peak_idx] = max(zscore_psth);
                zscore_peak_time_ms = data.resp_metrics.time_vector_s(zscore_peak_idx) * 1000;
                
                % Plot peak marker in z-score units
                plot(zscore_peak_time_ms, zscore_peak_val, ...
                    'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'HandleVisibility', 'off');
                text(zscore_peak_time_ms, zscore_peak_val, ...
                    sprintf(' Z_{max}: %.1f (%.1f Hz)', zscore_peak_val, data.resp_metrics.peak_firing_rate), ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
                    'Color', 'k', 'FontWeight', 'bold');
                
                % Calculate half-max in z-score units
                zscore_halfmax = zscore_peak_val / 2;
                zscore_halfmax_idx = find(zscore_psth(zscore_peak_idx:end) <= zscore_halfmax, 1, 'first');
                
                if ~isempty(zscore_halfmax_idx)
                    zscore_halfmax_idx = zscore_halfmax_idx + zscore_peak_idx - 1;
                    zscore_halfmax_time_ms = data.resp_metrics.time_vector_s(zscore_halfmax_idx) * 1000;
                    zscore_halfmax_val = zscore_psth(zscore_halfmax_idx);
                    
                    plot(zscore_halfmax_time_ms, zscore_halfmax_val, ...
                        'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'HandleVisibility', 'off');
                    text(zscore_halfmax_time_ms, zscore_halfmax_val, ...
                        sprintf(' Half Z_{max} @ %.1f ms', zscore_halfmax_time_ms), ...
                        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
                        'Color', 'k', 'FontWeight', 'bold');
                end
                
                ylabel('Z-score');
                xlabel('Time from stimulus (ms)');
                title('Z-score PSTH & Metrics');
                legend('Z-score PSTH', 'Location', 'northeast');
                grid on;
                
                linkaxes([ax1, ax2, ax3], 'x');
                xlim(psth_window_ms);
                
                % Save plot using pipeline parameters
                figName = sprintf('Individual_PSTH_and_Raster_channel_%d', data.channel_id);
                figPath = fullfile(patternFigFolder, figName);
                pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, fig);
                
                close(fig);
            end
        end
        
        % Save networkResponse data for this pattern
        if ~isempty(networkResponse)
            timestamp = datestr(now, 'ddmmmyyyy_HHMMSS');
            output_filename = fullfile(patternFigFolder, ...
                sprintf('networkResponse_pattern_%d_%s.mat', patternIdx, timestamp));
            save(output_filename, 'networkResponse');
        end
    end
    
    % Save consolidated average firing rate matrix (outside the pattern loop)
    % Matrix dimensions: [numTrialsTotal x numChannels] 
    % Values: Average firing rate (Hz) in post-stim window with artifact exclusion
    % NaN values indicate stimulated channels (excluded from analysis)
    avgFiringRateMatrix_info = struct();
    avgFiringRateMatrix_info.description = 'Consolidated average firing rates (Hz) in post-stimulus window across all patterns';
    avgFiringRateMatrix_info.dimensions = '[numTrialsTotal x numChannels]';
    avgFiringRateMatrix_info.analysis_window_s = psth_window_s;
    avgFiringRateMatrix_info.artifact_window_ms = artifact_window_ms;
    avgFiringRateMatrix_info.stimulated_channels_excluded = stimulatedChannels;
    avgFiringRateMatrix_info.allStimTimesConsolidated = allStimTimesConsolidated;
    avgFiringRateMatrix_info.stimPatternLabels = stimPatternLabels;
    avgFiringRateMatrix_info.numTrialsTotal = numTrialsTotal;
    
    % Create main output folder for consolidated data
    consolidatedFolder = fullfile(figFolder, 'Reservoir Computing Metrics');
    if ~exist(consolidatedFolder, 'dir')
        mkdir(consolidatedFolder);
    end
    
    matrix_filename = fullfile(consolidatedFolder, 'avgFiringRateMatrix_consolidated.mat');
    save(matrix_filename, 'avgFiringRateMatrix', 'avgFiringRateMatrix_info');
    
    fprintf('Saved consolidated firing rate matrix to: %s\n', matrix_filename);

end
