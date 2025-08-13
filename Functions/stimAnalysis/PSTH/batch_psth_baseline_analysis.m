% -------------------------------------------------------------------------
% BATCH PSTH ANALYSIS SCRIPT
% -------------------------------------------------------------------------
%
% Description:
% This script performs a PSTH-based analysis of evoked response for 
% multiple pairs of spike data and raw recording files. For each
% file pair, it detects stimulation events, calculates the PSTH response
% and baseline firing statistics for each channel, generates summary plots,
% and saves the quantitative results.
%
% Inputs:
% 1.  `filePairs` (cell array of structs): User-defined list of files to
%     analyze. Each struct must contain two fields:
%     - `spikeFile`: Path to the .mat file with spike times.
%     - `rawFile`: Path to the corresponding raw data .mat file.
% 2.  User-defined analysis parameters within the script, such as sampling
%     frequency (`fs`), window sizes, etc.
%
% Outputs:
% For each processed file pair, a new directory is created with the format:
% "PSTH_Analysis_[spike_file_basename]_[timestamp]"
% This directory will contain:
% 1.  `Analysis_Plot_Chan_[ID].png`: A plot for each
%     analyzed channel, showing the spike raster and smoothed PSTHs
% 2.  `networkResponse_all_channels_[timestamp].mat`: A .mat file containing a
%     struct array (`networkResponse`) with saved values of metrics for all valid
%     channels.
%
% Dependencies:
% 1.  `detect_stim_times.m`: Function to find stimulation artifact times from raw data.
% 2.  `calculate_psth_metrics.m`: Function to compute PSTH and related metrics.
%     (These functions must be in the MATLAB path).
%
% Last Modified:
% YA 
% 11Aug2025
%
% -------------------------------------------------------------------------

clear; clc;

%% --- User Parameters ---
% --- List of file pairs to analyze ---
% Add as many structs to this cell array as needed.
% Each struct must contain a 'spikeFile' and a 'rawFile'.
filePairs = { ...
    struct('spikeFile', 'OWT220207_1I_DIV63_HUB63_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_1I_DIV63_HUB63_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2B_DIV63_HUB24_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2B_DIV63_HUB24_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2D_DIV63_HUB73_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2D_DIV63_HUB73_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2E_DIV63_HUB36_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2E_DIV63_HUB36_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2I_DIV63_HUB24_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2I_DIV63_HUB24_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_1I_DIV63_PER37_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_1I_DIV63_PER37_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2B_DIV63_PER52_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2B_DIV63_PER52_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2D_DIV63_PER58_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2D_DIV63_PER58_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2E_DIV63_PER71_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2E_DIV63_PER71_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2I_DIV63_PER83_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2I_DIV63_PER_6UA.mat') ...
};

% --- Analysis Parameters ---
fs = 25000;                     % Sampling frequency of the recording in Hz.
spikeMethod = 'bior1p5';        % Name of the spike detection method used. This corresponds to a field in the loaded spike data struct.
numChannels = 60;               % Total number of channels in the recording array.
artifact_window_ms = [0, 2];    % Time window [start, end] in ms after a stimulus to treat as an artifact. Spikes in this window are ignored.
psth_window_s = [0, 0.02];      % Analysis window [start, end] in seconds for the post-stimulus time histogram (PSTH). E.g., 0 to 20ms post-stimulus.
psth_bin_width_s = 0.001;       % Width of each bin for the reference histogram in seconds (e.g., 1 ms).
num_baseline_psths = 100;       % Number of baseline PSTHs to construct backwards in time from each stimulus onset.
baseline_duration_s = 0.02;     % Duration of each individual baseline window in seconds (e.g., 20ms). This should match psth_window_s

% --- Stimulation Detection Parameters (for "longblank" method) ---
min_blanking_duration_ms = 4;   % Minimum duration of a signal blanking (voltage flatline) to be considered a stimulation onset.
min_interval_ms = 2500;         % Minimum time interval between two consecutive detected stimuli to be considered valid (like a refractory period).


% --- Channel Remapping ---
indices = [24 26 29 32 35 37, 21 22 25 30 31 36 39 40, 19 20 23 28 33 38 41 42, 16 17 18 27 34 43 44 45, 15 14 13 4 57 48 47 46, 12 11 8 3 58 53 50 49, 10 9 6 1 60 55 52 51, 7 5 2 59 56 54];
ids = [21 31 41 51 61 71, 12 22 32 42 52 62 72 82, 13 23 33 43 53 63 73 83, 14 24 34 44 54 64 74 84, 15 25 35 45 55 65 75 85, 16 26 36 46 56 66 76 86, 17 27 37 47 57 67 77 87, 28 38 48 58 68 78];
channelMap = containers.Map('KeyType','double','ValueType','double');
for i = 1:numel(indices), channelMap(indices(i)) = ids(i); end


%% --- Main Analysis ---
for k = 1:numel(filePairs)
    % --- Get current file pair ---
    spikeFile = filePairs{k}.spikeFile;
    rawFile = filePairs{k}.rawFile;

    %% Display progress information in the command window.
    %fprintf('\n\n============================================================\n');
    %fprintf('STARTING ANALYSIS FOR FILE PAIR %d of %d:\n', k, numel(filePairs));
    %fprintf('  Spike File: %s\n', spikeFile);
    %fprintf('  Raw File: %s\n', rawFile);
    %fprintf('============================================================\n');

    % --- Generate a unique output directory for this file pair ---
    [~, baseName, ~] = fileparts(spikeFile);
    timestamp = datestr(now, 'ddmmmyyyy_HH:MM');
    outputDir = sprintf('PSTH_Analysis_%s_%s', baseName, timestamp);
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    %fprintf('Saving analysis plots to folder: %s\n', outputDir);

    % --- Load Spike Times and Determine Stim Times ---
    %fprintf('Loading and processing data...\n');
    S = load(spikeFile);
    % Check for spike data in two possible formats for compatibility. This is particularly relevant when processing Mona's spike detection files.
    if isfield(S, 'spikeTimes'), spikeTimesConverted = S.spikeTimes;
    elseif isfield(S, 'spikes'), fprintf('Converting ''spikes'' matrix...\n');
        [row, col] = find(S.spikes); spikeTimesConverted = cell(1, numChannels);
        for ch = 1:numChannels, spike_samples = row(col == ch); spike_sec = spike_samples / fs; spikeTimesConverted{ch} = struct(spikeMethod, spike_sec); end
    else, error('Spike data not found in file: %s', spikeFile); end
    
    %fprintf('Detecting stimulation events using "longblank" method...\n');

    [stimTimes_ms, ~, ~] = detect_stim_times(rawFile, numChannels, fs, ...
        [], min_blanking_duration_ms, [], min_interval_ms, []);
    
    % Convert stim times to seconds 
    stimTimes = stimTimes_ms / 1000;
    stimTimes = sort(unique(stimTimes(:))); % Ensure all stimulus times are unique and sorted chronologically.
    
    %fprintf('Found %d stimulation events.\n', length(stimTimes)); 

    % --- Network Response Analysis Loop ---
    networkResponse = []; % Initialize an empty array to store results for each channel in the current file pair.
    valid_channel_count = 0; % Reset counter for each file pair
 
    for file_idx = 1:numChannels
        if ~isKey(channelMap, file_idx), continue; end
        channel_id = channelMap(file_idx);

        % Extract spike times for the current channel, handling cases where data might be missing or empty.
        if file_idx > numChannels || file_idx < 1 || isempty(spikeTimesConverted{file_idx}) || ~isfield(spikeTimesConverted{file_idx}, spikeMethod), all_spike_times_s = [];
        else, all_spike_times_s = spikeTimesConverted{file_idx}.(spikeMethod); end
        if isempty(all_spike_times_s), continue; end

        % Remove spikes that fall within the stimulation artifact window for each stimulus.
        spikeTimes_cleaned_s = all_spike_times_s;
        for stimIdx = 1:numel(stimTimes), stimTime = stimTimes(stimIdx);
            spikeTimes_cleaned_s = spikeTimes_cleaned_s(spikeTimes_cleaned_s < (stimTime + artifact_window_ms(1)/1000) | spikeTimes_cleaned_s >= (stimTime + artifact_window_ms(2)/1000));
        end
        spikeTimes_cleaned_s = sort(spikeTimes_cleaned_s(:));

        % --- PSTH Calculation for Response Window ---
        [response, resp_metrics] = calculate_psth_metrics(spikeTimes_cleaned_s, stimTimes, psth_window_s, psth_bin_width_s);
        if isempty(response.psth_samples), continue; end % If no spikes in the window, skip.
        
        % --- PSTH Calculation using Multiple Baselines ---
        %fprintf('Calculating %d baseline PSTHs...\n', num_baseline_psths);
        baseline_aucs = zeros(num_baseline_psths, 1);
        all_baseline_psth_smooth = [];
        for i = 1:num_baseline_psths
            % Define the time window for the current baseline period, moving backwards from the stimulus.
            start_s = -(i * baseline_duration_s);
            end_s = -((i-1) * baseline_duration_s);
            current_baseline_window_s = [start_s, end_s];

            % Apply artifact blanking within the baseline window as well.
            % This removes spikes that fall in the range [start_s + artifact_window_ms(1), start_s + artifact_window_ms(2))
            % relative to each stimulus time, so that when plotted as 0–20 ms the 0–artifact_window_ms region is blanked.
            spikeTimes_cleaned_baseline_s = spikeTimes_cleaned_s;
            blank_start_offset_s = start_s + artifact_window_ms(1)/1000;
            blank_end_offset_s   = start_s + artifact_window_ms(2)/1000;
            for stimIdx = 1:numel(stimTimes)
                stimTime = stimTimes(stimIdx);
                spikeTimes_cleaned_baseline_s = spikeTimes_cleaned_baseline_s( ...
                    spikeTimes_cleaned_baseline_s < (stimTime + blank_start_offset_s) | ...
                    spikeTimes_cleaned_baseline_s >= (stimTime + blank_end_offset_s) );
            end

            [~, base_metrics] = calculate_psth_metrics(spikeTimes_cleaned_baseline_s, stimTimes, current_baseline_window_s, psth_bin_width_s);
            baseline_aucs(i) = base_metrics.auc;
            if isempty(all_baseline_psth_smooth), all_baseline_psth_smooth = zeros(num_baseline_psths, length(base_metrics.psth_smooth)); end
            all_baseline_psth_smooth(i, :) = base_metrics.psth_smooth;
        end
        
        % --- Calculate Baseline-corrected AUC ---
        mean_baseline_auc = mean(baseline_aucs);
        auc_corrected = resp_metrics.auc - mean_baseline_auc; % Subtract the mean baseline AUC from the response AUC. TODO: also make a percentage value
        mean_baseline_psth = mean(all_baseline_psth_smooth, 1); % Calculate the average smoothed baseline PSTH.
        %fprintf('Response AUC: %.4f, Mean Baseline AUC: %.4f, Corrected AUC: %.4f\n', resp_metrics.auc, mean_baseline_auc, auc_corrected);

        % --- Decay Analysis using Half-Max from Peak ---
        [Rmax, Rmax_idx] = max(resp_metrics.psth_smooth);
        halfRmax = Rmax / 2;
        % Find the first point after the peak where the firing rate drops to or below half-max.
        halfRmax_idx = find(resp_metrics.psth_smooth(Rmax_idx:end) <= halfRmax, 1, 'first');
        if ~isempty(halfRmax_idx)
            halfRmax_idx = halfRmax_idx + Rmax_idx - 1;
            halfRmax_time_s = resp_metrics.time_vector_s(halfRmax_idx);
            halfRmax_val = resp_metrics.psth_smooth(halfRmax_idx);
        else % If no such point was found
            halfRmax_time_s = NaN;
            halfRmax_val = NaN;
        end

        % --- Plot Generation ---
        fig = figure('Position',[100 100 1200 900], 'Visible', 'off');
        psth_window_ms = psth_window_s * 1000; % Convert PSTH window to milliseconds for plotting.
        sgtitle(sprintf('Channel %d | Peak Rate: %.1f Hz | Corrected AUC: %.3f', channel_id, resp_metrics.peak_firing_rate, auc_corrected), 'FontWeight', 'bold');

        % Subplot 1: Spike Raster Plot
        ax1 = subplot(2,2,1); hold on;
        for trial_idx = 1:length(response.spikeTimes_byEvent), trial_spikes_s = response.spikeTimes_byEvent{trial_idx}; if ~isempty(trial_spikes_s), plot(trial_spikes_s * 1000, trial_idx * ones(size(trial_spikes_s)), 'r.', 'MarkerSize', 5); end; end
        hold off; set(gca, 'YDir', 'reverse'); xlim(psth_window_ms); ylim([0 length(stimTimes)+1]); ylabel('Trial Number'); title('Spike Raster (Response)'); grid on;

        % Subplot 2: Diagnostic plot comparing response PSTH to all baseline PSTHs.
        ax2 = subplot(2,2,2); hold on;
        baseline_time_ms = (base_metrics.time_vector_s - current_baseline_window_s(1)) * 1000;
        for i = 1:num_baseline_psths, plot(baseline_time_ms, all_baseline_psth_smooth(i, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); end
        p1_diag = plot(resp_metrics.time_vector_s*1000, resp_metrics.psth_smooth, 'r-', 'LineWidth', 2);
        p2_diag = plot(baseline_time_ms, mean_baseline_psth, 'k-', 'LineWidth', 2);
        hold off; title('Diagnostic: Response vs. Baselines'); ylabel('Firing Rate (spikes/s)'); xlabel('Time from stimulus (ms)'); legend([p1_diag, p2_diag], 'Response', 'Mean Baseline', 'Location', 'Best'); grid on;

        % Subplot 3: Smoothed response PSTH with metrics.
        ax3 = subplot(2,2,3); hold on;
        edges_s = psth_window_s(1):psth_bin_width_s:psth_window_s(2);
        bar(edges_s * 1000, response.psth_histogram, 1, 'FaceColor',.7*[1 1 1],'EdgeColor',.8*[1 1 1], 'HandleVisibility','off');
        plot(resp_metrics.time_vector_s * 1000, resp_metrics.psth_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed PSTH');
        plot(resp_metrics.peak_time_s*1000, resp_metrics.peak_firing_rate, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'R_{max}');
        text(resp_metrics.peak_time_s*1000, resp_metrics.peak_firing_rate, sprintf(' R_{max}: %.1f Hz', resp_metrics.peak_firing_rate), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');
        if ~isnan(halfRmax_time_s)
            plot(halfRmax_time_s*1000, halfRmax_val, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'Half R_{max}');
            text(halfRmax_time_s*1000, halfRmax_val, sprintf(' Half R_{max} @ %.1f ms', halfRmax_time_s*1000), 'VerticalAlignment', 'top', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');
        end
        hold off; ylabel('Firing Rate (spikes/s)'); xlabel('Time from stimulus (ms)'); title('Smoothed Response PSTH & Metrics'); legend('Location', 'northeast'); grid on;

        % Subplot 4: Adaptive kernel bandwidth used for smoothing. TODO: remove this plot
        ax4 = subplot(2,2,4);
        plot(resp_metrics.time_vector_s * 1000, resp_metrics.kernel_bandwidth_s * 1000, 'b-', 'LineWidth', 2);
        ylabel('Bandwidth (ms)'); xlabel('Time from stimulus (ms)'); title('Adaptive Kernel Bandwidth (Response)'); grid on;
        
        linkaxes([ax1, ax2, ax3, ax4], 'x'); xlim(psth_window_ms);
        plot_filename = fullfile(outputDir, sprintf('Analysis_Plot_Chan_%d.png', channel_id));
        saveas(fig, plot_filename); fprintf('Saved consolidated plot to %s\n', plot_filename); close(fig);

      % --- Save Results ---
        valid_channel_count = valid_channel_count + 1;
        % Store all calculated metrics for the current channel in the 'networkResponse' struct array.
        networkResponse(valid_channel_count).channel_id = channel_id;
        networkResponse(valid_channel_count).file_index = file_idx;
        networkResponse(valid_channel_count).auc_response = resp_metrics.auc;
        networkResponse(valid_channel_count).auc_baseline_mean = mean_baseline_auc;
        networkResponse(valid_channel_count).auc_corrected = auc_corrected;
        networkResponse(valid_channel_count).peak_firing_rate_hz = resp_metrics.peak_firing_rate;
        networkResponse(valid_channel_count).peak_time_ms = resp_metrics.peak_time_s * 1000;
        networkResponse(valid_channel_count).halfRmax_time_ms = halfRmax_time_s * 1000;
    end

    % --- Save to Output Directory ---
    if ~isempty(networkResponse)
        output_filename = fullfile(outputDir, sprintf('networkResponse_all_channels_%s.mat', timestamp));
        save(output_filename, 'networkResponse'); end
    
end

