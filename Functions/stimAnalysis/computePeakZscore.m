function peakZscore = computePeakZscore(all_spike_times_s, allStimTimes, psth_window_s, baseline_window_s, psth_gaussian_width_ms)
% COMPUTEPEAKZSCORE Computes the peak z-score from a smoothed PSTH.
%
% For a single electrode, this function:
% 1. Calculates the smoothed PSTH using calculate_psth_metrics.
% 2. Calculates baseline statistics from trial-by-trial firing rates in the pre-stimulus window.
% 3. Converts the smoothed PSTH to z-scores using the baseline mean and std.
% 4. Returns the peak (maximum) of the z-scored PSTH.
%
% This function's methodology matches that of stimActivityAnalysis.m.
%
% INPUTS
% ------
% all_spike_times_s : double vector
%     Vector of spike times in seconds for one electrode.
% allStimTimes : double vector
%     All stimulation event times in seconds to align to.
% psth_window_s : [1 x 2] double
%     PSTH analysis window [start, end] in seconds relative to stim times.
% baseline_window_s : [1 x 2] double
%     Window for baseline firing rate calculation [start, end] in seconds.
% psth_gaussian_width_ms : double
%     Width of the Gaussian smoothing kernel in milliseconds.
%
% OUTPUTS
% -------
% peakZscore : double
%     The peak z-score value in the post-stimulus period. Returns 0 if
%     insufficient data or an error occurs.

    if isempty(all_spike_times_s)
        peakZscore = 0;  % No spikes
        return;
    end

    numStimEvents = length(allStimTimes);

    % Calculate smoothed PSTH (matches stimActivityAnalysis.m)
    try
        psth_bin_width_s = 0.001; % Corresponds to default in calculate_psth_metrics
        [~, resp_metrics] = calculate_psth_metrics(...
            all_spike_times_s, allStimTimes, psth_window_s, psth_bin_width_s, ...
            'smoothing_method', 'gaussian', 'gaussian_width_ms', psth_gaussian_width_ms, ...
            'auc_start_s', 0);
    catch
        peakZscore = 0;  % Error in PSTH calculation
        return;
    end

    % Calculate baseline statistics from trial-by-trial firing rates (matches stimActivityAnalysis.m)
    baseline_firing_rates = zeros(1, numStimEvents);
    effective_window_duration = psth_window_s(2) - psth_window_s(1);  % Total analysis window (matches stimActivityAnalysis.m line 483)
    
    for stimIdx = 1:numStimEvents
        stimTime = allStimTimes(stimIdx);
        
        % Baseline period firing rate (pre-stimulus)
        baseline_start = stimTime + baseline_window_s(1);
        baseline_end = stimTime + baseline_window_s(2);
        
        baseline_spikes = all_spike_times_s(...
            all_spike_times_s >= baseline_start & all_spike_times_s < baseline_end);
        baseline_firing_rates(stimIdx) = length(baseline_spikes) / effective_window_duration;
    end
    
    % Calculate baseline statistics
    baseline_mean_hz = mean(baseline_firing_rates);
    baseline_std_hz = std(baseline_firing_rates);
    
    % Handle case of zero standard deviation
    if baseline_std_hz == 0
        baseline_std_hz = eps; % Avoid division by zero
    end
    
    % Convert smoothed PSTH to z-scores (matches stimActivityAnalysis.m)
    zscore_psth = (resp_metrics.psth_smooth - baseline_mean_hz) ./ baseline_std_hz;
    
    % Find peak z-score in the post-stimulus period
    post_stim_mask = resp_metrics.time_vector_s >= 0;
    if any(post_stim_mask)
        peakZscore = max(zscore_psth(post_stim_mask));
    else
        peakZscore = 0; % No post-stimulus bins
    end
end
