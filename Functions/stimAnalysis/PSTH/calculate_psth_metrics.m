function [psth_data, metrics] = calculate_psth_metrics(all_spike_times_s, stim_times_s, window_s, bin_width_s)
% CALCULATE_PSTH_METRICS Extracts spikes in a given window relative to
% stimulus times, computes a histogram, generates a smoothed PSTH using
% ssvkernel, and calculates key metrics.
%
% INPUTS:
%   all_spike_times_s - Vector of all spike times in seconds.
%   stim_times_s      - Vector of stimulus event times in seconds.
%   window_s          - A 1x2 vector defining the time window relative to
%                       each stimulus [start, end] in seconds.
%   bin_width_s       - The bin width for the raw PSTH histogram in seconds.
%
% OUTPUTS:
%   psth_data         - A struct containing intermediate data:
%     .spikeTimes_byEvent - Cell array of spike times per trial.
%     .psth_samples       - Pooled vector of all spike times relative to stim.
%     .psth_histogram     - The raw histogram of firing rates.
%
%   metrics           - A struct containing the final calculated metrics:
%     .time_vector_s      - Time vector for the smoothed PSTH.
%     .psth_smooth        - The smoothed PSTH (firing rate in spikes/s).
%     .kernel_bandwidth_s - The adaptive kernel bandwidth from ssvkernel.
%     .auc                - Area under the curve of the smoothed PSTH.
%     .peak_firing_rate   - Peak firing rate of the smoothed PSTH.
%     .peak_time_s        - Time of the peak firing rate.

    num_trials = length(stim_times_s);
    MIN_SPIKES_FOR_KDE = 5; % Set a threshold for minimum number of spikes to run ssvkernel

    % Extract spikes within the specified window for each trial
    out = WithinRanges(all_spike_times_s, stim_times_s + window_s, (1:num_trials)', 'matrix');
    spikeTimes_byEvent = arrayfun(@(n) all_spike_times_s(logical(out(:,n))) - stim_times_s(n), 1:num_trials, 'uni', 0)';
    psth_samples = cell2mat(spikeTimes_byEvent);

    % Create raw histogram
    edges_s = window_s(1):bin_width_s:window_s(2);
    if ~isempty(psth_samples)
        b = histc(psth_samples, edges_s);
        psth_histogram = b / (num_trials * bin_width_s);
    else
        psth_histogram = zeros(size(edges_s));
    end

    % Prepare outputs
    psth_data.spikeTimes_byEvent = spikeTimes_byEvent;
    psth_data.psth_samples = psth_samples;
    psth_data.psth_histogram = psth_histogram;

    % Smooth PSTH with ssvkernel and calculate metrics
    if numel(psth_samples) >= MIN_SPIKES_FOR_KDE
        L = 1000; % Number of points for smoothing
        t_s = linspace(window_s(1), window_s(2), L);
        [yv_pdf, tv_s, optw_variable_s] = ssvkernel(psth_samples, t_s);
        
        avg_spikes_per_trial = length(psth_samples) / num_trials;
        yv = yv_pdf * avg_spikes_per_trial; % Convert density to rate
        
        metrics.time_vector_s = tv_s;
        metrics.psth_smooth = yv;
        metrics.kernel_bandwidth_s = optw_variable_s;
        metrics.auc = trapz(tv_s, yv);
        [metrics.peak_firing_rate, max_idx] = max(yv);
        metrics.peak_time_s = tv_s(max_idx);
    else
        % Handle case with too few spikes for reliable KDE
        L = 1000;
        tv_s = linspace(window_s(1), window_s(2), L);
        metrics.time_vector_s = tv_s;
        metrics.psth_smooth = zeros(1, L);
        metrics.kernel_bandwidth_s = zeros(1, L);
        metrics.auc = 0;
        metrics.peak_firing_rate = 0;
        metrics.peak_time_s = NaN;
    end
end