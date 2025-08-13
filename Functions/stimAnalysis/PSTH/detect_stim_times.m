function [stimTimes, stimElectrodes, stimtimes] = detect_stim_times(rawFile, numChannels, fs, ~, min_blanking_duration_ms, ~, min_interval_ms, ~)
%DETECT_STIM_TIMES Detects stimulation times from raw data file using flat/blanking period detection.
%
%   This version uses a "longblank" detection method. It identifies periods where
%   the signal is flat for a minimum duration, considering these as stimulation
%   blanking artifacts.
%
%   [stimTimes, stimElectrodes, stimtimes] = detect_stim_times(rawFile, numChannels, fs, min_blanking_duration_ms, min_interval_ms)
%
% Inputs:
%   rawFile                 - Path to the .mat file containing raw data (variable 'dat')
%   numChannels             - Number of channels
%   fs                      - Sampling rate (Hz)
%   ~ (stimThreshold)       - Not used in this method.
%   min_blanking_duration_ms- Minimum duration of a flat signal to be considered a stim blanking event (ms).
%   ~ (flat_thresh)         - Not used in this method.
%   min_interval_ms         - Minimum refractory period between consecutive stims (ms).
%   ~ (flat_search_window_ms) - Not used in this method.
%
% Outputs:
%   stimTimes       - Array of stimulation times (ms)
%   stimElectrodes  - Array of electrode numbers for each stimulation
%   stimtimes       - Table with columns 'Time_ms' and 'Electrode'

    % Load raw data
    R = load(rawFile);
    if isfield(R, 'dat')
        dat = double(R.dat);
    else
        error('Raw data variable "dat" not found in %s', rawFile);
    end

    stim_detection_labels = {};

    % Convert parameters from ms to samples
    min_duration_samples = round(min_blanking_duration_ms * fs / 1000);
    min_interval_samples = round(min_interval_ms * fs / 1000);

    % Loop through each channel to detect stimulations
    for channel_idx = 1:numChannels
        channelDat = dat(:, channel_idx);

        % Find where the signal value changes
        change_points = [true; diff(channelDat) ~= 0];

        % Assign a group ID to each run of constant value
        group_id = cumsum(change_points);

        % Count the length of each run
        counts = accumarray(group_id, 1);

        % Find which groups represent long flat periods
        valid_groups = find(counts >= min_duration_samples);

        % Get the start index of each run
        [~, start_indices_all] = unique(group_id, 'first');

        % Filter for start indices of valid (long enough) flat periods
        start_indices = start_indices_all(ismember(group_id(start_indices_all), valid_groups));
        
        if isempty(start_indices), continue, end

        % Apply refractory period to remove closely spaced detections
        if length(start_indices) > 1
            keep = [true; diff(start_indices) > min_interval_samples];
            start_indices = start_indices(keep);
        end

        % Store detected stimulation times and the channel
        for i = 1:length(start_indices)
            stim_idx = start_indices(i);
            stim_time_ms = (stim_idx / fs) * 1000;
            stim_detection_labels = [stim_detection_labels; {stim_time_ms, channel_idx}]; %#ok<AGROW>
        end
    end

    % Convert cell array to a table and sort by time
    if ~isempty(stim_detection_labels)
        stimtimes = cell2table(stim_detection_labels, 'VariableNames', {'Time_ms', 'Electrode'});
        stimtimes = sortrows(stimtimes, 'Time_ms');
    else
        stimtimes = cell2table(cell(0,2), 'VariableNames', {'Time_ms', 'Electrode'});
    end

    % Final outputs
    stimTimes = stimtimes.Time_ms;
    stimElectrodes = stimtimes.Electrode;
end