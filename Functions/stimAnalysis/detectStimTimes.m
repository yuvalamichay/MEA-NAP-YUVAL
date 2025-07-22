function stimInfo = detectStimTimes(rawData, Params, channelNames, coords)
% Detects electrical sitmulation time from filtered voltage trace
% Parameters
% ----------
% filteredData : matrix
% matrix of shape (numTimeSamples, numChannels)
% Params : struct 
% structure with the following required fields 
%      stimDetectionMethod : str 
%           method of detecting electrical stimulation
% channelNames : vector 
% vector (or cell) of channel names, usually integer values
% Output
% ------
% stimInfo : cell 
% cell where each entry is a structure 
% is structure has the following fields
%      elecStimTimes : vector of electrical stimulation onset times (seconds)
%      elecStimDur : vector of electrical stimulation duration (seconds)
%      paddedChannels: vector of electrode indices where any stimulation time was 'padded' due to subthreshold stimulation artifacts

% Set parameters 
stimDetectionMethod = Params.stimDetectionMethod;
stimRefPeriod = Params.stimRefractoryPeriod; 
stimDur = Params.stimDuration;


numChannels = size(rawData, 2);
stimInfo = cell(numChannels, 1);
paddedChannels = [];

if strcmp(stimDetectionMethod, 'blanking')

    % IDEA: take the mode of the blanking times per electrode, and use
    % those as the stimulation times
    % NOTE: This is quite slow at the moment, not sure why...

    num_blanks_per_electrode = zeros(1, numChannels);
    electrode_blank_times = {};

    for channel_idx = 1:numChannels
        channelDat = rawData(:, channel_idx);
        
        % Minimum duration of a run (currently in samples)
        min_duration = round(Params.minBlankingDuration * Params.fs); % 25 works
        
        % Find where values change
        change_points = [1; diff(channelDat) ~= 0]; % [true (diff(channelDat) ~= 0)];
        
        % Assign group IDs to each run of constant value
        group_id = cumsum(change_points);
        
        % Count length of each run
        counts = accumarray(group_id(:), 1);
        
        % Find which groups meet the minimum duration
        valid_groups = find(counts >= min_duration);
        
        % Find start indices of each group
        [~, start_indices] = unique(group_id, 'first');
        
        % Select only those with valid length
        start_indices = start_indices(ismember(group_id(start_indices), valid_groups));

        num_blanks_per_electrode(channel_idx) = length(start_indices);
        electrode_blank_times{channel_idx} = start_indices / Params.fs;

    end
    
    blanks_count_mode = mode(num_blanks_per_electrode);
    electrodes_with_mode_count = find(num_blanks_per_electrode == blanks_count_mode);
    allStimTimesTemplate = median(horzcat(electrode_blank_times{electrodes_with_mode_count}), 2);

    if isempty(allStimTimesTemplate)
        fprintf('WARNING: NO BLANKING DETECTED')
    end

end


if strcmp(stimDetectionMethod, 'blanking')
    %{
    if strcmp(Params.stimProcessingMethod, 'filter')
        lowpass = Params.filterLowPass;
        highpass = Params.filterHighPass;
        wn = [lowpass highpass] / (Params.fs / 2);
        filterOrder = 3;
        [b, a] = butter(filterOrder, wn);
        % WIP: 
        filteredData = zeros(size(rawData));
        for channelIdx = 1:size(stimRawData.dat, 2)
            data = stimRawData.dat(:, channelIdx);
            trace = filtfilt(b, a, double(data));
            filteredData(:, channelIdx) = trace;
         end 
    else
    %}
        medianAbsDeviation = median(abs(rawData - mean(rawData, 1)), 1);
        medianZscore = abs(rawData - median(rawData, 1)) ./ medianAbsDeviation;
    % end 
end


% Do electrical stimulation detection (and assign stimulation duration)
for channel_idx = 1:numChannels
    traceMean = mean(rawData(:, channel_idx));
    traceStd = std(rawData(:, channel_idx));
    

    if strcmp(stimDetectionMethod, 'absPosThreshold')
        stimThreshold = Params.stimDetectionVal;
        elecStimTimes = find(rawData(:, channel_idx) > stimThreshold) / Params.fs;
    elseif strcmp(stimDetectionMethod, 'absNegThreshold')
        stimThreshold = Params.stimDetectionVal;
        elecStimTimes = find(rawData(:, channel_idx) < stimThreshold) / Params.fs;
    elseif strcmp(stimDetectionMethod, 'stdNeg')
        stimThreshold = traceMean - traceStd * Params.stimDetectionVal;
        elecStimTimes = find(rawData(:, channel_idx) < stimThreshold) / Params.fs;
    elseif strcmp(stimDetectionMethod, 'blanking')
        % First take the approximate stimulation time from the trace 
        % Then find the closest stimulation times according to the blanking
        % times
        stimThreshold = Params.stimDetectionVal;

        % TODO: this is not as good as the filtering approach... need to
        % have very long refractory period...
        approxElecStimTimes = find(medianZscore(:, channel_idx) > stimThreshold) / Params.fs;

        keepIdx = true(size(approxElecStimTimes)); % Logical mask for keeping elements
        lastValidIdx = 1; % Track last valid stim time
        
        for stimIdx = 2:length(approxElecStimTimes)
            if approxElecStimTimes(stimIdx) <= approxElecStimTimes(lastValidIdx) + stimRefPeriod
                keepIdx(stimIdx) = false; % Mark for removal
            else
                lastValidIdx = stimIdx; % Update last valid index
            end
        end
        approxElecStimTimes = approxElecStimTimes(keepIdx); % Keep only valid elements
        
        elecStimTimes = [];
        for stimIdx = 1:length(approxElecStimTimes)
            % NOTE: Here we assume that the blanking always occur before
            % the stimulation artifact, this resolves cases where two
            % blanking periods are equally close the the stimulation
            % artifact
            allStimTimesTemplateBeforeStim = allStimTimesTemplate(allStimTimesTemplate < approxElecStimTimes(stimIdx));
            [~, minIdx] = min(abs(approxElecStimTimes(stimIdx) - allStimTimesTemplateBeforeStim));
            if ~isempty(minIdx)
                elecStimTimes(stimIdx) = allStimTimesTemplateBeforeStim(minIdx);
            end 
        end

    else 
        error('No valid stimulus detection specified')
    end 
    
    % Remove stim within refractory period of each other 
    % V1 : Slow 
    %{
    for stimIdx = 1:length(elecStimTimes)
        
        stimTime = elecStimTimes(stimIdx);

        if ~isnan(stimTime)
            removeIndex = find( ...
                 (elecStimTimes > stimTime) & ...
                 (elecStimTimes <= stimTime + stimRefPeriod) ...
                );
            elecStimTimes(removeIndex) = nan;
        end

    end
    elecStimTimes = elecStimTimes(~isnan(elecStimTimes));
    %}

    % V2: Faster
    %
    keepIdx = true(size(elecStimTimes)); % Logical mask for keeping elements
    lastValidIdx = 1; % Track last valid stim time
    
    for stimIdx = 2:length(elecStimTimes)
        if elecStimTimes(stimIdx) <= elecStimTimes(lastValidIdx) + stimRefPeriod
            keepIdx(stimIdx) = false; % Mark for removal
        else
            lastValidIdx = stimIdx; % Update last valid index
        end
    end
    elecStimTimes = elecStimTimes(keepIdx); % Keep only valid elements
    %}

  %% PADDING 
    % --------------------------------------------------
    % This section identifies and 'pads' missing stimulation events
    % by predicting when they should have occurred, based on expected inter-stimulation intervals.
    % This helps to correct cases where some stimulation events are not detected
    % due to subthreshold stimulation artifacts.
    % Padding occurs only if:
    %   - stimDetectionMethod is 'blanking'
    %   - Params.padStimArtifacts is true
    % --------------------------------------------------
 
  if strcmp(stimDetectionMethod, 'blanking')
        originalStimTimes = elecStimTimes; % Saving elecStimTimes identified thus far for later comparison
        
  if length(elecStimTimes) > 5 % minimum number of stim times detected per channel for padding     
        stim_intervals = diff(elecStimTimes);  

        if isfield(Params, 'knownStimFrequency') && ~isempty(Params.knownStimFrequency)
            if numel(Params.knownStimFrequency) == 1 % Only one electrode is stimulated
                expected_interval = 1/Params.knownStimFrequency; % Global value
            else
                expected_interval = 1/Params.knownStimFrequency(channel_idx); % Individual value per electrode
            end
        else
            [counts, edges] = histcounts(stim_intervals, 'BinWidth', 0.001); % TODO: confirm this bin width?
            [~, max_count_idx] = max(counts); % Identifying most common interval (mode)
            expected_interval = (edges(max_count_idx) + edges(max_count_idx+1))/2; % Center of the bin with the highest count is taken as the expected interval
        end

        start_idx = round(min(elecStimTimes) * Params.fs); % Start of window is index of first stimulation time
        end_idx = size(rawData,1); % End of window is final index of recording

        expected_interval_samples = round(expected_interval * Params.fs);
        expected_stim_indices = start_idx:expected_interval_samples:end_idx;
        expected_stim_times = expected_stim_indices / Params.fs;

        half_interval_samples = round(expected_interval_samples / 2); % Constructing a window that is half the expected interval to pad only missing stimulation times
        originalStimIndices = round(originalStimTimes * Params.fs);

        unmatched_expected = [];
        for i = 1:length(expected_stim_indices)
            exp_idx = expected_stim_indices(i);
            in_any_window = any(abs(exp_idx - originalStimIndices) <= half_interval_samples);
            if ~in_any_window
                unmatched_expected(end+1) = exp_idx;
            end
        end

        padded_stim_times = unmatched_expected / Params.fs;
        allStimTimes = unique(sort([originalStimTimes padded_stim_times]));
        isOriginal = ismember(allStimTimes, originalStimTimes);

     % Annotate each stim event: type 'original' or 'padded'
        stimTimeInfo = repmat(struct('time', [], 'type', '', 'dur', []), length(allStimTimes), 1);
        for ii = 1:length(allStimTimes)
            stimTimeInfo(ii).time = allStimTimes(ii);
            stimTimeInfo(ii).type = isOriginal(ii) * "original" + ~isOriginal(ii) * "padded";
            stimTimeInfo(ii).dur = stimDur;
        end

        elecStimTimes = [stimTimeInfo.time];
        elecStimDur = [stimTimeInfo.dur];

        % Track padded channels
        if any(~isOriginal)
            paddedChannels(end+1) = channel_idx;
        end

        fprintf('Channel %d: Added %d unmatched stim times (total: %d, originally detected: %d)\n', ...
            channel_idx, length(elecStimTimes) - length(originalStimTimes), ...
            length(elecStimTimes), length(originalStimTimes));
    else
        % No padding, just annotate detected events
        allStimTimes = elecStimTimes;
        stimTimeInfo = repmat(struct('time', [], 'type', '', 'dur', []), length(allStimTimes), 1);
        for ii = 1:length(allStimTimes)
            stimTimeInfo(ii).time = allStimTimes(ii);
            stimTimeInfo(ii).type = "original";
            stimTimeInfo(ii).dur = stimDur;
        end
        elecStimDur = [stimTimeInfo.dur];

        fprintf('Channel %d: No padding needed (events detected: %d)\n', channel_idx, length(elecStimTimes));
    end

    % Assign output structure for this channel
    stimStruct = struct();
    stimStruct.elecStimTimes = [stimTimeInfo.time];
    stimStruct.elecStimDur = [stimTimeInfo.dur];
    stimStruct.stimTimeInfo = stimTimeInfo;
    stimStruct.channelName = channelNames(channel_idx);
    stimStruct.coords = coords(channel_idx, :);

    if exist('allStimTimesTemplate', 'var')
        stimStruct.allStimTimesTemplate = allStimTimesTemplate;
    end
    if isfield(Params, 'padStimArtifacts') && Params.padStimArtifacts
        stimStruct.originalStimTimes = originalStimTimes;
    end

    stimInfo{channel_idx} = stimStruct;
else
    % Not blanking or not padding: preserve existing structure
    stimTimeInfo = repmat(struct('time', [], 'type', '', 'dur', []), length(elecStimTimes), 1);
    for ii = 1:length(elecStimTimes)
        stimTimeInfo(ii).time = elecStimTimes(ii);
        stimTimeInfo(ii).type = "original";
        stimTimeInfo(ii).dur = stimDur;
    end

    stimStruct = struct();
    stimStruct.elecStimTimes = elecStimTimes;
    stimStruct.elecStimDur = repmat(stimDur, length(elecStimTimes), 1);
    stimStruct.stimTimeInfo = stimTimeInfo;
    stimStruct.channelName = channelNames(channel_idx);
    stimStruct.coords = coords(channel_idx, :);
    if strcmp(stimDetectionMethod, 'blanking') && exist('allStimTimesTemplate', 'var')
        stimStruct.allStimTimesTemplate = allStimTimesTemplate;
    end

    stimInfo{channel_idx} = stimStruct;
end
