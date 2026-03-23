function shuffleResults = stimShuffleTest(spikeData, allStimTimes, Params, Info)
% STIMSHUFFLETEST  Circular-shift shuffle test for post-stimulus peak z-score significance.
%
% Builds a null distribution of post-stimulus peak z-score values by circularly
% shifting each electrode's spike times by a random offset (with
% wrap-around within the recording duration), then recomputing the
% smoothed PSTH and its peak z-score relative to baseline.  Significance is
% assessed per electrode with a two-tailed test (p < 0.05) using the
% empirical 2.5th and 97.5th percentiles of the null distribution.
% No multiple-comparison correction is applied.
%
% INPUTS
% ------
% spikeData : struct
%     Must contain:
%       .spikeTimes  - cell array {1 x numChannels}, each entry is a struct
%                      with a field named Params.SpikesMethod containing a
%                      vector of spike times in seconds.
%       .stimInfo    - cell array {1 x numChannels} with stimulation info
%                      (used only for numChannels).
% allStimTimes : double vector
%     All stimulation event times in seconds to align to.
% Params : struct
%     Must contain:
%       .SpikesMethod          - string, field name for spike times
%       .stimAnalysisWindow    - [preStart, postEnd] in seconds (e.g. [-0.5 1])
%     Optional:
%       .Nshuffles             - number of circular-shift shuffles.
%                                Default: 500
%       .shuffleAlpha          - significance level for the two-tailed test.
%                                Default: 0.05
% Info : struct
%     Must contain:
%       .duration_s - recording duration in seconds (for wrap-around)
%
% OUTPUTS
% -------
% shuffleResults : struct with the following fields
%   .peakZscore_obs    - [numChannels x 1] observed peak z-score for each electrode
%   .peakZscore_null   - [numChannels x Nshuffles] null peak z-score distributions
%   .pctile_lo         - [numChannels x 1] lower percentile bound (2.5th)
%   .pctile_hi         - [numChannels x 1] upper percentile bound (97.5th)
%   .isSigLo           - [numChannels x 1] logical, true if peakZscore_obs < pctile_lo
%   .isSigHi           - [numChannels x 1] logical, true if peakZscore_obs > pctile_hi
%   .isSignificant     - [numChannels x 1] logical, true if significant (either tail)
%   .Nshuffles         - scalar, number of shuffles performed
%   .alpha             - scalar, significance level used
%   .postStimWindow    - [1 x 2] the post-stimulus window used for peak detection
%
% PROCEDURE
% ---------
% 1. Compute observed smoothed PSTH and peak z-score for each electrode.
% 2. For each shuffle iteration:
%      a. For each electrode, circularly shift its spike train by a random
%         offset uniformly drawn from (0, recordingDuration) with wrap-around.
%      b. Recompute the smoothed PSTH and peak z-score for every electrode
%         using the shifted spike times aligned to the *original* stim times.
% 3. For each electrode, determine the 2.5th and 97.5th percentiles of its
%    null peak z-score distribution.
% 4. Mark an electrode as significant if observed peak z-score falls outside the
%    null percentile interval.
%
% REFERENCE
% ---------
% Circular-shift approach inspired by adjM_thr_parallel.m (STTC thresholding).
% Peak z-score calculation matches stimActivityAnalysis.m methodology.
%
% Authors: GitHub Copilot & MEA-NAP team
% Date:    March 2026

%% Parse optional parameters
if ~isfield(Params, 'Nshuffles') || isempty(Params.Nshuffles)
    Nshuffles = 500;
else
    Nshuffles = Params.Nshuffles;
end

if ~isfield(Params, 'shuffleAlpha') || isempty(Params.shuffleAlpha)
    alpha = 0.05;
else
    alpha = Params.shuffleAlpha;
end

%% Setup
numChannels = length(spikeData.stimInfo);
duration_s = Info.duration_s;

% Analysis windows (matching stimActivityAnalysis.m)
psth_window_s = Params.stimAnalysisWindow;  % Full analysis window
poststim_duration_s = psth_window_s(2) - 0;  % Duration from stimulus to end of post-stim window
baseline_window_s = [-poststim_duration_s, 0];  % Baseline window for d-prime calculation (same duration as post-stim)

% PSTH parameters (matching stimActivityAnalysis.m)
psth_gaussian_width_ms = 1;  % Gaussian kernel width (matching stimActivityAnalysis.m line 416)

%% 1. Compute observed peak z-score
spikeMethod = Params.SpikesMethod;
peakZscore_obs = computePeakZscoreForEachChannel(spikeData.spikeTimes, allStimTimes, ...
    psth_window_s, baseline_window_s, psth_gaussian_width_ms, ...
    numChannels, spikeMethod);

%% 2. Build null distribution via circular shift
peakZscore_null = zeros(numChannels, Nshuffles);

% Check for Parallel Computing Toolbox
matlabInstallation = ver;
toolboxNames = {matlabInstallation.Name};
parallelToolboxInstalled = any(strcmp(toolboxNames, 'Parallel Computing Toolbox'));

if parallelToolboxInstalled
    parfor shuffleIdx = 1:Nshuffles
        % Create shifted copy of spike times
        shuffledSpikeTimes = spikeData.spikeTimes;  %#ok<PFBNS>
        for chIdx = 1:numChannels
            chSpikes = shuffledSpikeTimes{chIdx}.(spikeMethod);  %#ok<PFBNS>
            % Random circular shift: uniform in (0, duration_s)
            delta = rand() * duration_s;  %#ok<PFBNS>
            shiftedSpikes = chSpikes + delta;
            % Wrap around
            shiftedSpikes(shiftedSpikes > duration_s) = shiftedSpikes(shiftedSpikes > duration_s) - duration_s;
            shiftedSpikes = sort(shiftedSpikes);
            shuffledSpikeTimes{chIdx}.(spikeMethod) = shiftedSpikes;
        end
        peakZscore_null(:, shuffleIdx) = computePeakZscoreForEachChannel(shuffledSpikeTimes, ...
            allStimTimes, psth_window_s, baseline_window_s, psth_gaussian_width_ms, ...
            numChannels, spikeMethod);  %#ok<PFBNS>
    end
else
    for shuffleIdx = 1:Nshuffles
        shuffledSpikeTimes = spikeData.spikeTimes;
        for chIdx = 1:numChannels
            chSpikes = shuffledSpikeTimes{chIdx}.(spikeMethod);
            delta = rand() * duration_s;
            shiftedSpikes = chSpikes + delta;
            shiftedSpikes(shiftedSpikes > duration_s) = shiftedSpikes(shiftedSpikes > duration_s) - duration_s;
            shiftedSpikes = sort(shiftedSpikes);
            shuffledSpikeTimes{chIdx}.(spikeMethod) = shiftedSpikes;
        end
        peakZscore_null(:, shuffleIdx) = computePeakZscoreForEachChannel(shuffledSpikeTimes, ...
            allStimTimes, psth_window_s, baseline_window_s, psth_gaussian_width_ms, ...
            numChannels, spikeMethod);
    end
end

%% 3. Determine significance per electrode (two-tailed)
lo_pctile = (alpha / 2) * 100;          % 2.5
hi_pctile = (1 - alpha / 2) * 100;      % 97.5

pctile_lo = prctile(peakZscore_null, lo_pctile, 2);  % [numChannels x 1]
pctile_hi = prctile(peakZscore_null, hi_pctile, 2);

isSigLo = peakZscore_obs < pctile_lo;
isSigHi = peakZscore_obs > pctile_hi;
isSignificant = isSigLo | isSigHi;

%% 4. Package output
shuffleResults.peakZscore_obs   = peakZscore_obs;
shuffleResults.peakZscore_null  = peakZscore_null;
shuffleResults.pctile_lo        = pctile_lo;
shuffleResults.pctile_hi        = pctile_hi;
shuffleResults.isSigLo          = isSigLo;
shuffleResults.isSigHi          = isSigHi;
shuffleResults.isSignificant    = isSignificant;
shuffleResults.Nshuffles        = Nshuffles;
shuffleResults.alpha            = alpha;
shuffleResults.postStimWindow   = [0, psth_window_s(2)];

end
