function shuffleResults = stimShuffleTest(spikeData, allStimTimes, Params, Info)
% STIMSHUFFLETEST  Circular-shift shuffle test for post-stimulus AUC significance.
%
% Builds a null distribution of post-stimulus AUC values by circularly
% shifting each electrode's spike times by a random offset (with
% wrap-around within the recording duration), then recomputing the
% trial-averaged post-stimulus firing rate and its AUC.  Significance is
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
%       .shuffleBinWidth       - bin width in seconds for the shuffle PSTH
%                                histogram.  Default: 0.002 (2 ms).
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
%   .AUC_obs           - [numChannels x 1] observed AUC for each electrode
%   .AUC_null          - [numChannels x Nshuffles] null AUC distributions
%   .pctile_lo         - [numChannels x 1] lower percentile bound (2.5th)
%   .pctile_hi         - [numChannels x 1] upper percentile bound (97.5th)
%   .isSigLo           - [numChannels x 1] logical, true if AUC_obs < pctile_lo
%   .isSigHi           - [numChannels x 1] logical, true if AUC_obs > pctile_hi
%   .isSignificant     - [numChannels x 1] logical, true if significant (either tail)
%   .Nshuffles         - scalar, number of shuffles performed
%   .alpha             - scalar, significance level used
%   .postStimWindow    - [1 x 2] the post-stimulus window used for AUC
%   .binWidth          - scalar, bin width used (s)
%
% PROCEDURE
% ---------
% 1. Compute observed post-stimulus PSTH and AUC for each electrode.
% 2. For each shuffle iteration:
%      a. For each electrode, circularly shift its spike train by a random
%         offset uniformly drawn from (0, recordingDuration) with wrap-around.
%      b. Recompute the post-stimulus PSTH and AUC for every electrode
%         using the shifted spike times aligned to the *original* stim times.
% 3. For each electrode, determine the 2.5th and 97.5th percentiles of its
%    null AUC distribution.
% 4. Mark an electrode as significant if observed AUC falls outside the
%    null percentile interval.
%
% REFERENCE
% ---------
% Circular-shift approach inspired by adjM_thr_parallel.m (STTC thresholding).
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
numStimEvents = length(allStimTimes);
duration_s = Info.duration_s;

% Post-stimulus window: from 0 to the end of stimAnalysisWindow
postStimWindow = [0, Params.stimAnalysisWindow(2)];

% Bin width for the shuffle PSTH (independent of Params.rasterBinWidth)
if isfield(Params, 'shuffleBinWidth') && ~isempty(Params.shuffleBinWidth)
    binWidth = Params.shuffleBinWidth;
else
    binWidth = 0.002;  % default: 2 ms
end

rasterBins = postStimWindow(1):binWidth:postStimWindow(2);
numBins = length(rasterBins) - 1;

%% 1. Compute observed AUC
spikeMethod = Params.SpikesMethod;
AUC_obs = computePostStimAUC(spikeData.spikeTimes, allStimTimes, ...
    rasterBins, binWidth, numChannels, numStimEvents, numBins, spikeMethod);

%% 2. Build null distribution via circular shift
AUC_null = zeros(numChannels, Nshuffles);

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
        AUC_null(:, shuffleIdx) = computePostStimAUC(shuffledSpikeTimes, ...
            allStimTimes, rasterBins, binWidth, numChannels, numStimEvents, numBins, spikeMethod);  %#ok<PFBNS>
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
        AUC_null(:, shuffleIdx) = computePostStimAUC(shuffledSpikeTimes, ...
            allStimTimes, rasterBins, binWidth, numChannels, numStimEvents, numBins, spikeMethod);
    end
end

%% 3. Determine significance per electrode (two-tailed)
lo_pctile = (alpha / 2) * 100;          % 2.5
hi_pctile = (1 - alpha / 2) * 100;      % 97.5

pctile_lo = prctile(AUC_null, lo_pctile, 2);  % [numChannels x 1]
pctile_hi = prctile(AUC_null, hi_pctile, 2);

isSigLo = AUC_obs < pctile_lo;
isSigHi = AUC_obs > pctile_hi;
isSignificant = isSigLo | isSigHi;

%% 4. Package output
shuffleResults.AUC_obs       = AUC_obs;
shuffleResults.AUC_null      = AUC_null;
shuffleResults.pctile_lo     = pctile_lo;
shuffleResults.pctile_hi     = pctile_hi;
shuffleResults.isSigLo       = isSigLo;
shuffleResults.isSigHi       = isSigHi;
shuffleResults.isSignificant = isSignificant;
shuffleResults.Nshuffles     = Nshuffles;
shuffleResults.alpha         = alpha;
shuffleResults.postStimWindow = postStimWindow;
shuffleResults.binWidth      = binWidth;

end


%% ========================================================================
%  Local function: compute per-electrode AUC from spike times
%  ========================================================================
function aucVec = computePostStimAUC(spikeTimes_cell, allStimTimes, ...
    rasterBins, binWidth, numChannels, numStimEvents, numBins, spikeMethod)
% Computes the area under the mean PSTH curve for each electrode.
%
% For each electrode, the firing rate is binned around each stimulus event,
% averaged across trials, and integrated (trapz) to give the AUC.

    aucVec = zeros(numChannels, 1);
    binCentres = rasterBins(1:end-1) + binWidth / 2;

    for chIdx = 1:numChannels
        chSpikes = spikeTimes_cell{chIdx}.(spikeMethod);
        % Firing rate aligned to each stim event: (numStimEvents x numBins)
        frMat = zeros(numStimEvents, numBins);
        for evIdx = 1:numStimEvents
            frMat(evIdx, :) = histcounts(chSpikes - allStimTimes(evIdx), rasterBins) / binWidth;
        end
        meanFR = mean(frMat, 1);  % average across trials
        aucVec(chIdx) = trapz(binCentres, meanFR);
    end
end
