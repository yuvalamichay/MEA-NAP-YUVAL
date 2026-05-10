function plotStimShuffleResults(shuffleResults, spikeData, Info, Params, figFolder, patternIdx)
% PLOTSTIMSHUFFLERESULTS  Visualise results of the circular-shift shuffle test.
%
% Generates two figures per pattern:
%   1. Null distribution heatmap with observed AUC overlay and errorbar plot.
%   2. Electrode layout heatmap coloured by significance.
%
% INPUTS
% ------
% shuffleResults : struct   - output of stimShuffleTest
% spikeData      : struct   - spike data (for electrode count / layout info)
% Info           : struct   - experiment info (for title)
% Params         : struct   - pipeline parameters
% figFolder      : char     - path to save figures
% patternIdx     : scalar   - (optional) pattern index for labelling.
%                             If empty or not provided, labelled as 'all stim'.
%
% See also: stimShuffleTest, pipelineSaveFig

if nargin < 6
    patternIdx = [];
end

if isempty(patternIdx)
    nameSuffix = '_allStim';
    titleSuffix = 'all stim';
else
    nameSuffix = sprintf('_pattern_%.f', patternIdx);
    titleSuffix = sprintf('pattern %.f', patternIdx);
end

numChannels = length(shuffleResults.peakZscore_obs);
peakZscore_obs = shuffleResults.peakZscore_obs;
pctile_thresh = shuffleResults.pctile_thresh;
isSignificant = shuffleResults.isSignificant;
xVec = 1:numChannels;

%% Figure 1: Null distribution heatmap with observed peak z-score overlay
fig1 = figure('Position', [100, 100, 800, 500], 'Color', 'w', 'Visible', 'off');

% Sort null distributions for each electrode for visualisation
peakZscore_null_sorted = sort(shuffleResults.peakZscore_null, 2);

subplot(2, 1, 1)
imagesc(1:shuffleResults.Nshuffles, 1:numChannels, peakZscore_null_sorted)
colorbar
ylabel('Channel Index')
xlabel('Shuffle (sorted)')
title(sprintf('Null peak z-score distribution – %s  |  %s', titleSuffix, char(Info.FN{1})), ...
    'Interpreter', 'none')
set(gca, 'TickDir', 'out')
box off

subplot(2, 1, 2)
hold on
% Plot null median and 95th percentile threshold per electrode
nullMedian = median(shuffleResults.peakZscore_null, 2);
errorbar(xVec, nullMedian, zeros(size(nullMedian)), pctile_thresh - nullMedian, ...
    '.', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1, 'DisplayName', 'Null median \rightarrow 95th pctile');
% Overlay observed
scatter(xVec(~isSignificant), peakZscore_obs(~isSignificant), 25, ...
    [0.3, 0.3, 0.7], 'filled', 'DisplayName', 'Not significant');
scatter(xVec(isSignificant), peakZscore_obs(isSignificant), 40, ...
    [0.9, 0.2, 0.2], 'filled', 'DisplayName', 'Significant');
hold off
xlabel('Channel Index')
ylabel('Peak Z-score')
legend('Location', 'best')
title(sprintf('Observed vs null – %s', titleSuffix))
box off
set(gca, 'TickDir', 'out')

figName = sprintf('12_shuffle_test_null_dist%s', nameSuffix);
pipelineSaveFig(fullfile(figFolder, figName), Params.figExt, Params.fullSVG, fig1);
close(fig1);

%% Figure 2: Electrode heatmap of significance (if plotStimHeatmapWmetric exists)
if exist('plotStimHeatmapWmetric', 'file') && ~isempty(patternIdx)
    fig2 = figure('Visible', 'off');
    nodeMetric = double(isSignificant);
    vrange = [0, 1];
    cmap = [0.85, 0.85, 0.85; 0.9, 0.2, 0.2];  % grey = not sig, red = sig
    cmapLabel = 'Significant';
    fig2 = plotStimHeatmapWmetric(nodeMetric, vrange, cmap, cmapLabel, ...
        spikeData.stimInfo, patternIdx, fig2);
    title(sprintf('Significant electrodes – %s', titleSuffix));
    figName = sprintf('12_shuffle_test_sig_heatmap%s', nameSuffix);
    figPath = fullfile(figFolder, figName);
    pipelineSaveFig(figPath, Params.figExt, Params.fullSVG, fig2);
    close(fig2);
end

end
