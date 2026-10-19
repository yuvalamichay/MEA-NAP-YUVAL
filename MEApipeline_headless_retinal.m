function MEApipeline_headless_retinal(config)
%MEAPIPELINE_HEADLESS_RETINAL Run MATLAB Steps 2-4 on converted hdMEA spikes.
% config is a scalar struct, or a MAT file containing config or Params.
% Required fields: spikeDetectedData and spreadSheetFileName (batch CSV).
% Optional fields include customCoordCsv, rawData, outputDataFolder,
% outputDataFolderName, FuncConLagval, cartographyLagVal and AdvancedSettings
% overrides. Paths should be absolute. See docs/hdmea-headless.md.

if nargin == 0
    error('MEANAP:headlessConfigRequired', ...
        'Supply a configuration struct or MAT file; see docs/hdmea-headless.md.');
end
if ischar(config) || (isstring(config) && isscalar(config))
    saved = load(config);
    if isfield(saved, 'config')
        config = saved.config;
    elseif isfield(saved, 'Params')
        config = saved.Params;
    else
        error('MEANAP:headlessConfigMissing', 'MAT file must contain config or Params.');
    end
end
if ~isstruct(config) || ~isscalar(config)
    error('MEANAP:headlessConfigType', 'Configuration must be a scalar struct.');
end
required = {'spikeDetectedData', 'spreadSheetFileName'};
for index = 1:numel(required)
    if ~isfield(config, required{index}) || isempty(config.(required{index}))
        error('MEANAP:headlessConfigMissing', 'Set config.%s.', required{index});
    end
end
if ~isfolder(config.spikeDetectedData) || ~isfile(config.spreadSheetFileName)
    error('MEANAP:headlessInputMissing', 'Spike folder and batch CSV must exist.');
end
[~, ~, extension] = fileparts(config.spreadSheetFileName);
if ~strcmpi(extension, '.csv')
    error('MEANAP:headlessBatchType', 'The hdMEA headless runner requires a batch CSV.');
end
spikeFiles = dir(fullfile(config.spikeDetectedData, '*_spikes.mat'));
if isempty(spikeFiles)
    error('MEANAP:noSpikeFiles', 'No *_spikes.mat files found in %s.', config.spikeDetectedData);
end
metadata = load(fullfile(config.spikeDetectedData, spikeFiles(1).name), 'spikeDetectionResult');
if isfield(metadata.spikeDetectionResult, 'params')
    spikeParams = metadata.spikeDetectionResult.params;
else
    spikeParams = metadata.spikeDetectionResult.Params;
end
validateattributes(spikeParams.fs, {'numeric'}, {'scalar', 'finite', 'positive'});

Params = struct();
Params.outputDataFolder = fullfile(fileparts(config.spikeDetectedData), 'MEANAP_output');
Params.outputDataFolderName = ['OutputData_hdMEA_', datestr(now, 'yyyymmdd_HHMMSS')];
Params.customCoordCsv = fullfile(config.spikeDetectedData, 'meanap_custom_coords.csv');
% A per-recording CustomCoords column can supply the initial layout too.
if ~isfield(config, 'customCoordCsv') && ~isfile(Params.customCoordCsv)
    importOptions = detectImportOptions(config.spreadSheetFileName);
    importOptions.VariableNamingRule = 'modify';
    if isfield(config, 'csvRange')
        importOptions.DataLines = config.csvRange;
    end
    batch = readtable(config.spreadSheetFileName, importOptions);
    if ismember('CustomCoords', batch.Properties.VariableNames) && height(batch) > 0
        Params.customCoordCsv = char(string(batch.CustomCoords(1)));
        if ~isfile(Params.customCoordCsv)
            Params.customCoordCsv = fullfile(fileparts(config.spreadSheetFileName), Params.customCoordCsv);
        end
    end
end
Params.channelLayout = 'Custom';
Params.fs = spikeParams.fs;
Params.dSampF = spikeParams.fs;
Params.potentialDifferenceUnit = 'V';
Params.SpikesMethod = 'spykingcircus';
Params.FuncConLagval = 25;
Params.cartographyLagVal = 25;
Params.ProbThreshPlotChecks = 0;
Params.figExt = {'.png'};
Params.fullSVG = 0;
Params.showOneFig = 1;
Params.timeProcesses = 1;
Params.minActivityLevel = 0.01;
Params.singleChannelIsiThreshold = 'automatic';
Params.use_theoretical_bounds = 0;
Params.minNodeSize = 0.01;
Params.maxNodeSize = 1;
Params.nodeScalingMethod = 'Linear';
Params.nodeScalingPower = 1;
Params.networkPlotEdgeThresholdMethod = 'Absolute Value';
Params.networkPlotEdgeThresholdPercentile = 95;
Params.networkPlotEdgeThreshold = 0.001;
Params.maxNumEdgesToPlot = 10000;
Params.edgeSubsamplingMethod = 'highToLow';
Params.nodeLayout = 'Original';
% The default AdvancedSettings list requests PC_raw/Cmcblty even though
% their implementation is commented out in ExtractNetMet.
Params.netMetToCal = {'aN', 'Dens', 'NDmean', 'NDtop25', 'sigEdgesMean', ...
    'sigEdgesTop10', 'CC', 'nMod', 'Q', 'PL', 'Eglob', 'SW', 'SWw', ...
    'effRank', 'num_nnmf_components', 'nComponentsRelNS', 'NSmean', ...
    'ElocMean', 'PCmean', 'PCmeanTop10', 'PCmeanBottom10', ...
    'percentZscoreGreaterThanZero', 'percentZscoreLessThanZero', ...
    'NCpn1', 'NCpn2', 'NCpn3', 'NCpn4', 'NCpn5', 'NCpn6', 'ND', 'NS', ...
    'Z', 'Eloc', 'PC', 'BC', 'MEW', 'aveControl', 'modalControl', 'aveControlMean'};

fields = fieldnames(config);
for index = 1:numel(fields)
    Params.(fields{index}) = config.(fields{index});
end
% These define this entry point. The branch's current MEApipeline performs
% the actual analysis, including its current burst-detection improvements.
Params.HomeDir = fileparts(mfilename('fullpath'));
Params.guiMode = 0;
Params.detectSpikes = 0;
Params.runSpikeCheckOnPrevSpikeData = 0;
Params.suite2pMode = 0;
Params.stimulationMode = 0;
Params.priorAnalysis = 0;
Params.startAnalysisStep = 2;
if ~isfield(Params, 'rawData') || isempty(Params.rawData)
    Params.rawData = fullfile(Params.outputDataFolder, ...
        Params.outputDataFolderName, 'ChannelPlaceholders');
    make_channel_placeholder_mats(Params.spikeDetectedData, Params.rawData);
end

originalFolder = pwd;
originalPath = path;
originalRenderer = get(groot, 'DefaultFigureRenderer');
originalVisibility = get(groot, 'DefaultFigureVisible');
cleanup = onCleanup(@() restoreSession(originalFolder, originalPath, ...
    originalRenderer, originalVisibility)); %#ok<NASGU>
set(groot, 'DefaultFigureRenderer', 'painters', 'DefaultFigureVisible', 'off');
cd(Params.HomeDir);
MEApipeline(Params);
end

function restoreSession(folder, originalPath, renderer, visibility)
cd(folder);
path(originalPath);
set(groot, 'DefaultFigureRenderer', renderer, 'DefaultFigureVisible', visibility);
end
