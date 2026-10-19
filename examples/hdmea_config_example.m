function config = hdmea_config_example(dataFolder, outputFolder)
%HDMEA_CONFIG_EXAMPLE Configuration for converted retinal hdMEA recordings.
% Add the repository root and examples/ to your MATLAB path, then:
%   config = hdmea_config_example('/path/to/Retinal MEA analysis');
%   MEApipeline_headless_retinal(config);

if nargin < 2
    outputFolder = fullfile(dataFolder, 'MEANAP_output');
end
config = struct();
config.spikeDetectedData = fullfile(dataFolder, 'MEANAP_spikes');
config.spreadSheetFileName = fullfile(config.spikeDetectedData, 'MEANAP_batch.csv');
config.customCoordCsv = fullfile(config.spikeDetectedData, 'meanap_custom_coords.csv');
config.outputDataFolder = outputFolder;
config.outputDataFolderName = ['OutputData_hdMEA_', datestr(now, 'yyyymmdd_HHMMSS')];
config.FuncConLagval = 25;
config.cartographyLagVal = config.FuncConLagval;
config.ProbThreshRepNum = 200;
config.ProbThreshTail = 0.05;
% The spike file supplies the matrix sampling rate. Do not set acquisition
% frequency here when the converter used a separate AnalysisFs.
end
