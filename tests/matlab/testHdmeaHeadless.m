function tests = testHdmeaHeadless
% Focused regression tests for the imported hdMEA/headless workflow.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.originalPath = path;
addpath(root, genpath(fullfile(root, 'Functions')));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function setup(testCase)
testCase.TestData.folder = tempname;
mkdir(testCase.TestData.folder);
end

function teardown(testCase)
rmdir(testCase.TestData.folder, 's');
end

function testBuiltInLayoutAndMarkerSize(testCase)
[channels, coords] = getCoordsFromLayout('Axion64');
verifyEqual(testCase, numel(channels), 64);
verifyEqual(testCase, size(coords), [64 2]);
verifyEqual(testCase, getNodeScaleFactor(coords(:,1), coords(:,2), 2/3, struct()), 2/3);
end

function testCustomLayoutUsesRealCoordinates(testCase)
file = fullfile(testCase.TestData.folder, 'layout.csv');
writetable(table([9; 2], [0; .02], [.1; .2], ...
    'VariableNames', {'channel', 'x_norm', 'y_norm'}), file);
[channels, coords] = getCoordsFromLayout('Custom', file);
verifyEqual(testCase, channels, [9 2]);
verifyEqual(testCase, coords, [0 .8; .16 1.6], 'AbsTol', 1e-12);
end

function testInvalidCustomLayoutIsRejected(testCase)
verifyError(testCase, @() getCoordsFromLayout('Custom'), 'MEANAP:customLayoutMissing');
file = fullfile(testCase.TestData.folder, 'layout.csv');
writetable(table([1; 1], [0; .5], [0; .5], ...
    'VariableNames', {'channel', 'x_norm', 'y_norm'}), file);
verifyError(testCase, @() getCoordsFromLayout('Custom', file), 'MEANAP:customLayoutValues');
end

function testMarkerSpacingHandlesCoincidentAndVerticalNodes(testCase)
verifyEqual(testCase, getNodeScaleFactor([0; 0; .1; .1], [0; 0; 0; 0], 2/3, struct()), ...
    .08, 'AbsTol', 1e-12);
verifyEqual(testCase, getNodeScaleFactor([0; 0], [0; .05], 2/3, struct()), ...
    .04, 'AbsTol', 1e-12);
verifyEqual(testCase, getNodeScaleFactor([0; 0], [0; 0], 2/3, struct()), 2/3);
verifyEqual(testCase, getNodeScaleFactor([0; 1], [0; 1], 2/3, ...
    struct('nodeScaleFactor', .03)), .03);
end

function testSpreadsheetPreservesRecordingChannelOrder(testCase)
folder = testCase.TestData.folder;
customCoordCsv = fullfile(folder, 'layout.csv');
writetable(table([9; 2], [0; .02], [.1; .2], ...
    'VariableNames', {'channel', 'x_norm', 'y_norm'}), customCoordCsv);
channels = [2; 9];
save(fullfile(folder, 'recording.mat'), 'channels');
spreadsheet_filename = fullfile(folder, 'batch.csv');
writetable(table({'recording'}, 6, {'SOHU'}, ...
    'VariableNames', {'Recording', 'Age', 'Group'}), spreadsheet_filename);
spreadsheet_file_type = 'csv';
csvRange = [2 Inf];
rawData = folder;
Params = struct('priorAnalysis', 0, 'suite2pMode', 0, ...
    'channelLayout', 'Custom', 'customCoordCsv', customCoordCsv);
setUpSpreadSheet;
verifyEqual(testCase, Params.channels{1}, [2; 9]);
verifyEqual(testCase, Params.coords{1}, [.16 1.6; 0 .8], 'AbsTol', 1e-12);
verifyEqual(testCase, Params.electrodesToGroundPerRecording, {''});
end

function testPerRecordingCoordinateFiles(testCase)
folder = testCase.TestData.folder;
channels = 1;
save(fullfile(folder, 'first.mat'), 'channels');
save(fullfile(folder, 'second.mat'), 'channels');
for index = 1:2
    writetable(table(1, index/10, index/5, ...
        'VariableNames', {'channel', 'x_norm', 'y_norm'}), ...
        fullfile(folder, sprintf('layout%d.csv', index)));
end
spreadsheet_filename = fullfile(folder, 'batch.csv');
writetable(table({'first'; 'second'}, [6; 6], {'SOHU'; 'SOHU'}, ...
    {''; ''}, {'Custom'; 'Custom'}, {'layout1.csv'; 'layout2.csv'}, ...
    'VariableNames', {'Recording', 'Age', 'Group', 'Ground', 'ChannelLayout', 'CustomCoords'}), ...
    spreadsheet_filename);
spreadsheet_file_type = 'csv';
csvRange = [2 Inf];
rawData = folder;
Params = struct('priorAnalysis', 0, 'suite2pMode', 0, 'channelLayout', 'Custom');
setUpSpreadSheet;
verifyEqual(testCase, Params.coords{1}, [.8 1.6], 'AbsTol', 1e-12);
verifyEqual(testCase, Params.coords{2}, [1.6 3.2], 'AbsTol', 1e-12);
end

function testCustomOutputNameSelectsCurrentRun(testCase)
folder = testCase.TestData.folder;
Info = struct('Grp', {{'SOHU'}}, 'DIV', 6);
NetMet.adjM25mslag.aN = 17;
save(fullfile(folder, 'recording_RetinalRun.mat'), 'Info', 'NetMet');
NetMet.adjM25mslag.aN = 99;
save(fullfile(folder, 'recording_OutputDataOld.mat'), 'Info', 'NetMet');
Params = struct('customGrpOrder', {{}}, 'GrpNm', {{'SOHU'}}, 'DivNm', 6, ...
    'networkLevelNetMetToPlot', {{'aN'}}, 'FuncConLagval', 25, ...
    'lagIndependentMets', {{}}, 'outputDataFolderName', 'RetinalRun');
data = compileAllExpData({'recording'}, folder, Params, 'network');
verifyEqual(testCase, data.SOHU.TP1.aN, 17);
delete(fullfile(folder, 'recording_RetinalRun.mat'));
save(fullfile(folder, 'recording_AnotherRun.mat'), 'Info', 'NetMet');
verifyError(testCase, @() compileAllExpData({'recording'}, folder, Params, 'network'), ...
    'compileAllExpData:multipleExpFiles');
end

function testChannelPlaceholdersRejectStaleMetadata(testCase)
folder = testCase.TestData.folder;
channels = [9; 2];
save(fullfile(folder, 'recording_spikes.mat'), 'channels');
outputFolder = fullfile(folder, 'placeholders');
make_channel_placeholder_mats(folder, outputFolder);
saved = load(fullfile(outputFolder, 'recording.mat'));
verifyEqual(testCase, fieldnames(saved), {'channels'});
verifyEqual(testCase, saved.channels, channels);
make_channel_placeholder_mats(folder, outputFolder);
channels = [2; 9];
save(fullfile(folder, 'recording_spikes.mat'), 'channels');
verifyError(testCase, @() make_channel_placeholder_mats(folder, outputFolder), ...
    'MEANAP:channelPlaceholderMismatch');
end

function testConverterPreservesTimesAndRecordsAnalysisRate(testCase)
folder = testCase.TestData.folder;
resultFile = fullfile(folder, 'fixture.raw.result.hdf5');
clustersFile = fullfile(folder, 'fixture.raw.clusters.hdf5');
h5create(resultFile, '/spiketimes/temp_0', [4 1]);
h5write(resultFile, '/spiketimes/temp_0', [1; 120; 35000; 51000]);
h5create(resultFile, '/spiketimes/temp_1', [3 1]);
h5write(resultFile, '/spiketimes/temp_1', [2; 1700; 33000]);
h5create(clustersFile, '/electrodes', [2 1]);
h5write(clustersFile, '/electrodes', [0; 1]);
out = spykingcircus_to_meanap_v3(folder, 'Fs', 20000, 'Duration', 2.75, ...
    'AgeGroup', 8, 'Group', 'Control');
saved = load(out.saveName);
verifyEqual(testCase, saved.spikeTimes{1}.spykingcircus, [1; 120; 35000]/20000);
verifyEqual(testCase, saved.spikeDetectionResult.params.fs, 1000);
verifyEqual(testCase, saved.spikeDetectionResult.params.acquisitionFs, 20000);
verifyEqual(testCase, saved.spikeDetectionResult.params.duration, 2);
batch = readtable(out.batchCsv);
verifyEqual(testCase, batch{1,2}, 8);
verifyEqual(testCase, string(batch{1,3}), "Control");
verifyTrue(testCase, verify_meanap_spikefile(out.saveName));
end

function testHeadlessConfigurationErrorsAreClear(testCase)
verifyError(testCase, @() MEApipeline_headless_retinal(), 'MEANAP:headlessConfigRequired');
verifyError(testCase, @() MEApipeline_headless_retinal(struct()), 'MEANAP:headlessConfigMissing');
end
