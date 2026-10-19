function make_channel_placeholder_mats(spikeDataFolder, outputFolder)
%MAKE_CHANNEL_PLACEHOLDER_MATS Supply channel metadata without raw voltages.
% Creates <recording>.mat containing only channels from <recording>_spikes.mat.
% Existing files are reused only when their channel lists agree.

spikeFiles = dir(fullfile(spikeDataFolder, '*_spikes.mat'));
if isempty(spikeFiles)
    error('MEANAP:noSpikeFiles', 'No *_spikes.mat files found in %s.', spikeDataFolder);
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
for index = 1:numel(spikeFiles)
    source = load(fullfile(spikeDataFolder, spikeFiles(index).name), 'channels');
    if ~isfield(source, 'channels')
        error('MEANAP:channelsMissing', '%s has no channels variable.', spikeFiles(index).name);
    end
    recordingName = erase(spikeFiles(index).name, '_spikes.mat');
    outputPath = fullfile(outputFolder, [recordingName, '.mat']);
    if isfile(outputPath)
        existing = load(outputPath, 'channels');
        if ~isfield(existing, 'channels') || ~isequal(existing.channels(:), source.channels(:))
            error('MEANAP:channelPlaceholderMismatch', ...
                'Channel metadata differs in %s; choose a new placeholder folder.', outputPath);
        end
    else
        channels = source.channels;
        save(outputPath, 'channels');
    end
end
end
