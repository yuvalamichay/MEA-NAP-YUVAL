function [channels,coords] = getCoordsFromLayout(channelLayout, customCoordCsv)
%GETCOORDSFROMLAYOUT Specifies mapping from channel number to spatial
%coordinates
% We provide a few default coordinates from 
% the Multichannel systems (MCS) (60 electrode In Vitro MEA set up) and
% the Axion system. Alternatively, you can specify your own layout 
% by specifying giving the input arugment 'Custom' and editing 
% the line after: 'elseif strcmp(channelLayout, 'Custom')'
% Parameters 
% ----------
% channelLayout : str 
%      channel coordinate to use 
%      'MCS60old' : old version of MCS 60 electrode MEA 
%      'MCS60' : lateset version of MCS 60 electrode MEA 
%      'MCS59' : MCS 60 electrode MEA with grounded electrode removed 
%      'Axion64' : Axion 64 electrode MEA (6-well plate)
%      'Axion16' : Axion 16 electrode MEA (48-well plate)
%      'Custom' : specify your own channel layout
if strcmp(channelLayout, 'MCS60old')

    channels = [11, 12, 13, 14, 15, 16, 17, 18, ... 
            21, 22, 23, 24, 25, 26, 27, 28, ...
            31, 32, 33, 34, 35, 36, 37, 38, ...
            41, 42, 43, 44, 45, 46, 47, 48, ...
            51, 52, 53, 54, 55, 56, 57, 58, ...
            61, 62, 63, 64, 65, 66, 67, 68, ...,
            71, 72, 73, 74, 75, 76, 77, 78, ...,
            81, 82, 83, 84, 85, 86, 87, 88];

    channelsOrdering = ...
    [21 31 41 51 61 71 12 22 32 42 52 62 72 82 13 23 33 43 53 63 ... 
    73 83 14 24 34 44 54 64 74 84 15 25 35 45 55 65 75 85 16 26 ...
    36 46 56 66 76 86 17 27 37 47 57 67 77 87 28 38 48 58 68 78];

    Params.coords = zeros(length(channels), 2);
    Params.coords(:, 2) = repmat(linspace(1, 0, 8), 1, 8);
    Params.coords(:, 1) = repelem(linspace(0, 1, 8), 1, 8);

    subset_idx = find(~ismember(channels, [11, 81, 18, 88]));
    channels = channels(subset_idx);

    reorderingIdx = zeros(length(channels), 1);
    for n = 1:length(channels)
        reorderingIdx(n) = find(channelsOrdering(n) == channels);
    end 
    
    Params.coords = Params.coords(subset_idx, :);

    % Re-order the channel IDs and coordinates to match the original
    % ordering
    channels = channels(reorderingIdx); 
    coords = Params.coords(reorderingIdx, :);

elseif strcmp(channelLayout, 'MCS60')

     channels = [11, 12, 13, 14, 15, 16, 17, 18, ... 
            21, 22, 23, 24, 25, 26, 27, 28, ...
            31, 32, 33, 34, 35, 36, 37, 38, ...
            41, 42, 43, 44, 45, 46, 47, 48, ...
            51, 52, 53, 54, 55, 56, 57, 58, ...
            61, 62, 63, 64, 65, 66, 67, 68, ...,
            71, 72, 73, 74, 75, 76, 77, 78, ...,
            81, 82, 83, 84, 85, 86, 87, 88];

    channelsOrdering = ...
    [47, 48, 46, 45, 38, 37, 28, 36, 27, 17, 26, 16, 35, 25, ...
    15, 14, 24, 34, 13, 23, 12, 22, 33, 21, 32, 31, 44, 43, 41, 42, ...
    52, 51, 53, 54, 61, 62, 71, 63, 72, 82, 73, 83, 64, 74, 84, 85, 75, ...
    65, 86, 76, 87, 77, 66, 78, 67, 68, 55, 56, 58, 57];

    Params.coords = zeros(length(channels), 2);
    Params.coords(:, 2) = repmat(linspace(1, 0, 8), 1, 8);
    Params.coords(:, 1) = repelem(linspace(0, 1, 8), 1, 8);

    subset_idx = find(~ismember(channels, [11, 81, 18, 88]));
    channels = channels(subset_idx);

    reorderingIdx = zeros(length(channels), 1);
    for n = 1:length(channels)
        reorderingIdx(n) = find(channelsOrdering(n) == channels);
    end 
    
    Params.coords = Params.coords(subset_idx, :);

    % Re-order the channel IDs and coordinates to match the original
    % ordering
    channels = channels(reorderingIdx);
    coords = Params.coords(reorderingIdx, :);
    % Params.reorderingIdx = reorderingIdx;

elseif strcmp(channelLayout, 'MCS59')

    channels = [11, 12, 13, 14, 15, 16, 17, 18, ... 
            21, 22, 23, 24, 25, 26, 27, 28, ...
            31, 32, 33, 34, 35, 36, 37, 38, ...
            41, 42, 43, 44, 45, 46, 47, 48, ...
            51, 52, 53, 54, 55, 56, 57, 58, ...
            61, 62, 63, 64, 65, 66, 67, 68, ...,
            71, 72, 73, 74, 75, 76, 77, 78, ...,
            81, 82, 83, 84, 85, 86, 87, 88];

    channelsOrdering = ...
    [47, 48, 46, 45, 38, 37, 28, 36, 27, 17, 26, 16, 35, 25, ...
    15, 14, 24, 34, 13, 23, 12, 22, 33, 21, 32, 31, 44, 43, 41, 42, ...
    52, 51, 53, 54, 61, 62, 71, 63, 72, 82, 73, 83, 64, 74, 84, 85, 75, ...
    65, 86, 76, 87, 77, 66, 78, 67, 68, 55, 56, 58, 57];

    Params.coords = zeros(length(channels), 2);
    Params.coords(:, 2) = repmat(linspace(1, 0, 8), 1, 8);
    Params.coords(:, 1) = repelem(linspace(0, 1, 8), 1, 8);

    subset_idx = find(~ismember(channels, [11, 81, 18, 88]));
    channels = channels(subset_idx);

    reorderingIdx = zeros(length(channels), 1);
    for n = 1:length(channels)
        reorderingIdx(n) = find(channelsOrdering(n) == channels);
    end 
    
    Params.coords = Params.coords(subset_idx, :);

    % Re-order the channel IDs and coordinates to match the original
    % ordering
    channels = channels(reorderingIdx);
    Params.channels = channels; 
    Params.coords = Params.coords(reorderingIdx, :);

    inclusionIndex = find(channelsOrdering ~= 82);
    channels = channels(inclusionIndex);
    coords = Params.coords(inclusionIndex, :);
    % Params.reorderingIdx = reorderingIdx; % (inclusionIndex)


elseif strcmp(channelLayout, 'Axion64')
    
    channels = []; 
    coords = zeros(length(channels), 2);
    
    numRow = 8;
    numCol = 8;
    spacing_vals = linspace(0, 1, numRow);
    
    channelIdx = 1;
    for rowIdx = 1:numRow 
        for colIdx = 1:numCol 
           channels(channelIdx) = colIdx * 10 + rowIdx;
           coords(channelIdx, 1) = spacing_vals(colIdx);
           coords(channelIdx, 2) = spacing_vals(rowIdx);
           channelIdx = channelIdx + 1;
       end 
    end 

elseif strcmp(channelLayout, 'Axion16')    
    
    channels = []; 
    coords = zeros(length(channels), 2);
    
    numRow = 4;
    numCol = 4;
    spacing_vals = linspace(0, 1, numRow);
    
    channelIdx = 1;
    for rowIdx = 1:numRow 
        for colIdx = 1:numCol 
           channels(channelIdx) = colIdx * 10 + rowIdx;
           coords(channelIdx, 1) = spacing_vals(colIdx);
           coords(channelIdx, 2) = spacing_vals(rowIdx);
           channelIdx = channelIdx + 1;
       end 
    end 

elseif strcmp(channelLayout, 'Custom')
    % Coordinates are normalized to 0-1 here; the common x8 scaling below
    % converts them to MEA-NAP plotting units.
    if nargin < 2 || isempty(customCoordCsv) || ~isfile(customCoordCsv)
        error('MEANAP:customLayoutMissing', ...
            'Set Params.customCoordCsv to a CSV with channel, x_norm and y_norm columns.');
    end
    coordTable = readtable(customCoordCsv);
    required = {'channel', 'x_norm', 'y_norm'};
    if ~all(ismember(required, coordTable.Properties.VariableNames))
        error('MEANAP:customLayoutColumns', ...
            'Custom layout CSV must contain channel, x_norm and y_norm columns.');
    end
    if ~isnumeric(coordTable.channel) || ~isnumeric(coordTable.x_norm) ...
            || ~isnumeric(coordTable.y_norm)
        error('MEANAP:customLayoutValues', 'Custom layout columns must be numeric.');
    end
    channels = coordTable.channel(:)';
    coords = [coordTable.x_norm, coordTable.y_norm];
    if isempty(channels) ...
            || any(~isfinite(channels)) || any(~isfinite(coords(:))) ...
            || numel(unique(channels)) ~= numel(channels)
        error('MEANAP:customLayoutValues', ...
            'Custom layout needs unique finite numeric channel IDs and finite coordinates.');
    end
    if any(coords(:) < 0 | coords(:) > 1)
        error('MEANAP:customLayoutRange', 'Custom x_norm and y_norm must be between 0 and 1.');
    end

end 

coords  = coords * 8;  % Do not remove this line after specifying coordinate positions in (0 - 1 format)

end

