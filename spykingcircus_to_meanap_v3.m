function out = spykingcircus_to_meanap_v3(folder, varargin)
% SPYKINGCIRCUS_TO_MEANAP_V3  Convert SpyKING CIRCUS output into a MEA-NAP spike file.
%
% v3 changes (vs v2), all aimed at how MEA-NAP itself handles the file:
%   - 'AnalysisFs' (default 1000 Hz) is written as the file's fs. MEA-NAP
%     bins spikes into a (samples x nodes) matrix at this rate and converts
%     it to a full matrix in several places (raster, effective rank, NMF).
%     At 20 kHz with 306 nodes that is ~6.2 GB each time; at 1 kHz ~0.3 GB.
%     Spike times keep full 0.05 ms precision, and STTC uses the times
%     themselves. The analysis rate also sets the surrogate circular-shift
%     grid, so thresholded connectivity need not be numerically identical.
%   - 'WholeSeconds' (default true): MEA-NAP floors the duration to whole
%     seconds; spikes after that point are dropped here too, so spike times
%     and spike matrix agree (here: the last 0.68 s, 5,758 spikes).
%   - 'FlipY' (default true): MEA-NAP draws y upwards, the MaxWell chip
%     counts y downwards; flipping keeps MEA-NAP maps the same way up as
%     MaxLab's electrode view.
%
% v2 changes (vs spykingcircus_to_meanap.m):
%   - duration taken from the raw frame count (126.68 s here); /info/duration
%     in the result file is rounded down to whole seconds (126)
%   - node coordinates normalised with ONE scale factor for x and y, so the
%     electrode layout keeps its true aspect ratio in MEA-NAP plots
%   - coordinates already in MEA-NAP units (0-8) are saved as coordsMEANAP
%   - stimulus TTL events (/bits) are exported as stimEvents, in seconds on
%     the same clock as the spike times (needs only deflate, no MaxWell plugin)
%
% Written against the actual layout of this dataset (MaxWell MaxOne/MaxTwo
% recording, 941 routed channels, sorted with SpyKING CIRCUS v1, 306 units):
%
%   *.result.hdf5    /spiketimes/temp_K   spike times, in samples
%                    /amplitudes/temp_K   2 x nSpikes fitted amplitudes
%                    /info/duration       recording duration
%   *.clusters.hdf5  /electrodes          peak channel per unit (0-based)
%   *.raw.h5         /data_store/data0000/settings/sampling   sampling rate
%                    /data_store/data0000/settings/mapping    channel -> electrode, x, y
%                    /data_store/data0000/groups/routed/channels
%
% Usage:
%   out = spykingcircus_to_meanap_v3(pwd)
%   out = spykingcircus_to_meanap_v3(pwd, 'NodeMode', 'electrode')
%   out = spykingcircus_to_meanap_v3(pwd, 'MinSpikes', 50)
%   out = spykingcircus_to_meanap_v3(pwd, 'AnalysisFs', 20000)  % v2 behaviour
%
% Options (name/value):
%   'Base'       base recording name; default: from the result file, '.raw' removed
%   'NodeMode'   'unit' (default) = one MEA-NAP node per sorted unit
%                'electrode'      = pool units that share a peak channel
%   'Method'     spike-times field name; default 'spykingcircus'. Set MEA-NAP's
%                "Spike method for analysis" to this exact string.
%   'MinSpikes'  drop units with fewer than this many spikes; default 0
%   'Fs'         sampling rate in Hz; default: read from the raw file
%   'Duration'   duration in seconds; default: raw frame count / Fs
%   'OutFolder'  output folder; default: a 'MEANAP_spikes' subfolder
%   'AnalysisFs' sampling rate written for MEA-NAP's spike matrix; default 1000
%   'WholeSeconds' trim to floor(duration) as MEA-NAP does; default true
%   'FlipY'      flip y so MEA-NAP plots match MaxLab's view; default true
%   'AgeGroup'   batch metadata; default 6 (original retinal dataset)
%   'Group'      batch group name; default 'SOHU' (original retinal dataset)
%
% Writes:
%   <base>_spikes.mat            spikeTimes, channels, spikeDetectionResult,
%                                spikeWaveforms, spikeAmps, thresholds, coords
%   meanap_custom_coords.csv     node coordinates for the custom channel layout
%   meanap_custom_layout.m       drop-in replacement for the Custom branch of
%                                MEA-NAP's getCoordsFromLayout.m
%   MEANAP_batch.csv             batch spreadsheet with one row for this recording

%% ---------------- options ----------------
if nargin < 1 || isempty(folder)
    folder = pwd;
end

p = inputParser;
addParameter(p, 'Base', '');
addParameter(p, 'NodeMode', 'unit');
addParameter(p, 'Method', 'spykingcircus');
addParameter(p, 'MinSpikes', 0);
addParameter(p, 'Fs', []);
addParameter(p, 'Duration', []);
addParameter(p, 'OutFolder', '');
addParameter(p, 'AnalysisFs', 1000);
addParameter(p, 'WholeSeconds', true);
addParameter(p, 'FlipY', true);
addParameter(p, 'AgeGroup', 6);
addParameter(p, 'Group', 'SOHU');
parse(p, varargin{:});
opt = p.Results;
validateattributes(opt.AnalysisFs, {'numeric'}, {'scalar', 'finite', 'positive', 'integer'});
validateattributes(opt.AgeGroup, {'numeric'}, {'scalar', 'finite'});
if ~isvarname(char(opt.Group))
    error('MEANAP:invalidGroup', 'Group must be a valid MATLAB field name.');
end

method = opt.Method;
if ~isvarname(method)
    error('Method must be a valid MATLAB field name.');
end

%% ---------------- locate files ----------------
% Prefer the merged variants, but fall back when they are absent.
resultFile   = pickFile(folder, {'*.result-merged.hdf5', '*.result.hdf5'},   'result');
clustersFile = pickFile(folder, {'*.clusters-merged.hdf5', '*.clusters.hdf5'}, 'clusters');
rawFiles     = dir(fullfile(folder, '*.raw.h5'));
if isempty(rawFiles)
    rawFile = '';
else
    rawFile = fullfile(rawFiles(1).folder, rawFiles(1).name);
end

base = opt.Base;
if isempty(base)
    [~, n] = fileparts(resultFile);
    n = regexprep(n, '\.result(-merged)?$', '');
    n = regexprep(n, '\.raw$', '');
    base = n;
end

outFolder = opt.OutFolder;
if isempty(outFolder)
    outFolder = fullfile(folder, 'MEANAP_spikes');
end
if ~isfolder(outFolder)
    mkdir(outFolder);
end

fprintf('Result file    : %s\n', shortName(resultFile));
fprintf('Clusters file  : %s\n', shortName(clustersFile));
if isempty(rawFile)
    fprintf('Raw file       : (none - pass Fs and Duration yourself)\n');
else
    fprintf('Raw file       : %s\n', shortName(rawFile));
end
fprintf('Base name      : %s\n', base);

%% ---------------- sampling rate ----------------
fs = opt.Fs;
if isempty(fs) && ~isempty(rawFile)
    fs = readFirstAvailable(rawFile, { ...
        '/data_store/data0000/settings/sampling', ...
        '/recordings/rec0000/well000/settings/sampling', ...
        '/settings/sampling'});
end
if isempty(fs)
    error('Could not read the sampling rate; pass it with ''Fs''.');
end
fs = double(fs(1));
fprintf('Sampling rate  : %g Hz\n', fs);

%% ---------------- spike times ----------------
stInfo = h5info(resultFile, '/spiketimes');
tempNames = sortTempNames({stInfo.Datasets.Name});
nUnits = numel(tempNames);
fprintf('Units          : %d\n', nUnits);

unitTimes = cell(1, nUnits);
unitAmps  = cell(1, nUnits);
for i = 1:nUnits
    frames = double(h5read(resultFile, ['/spiketimes/' tempNames{i}]));
    [t, iSort] = sort(frames(:) / fs);
    unitTimes{i} = t;

    a = [];
    try
        raw = double(h5read(resultFile, ['/amplitudes/' tempNames{i}]));
        % stored 2 x nSpikes: first row is the fitted amplitude
        if size(raw, 1) == 2 && size(raw, 2) == numel(t)
            a = raw(1, :)';
        elseif size(raw, 2) == 2 && size(raw, 1) == numel(t)
            a = raw(:, 1);
        end
        if ~isempty(a)
            a = a(iSort);
        end
    catch
        % amplitudes are optional
    end
    unitAmps{i} = a;
end

fprintf('Total spikes   : %d\n', sum(cellfun(@numel, unitTimes)));

%% ---------------- duration ----------------
duration_s = opt.Duration;
% Prefer the exact raw frame count: /info/duration is whole seconds only.
% (h5info reads metadata only, so this works without the MaxWell plugin.)
if isempty(duration_s) && ~isempty(rawFile)
    try
        ri = h5info(rawFile, '/data_store/data0000/groups/routed/raw');
        duration_s = max(ri.Dataspace.Size) / fs;
    catch
        duration_s = [];
    end
end
if isempty(duration_s)
    try
        d = double(h5read(resultFile, '/info/duration'));
        d = d(1);
        % SpyKING CIRCUS writes this in samples on some versions
        if d > 1e5
            duration_s = d / fs;
        else
            duration_s = d;
        end
    catch
        duration_s = [];
    end
end
if isempty(duration_s)
    allT = vertcat(unitTimes{:});
    duration_s = ceil(max(allT));
    warning('Duration taken from the last spike (%.1f s).', duration_s);
end
duration_s = double(duration_s);
fprintf('Duration       : %.2f s\n', duration_s);
if opt.WholeSeconds && duration_s > floor(duration_s)
    nLate = sum(cellfun(@(t) sum(t > floor(duration_s)), unitTimes));
    duration_s = floor(duration_s);
    fprintf('Trimmed to     : %d s (MEA-NAP uses whole seconds; %d later spikes dropped)\n', ...
        duration_s, nLate);
end

analysisFs = double(opt.AnalysisFs);
if isempty(analysisFs)
    analysisFs = fs;
end
fprintf('MEA-NAP fs     : %g Hz (spike-matrix bin = %.3g ms; spike times keep full precision)\n', ...
    analysisFs, 1000 / analysisFs);

%% ---------------- peak channel per unit ----------------
peakChan = double(h5read(clustersFile, '/electrodes'));   % 0-based
peakChan = peakChan(:);
if numel(peakChan) ~= nUnits
    error(['/electrodes in %s has %d entries but the result file has %d units. ' ...
           'The two files are from different runs - use the matching merged/unmerged pair.'], ...
           shortName(clustersFile), numel(peakChan), nUnits);
end

%% ---------------- electrode mapping from the raw file ----------------
mapElectrode = [];
mapX = [];
mapY = [];
if ~isempty(rawFile)
    [mapElectrode, mapX, mapY] = readMaxwellMapping(rawFile);
end

%% ---------------- drop sparse units ----------------
keepUnit = cellfun(@numel, unitTimes) >= opt.MinSpikes;
if any(~keepUnit)
    fprintf('Dropping %d unit(s) with < %d spikes\n', sum(~keepUnit), opt.MinSpikes);
end
unitTimes = unitTimes(keepUnit);
unitAmps  = unitAmps(keepUnit);
peakChan  = peakChan(keepUnit);
unitIds   = find(keepUnit) - 1;        % original temp_K numbers
nUnits    = numel(unitTimes);

%% ---------------- assemble MEA-NAP nodes ----------------
mapRow = peakChan + 1;                 % 1-based row into the mapping table

switch lower(opt.NodeMode)
    case 'unit'
        nodeTimes  = unitTimes;
        nodeAmps   = unitAmps;
        nodeMapRow = mapRow;
        channels   = (1:nUnits)';
        nodeLabel  = unitIds;

    case 'electrode'
        [uRows, ~, grp] = unique(mapRow);
        nNodes = numel(uRows);
        nodeTimes = cell(1, nNodes);
        nodeAmps  = cell(1, nNodes);
        for g = 1:nNodes
            members = find(grp == g);
            t = [];
            a = [];
            aOk = true;
            for m = members(:)'
                t = [t; unitTimes{m}];                       %#ok<AGROW>
                if numel(unitAmps{m}) == numel(unitTimes{m})
                    a = [a; unitAmps{m}];                    %#ok<AGROW>
                else
                    aOk = false;
                end
            end
            [t, iSort] = sort(t);
            nodeTimes{g} = t;
            if aOk && numel(a) == numel(t)
                nodeAmps{g} = a(iSort);
            else
                nodeAmps{g} = [];
            end
        end
        nodeMapRow = uRows;
        if ~isempty(mapElectrode)
            channels = mapElectrode(uRows);
        else
            channels = uRows;
        end
        channels  = channels(:);
        nodeLabel = uRows - 1;

    otherwise
        error('NodeMode must be ''unit'' or ''electrode''.');
end

nNodes = numel(nodeTimes);
fprintf('MEA-NAP nodes  : %d (%s mode)\n', nNodes, lower(opt.NodeMode));

%% ---------------- coordinates ----------------
coords = nan(nNodes, 2);
if ~isempty(mapX)
    for i = 1:nNodes
        r = nodeMapRow(i);
        if r >= 1 && r <= numel(mapX)
            coords(i, :) = [mapX(r), mapY(r)];
        end
    end
end

if all(isnan(coords(:, 1)))
    warning('No electrode coordinates found; laying nodes out on a grid.');
    side = ceil(sqrt(nNodes));
    [gx, gy] = meshgrid(1:side, 1:side);
    coords = [gx(1:nNodes)', gy(1:nNodes)'];
else
    bad = isnan(coords(:, 1));
    if any(bad)
        warning('%d node(s) had no coordinate; placed at the centre.', sum(bad));
        coords(bad, 1) = mean(coords(~bad, 1));
        coords(bad, 2) = mean(coords(~bad, 2));
    end
end

% Normalised 0-1 copy: what MEA-NAP's custom layout expects before its own x8.
% One scale factor for both axes, so distances (and the aspect ratio of the
% electrode patches) are preserved; the longer side spans 0-1.
lo = min(coords, [], 1);
span = max(max(coords, [], 1) - lo);
if span > 0
    coordsNorm = (coords - lo) / span;
else
    coordsNorm = 0.5 * ones(size(coords));
end
if opt.FlipY
    % MEA-NAP plots with y pointing up; the chip's y counts down the array
    coordsNorm(:, 2) = max(coordsNorm(:, 2)) - coordsNorm(:, 2);
end
coordsMEANAP = coordsNorm * 8;         % = what getCoordsFromLayout returns

%% ---------------- stimulus TTL events (optional) ----------------
% /bits/0000 holds every change of the digital input word with its frame
% number. In this recording bit 5 (value 32) is a ~60 Hz display frame
% sync and bit 7 (value 128) drops briefly at each stimulus epoch change.
% Frame numbers are absolute, so subtract the first frame of the recording.
stimEvents = struct();
if ~isempty(rawFile)
    try
        b = h5read(rawFile, '/bits/0000');
        f0 = double(h5read(rawFile, '/data_store/data0000/groups/routed/frame_nos', 1, 1));
        fr = double(b.frameno(:));
        w  = double(b.bits(:));
        stimEvents.firstFrame  = f0;
        stimEvents.t_s         = (fr - f0) / fs;     % same clock as spikeTimes
        stimEvents.word        = w;
        bit5 = bitand(w, 32) > 0;
        bit7 = bitand(w, 128) > 0;
        stimEvents.frameSyncOn_s = stimEvents.t_s([bit5(1); diff(bit5) > 0]);
        stimEvents.marker_s      = stimEvents.t_s([~bit7(1); diff(bit7) < 0]);
        fprintf('Stim events    : %d TTL changes, %d epoch markers (first frame %d)\n', ...
            numel(fr), numel(stimEvents.marker_s), f0);
    catch ME
        fprintf('Stim events    : not read (%s)\n', ME.message);
    end
end

%% ---------------- build the MEA-NAP variables ----------------
spikeTimes     = cell(1, nNodes);
spikeWaveforms = cell(1, nNodes);
spikeAmps      = cell(1, nNodes);
thresholds     = cell(1, nNodes);

for i = 1:nNodes
    t = nodeTimes{i};
    inWindow = t >= 0 & t < duration_s;      % MEA-NAP bins 0 .. duration_s

    s = struct(); s.(method) = double(t(inWindow));
    spikeTimes{i} = s;

    w = struct(); w.(method) = [];           % waveforms not carried across
    spikeWaveforms{i} = w;

    a = struct();
    if numel(nodeAmps{i}) == numel(t)
        a.(method) = double(nodeAmps{i}(inWindow));
    else
        a.(method) = [];
    end
    spikeAmps{i} = a;

    th = struct(); th.(method) = NaN;
    thresholds{i} = th;
end

Params = struct();
Params.fs                  = analysisFs;   % what MEA-NAP bins the spike matrix at
Params.acquisitionFs       = fs;           % the chip's real sampling rate
Params.duration            = duration_s;
Params.unit                = 's';
Params.wnameList           = {method};
Params.SpikesMethod        = method;
Params.source              = 'SpyKING CIRCUS';
Params.sourceResultFile    = resultFile;
Params.sourceClustersFile  = clustersFile;
Params.sourceRawFile       = rawFile;
Params.nodeMode            = lower(opt.NodeMode);
Params.nUnitsSorted        = numel(unitIds);
Params.peakChannel         = nodeMapRow - 1;
Params.convertedOn         = datestr(now);

spikeDetectionResult = struct();
spikeDetectionResult.params = Params;
spikeDetectionResult.method = method;

%% ---------------- save ----------------
saveName = fullfile(outFolder, [base '_spikes.mat']);
save(saveName, 'spikeTimes', 'channels', 'spikeDetectionResult', ...
    'spikeWaveforms', 'spikeAmps', 'thresholds', 'coords', 'coordsMEANAP', ...
    'stimEvents', '-v7.3');

coordsCsv = fullfile(outFolder, 'meanap_custom_coords.csv');
T = table((1:nNodes)', channels(:), nodeLabel(:), nodeMapRow(:) - 1, ...
    coordsNorm(:, 1), coordsNorm(:, 2), coords(:, 1), coords(:, 2), ...
    cellfun(@numel, nodeTimes)', ...
    'VariableNames', {'node', 'channel', 'source_id', 'peak_channel', ...
                      'x_norm', 'y_norm', 'x_um', 'y_um', 'n_spikes'});
writetable(T, coordsCsv);

% A ready-made replacement for the Custom branch of getCoordsFromLayout.m
layoutM = fullfile(outFolder, 'meanap_custom_layout.m');
fid = fopen(layoutM, 'w');
fprintf(fid, '%% Paste this in place of the ''Custom'' branch of\n');
fprintf(fid, '%% MEA-NAP/Functions/getCoordsFromLayout.m\n');
fprintf(fid, '%% (the shipped branch references an undefined Params and errors).\n');
fprintf(fid, '%% Works for ONE layout. For several MaxWell recordings with different\n');
fprintf(fid, '%% electrode configurations, use the per-recording setUpSpreadSheet patch.\n');
fprintf(fid, 'elseif strcmp(channelLayout, ''Custom'')\n');
fprintf(fid, '    coordTable = readtable(''%s'');\n', coordsCsv);
fprintf(fid, '    channels = coordTable.channel'';\n');
fprintf(fid, '    coords = [coordTable.x_norm, coordTable.y_norm];\n');
fclose(fid);

batchCsv = fullfile(outFolder, 'MEANAP_batch.csv');
batch = table({base}, opt.AgeGroup, {char(opt.Group)}, {''}, ...
    'VariableNames', {'Recording filename', 'Age group', 'Group', 'Ground'});
writetable(batch, batchCsv);

fprintf('\nWrote:\n  %s\n  %s\n  %s\n  %s\n', ...
    saveName, coordsCsv, layoutM, batchCsv);

fprintf(['\nIn MEA-NAP:\n' ...
         '  Spike data folder        : %s\n' ...
         '  Spike method for analysis: %s\n' ...
         '  Detect spikes            : 0 (off)\n' ...
         '  Start analysis at step   : 2\n' ...
         '  Channel layout           : Custom (%d nodes)\n' ...
         '  Sampling frequency       : %g (MEA-NAP takes this from the spike file)\n' ...
         'Check the Age group and Group columns in MEANAP_batch.csv first.\n'], ...
         outFolder, method, nNodes, analysisFs);

out = struct();
out.saveName   = saveName;
out.coordsCsv  = coordsCsv;
out.layoutM    = layoutM;
out.batchCsv   = batchCsv;
out.nNodes     = nNodes;
out.fs         = analysisFs;
out.acquisitionFs = fs;
out.duration_s = duration_s;
out.channels   = channels;
out.coords     = coords;
out.coordsMEANAP = coordsMEANAP;
out.stimEvents = stimEvents;
out.table      = T;

end


%% ======================= helpers =======================

function f = pickFile(folder, patterns, label)
% Return the first pattern that matches, so merged files win when present.
for i = 1:numel(patterns)
    d = dir(fullfile(folder, patterns{i}));
    % '*.result.hdf5' also matches '*.result-merged.hdf5' on some platforms
    if ~isempty(d)
        keep = ~contains({d.name}, '-merged') | contains(patterns{i}, '-merged');
        d = d(keep);
    end
    if ~isempty(d)
        f = fullfile(d(1).folder, d(1).name);
        return
    end
end
error('No %s file found in %s (looked for: %s)', label, folder, strjoin(patterns, ', '));
end


function names = sortTempNames(names)
% Order temp_0, temp_1, temp_2 ... numerically, not alphabetically.
num = nan(1, numel(names));
for i = 1:numel(names)
    tok = regexp(names{i}, '(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        num(i) = str2double(tok{1});
    end
end
[~, order] = sort(num);
names = names(order);
end


function v = readFirstAvailable(file, candidates)
v = [];
for i = 1:numel(candidates)
    try
        v = h5read(file, candidates{i});
        if ~isempty(v)
            return
        end
    catch
        % try the next location
    end
end
end


function [electrode, x, y] = readMaxwellMapping(rawFile)
% Read the channel -> electrode / position table from a MaxWell file.
electrode = [];
x = [];
y = [];

candidates = { ...
    '/data_store/data0000/settings/mapping', ...
    '/recordings/rec0000/well000/settings/mapping', ...
    '/settings/mapping', ...
    '/mapping'};

for i = 1:numel(candidates)
    try
        m = h5read(rawFile, candidates{i});
    catch
        continue
    end
    if ~isstruct(m)
        continue
    end
    fn = fieldnames(m);
    xName = firstMatch(fn, {'x', 'X', 'x_um'});
    yName = firstMatch(fn, {'y', 'Y', 'y_um'});
    if isempty(xName) || isempty(yName)
        continue
    end
    x = double(m.(xName)(:));
    y = double(m.(yName)(:));

    eName = firstMatch(fn, {'electrode', 'Electrode', 'channel', 'Channel'});
    if isempty(eName)
        electrode = (1:numel(x))';
    else
        electrode = double(m.(eName)(:));
    end

    fprintf('Mapping        : %s (%d electrodes)\n', candidates{i}, numel(x));
    return
end

warning('No channel mapping found in the raw file; coordinates will be a grid.');
end


function name = firstMatch(available, wanted)
name = '';
for i = 1:numel(wanted)
    if ismember(wanted{i}, available)
        name = wanted{i};
        return
    end
end
end


function s = shortName(f)
[~, n, e] = fileparts(f);
s = [n e];
end
