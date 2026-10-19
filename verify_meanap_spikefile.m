function ok = verify_meanap_spikefile(matFile)
% VERIFY_MEANAP_SPIKEFILE  Check a converted *_spikes.mat before running MEA-NAP.
%
% Checks the converted file fields consumed by MEA-NAP,
% in the way MEA-NAP reads it, and reports the memory the pipeline will try to
% allocate for the spike matrix.
%
% Usage:
%   verify_meanap_spikefile                       % finds the *_spikes.mat below pwd
%   verify_meanap_spikefile('/path/to/X_spikes.mat')
%
% Returns true if no FAIL-level problem was found.

if nargin < 1 || isempty(matFile)
    d = dir(fullfile(pwd, '**', '*_spikes.mat'));
    if isempty(d)
        error('No *_spikes.mat found below %s', pwd);
    end
    matFile = fullfile(d(1).folder, d(1).name);
end

fprintf('File: %s\n', matFile);
info = dir(matFile);
fprintf('Size: %.1f MB\n\n', info.bytes / 1e6);

S = load(matFile);
nFail = 0;
nWarn = 0;

%% ---- required variables -------------------------------------------------
required = {'spikeTimes', 'channels', 'spikeDetectionResult'};
for i = 1:numel(required)
    if ~isfield(S, required{i})
        nFail = nFail + 1;
        fprintf('FAIL  missing variable: %s\n', required{i});
    end
end
if nFail > 0
    fprintf('\nCannot continue without the required variables.\n');
    ok = false;
    return
end

optional = {'spikeWaveforms', 'spikeAmps', 'thresholds', 'coords', ...
            'coordsMEANAP', 'stimEvents'};
present = optional(isfield(S, optional));
fprintf('Variables present : %s\n', strjoin([required present], ', '));

%% ---- method name --------------------------------------------------------
spikeTimes = S.spikeTimes;
nNodes = numel(spikeTimes);
if ~iscell(spikeTimes)
    nFail = nFail + 1;
    fprintf('FAIL  spikeTimes is %s, must be a cell array\n', class(spikeTimes));
    ok = false; return
end

methods = fieldnames(spikeTimes{1});
fprintf('Nodes             : %d\n', nNodes);
fprintf('Method field(s)   : %s\n', strjoin(methods', ', '));
if numel(methods) ~= 1
    nWarn = nWarn + 1;
    fprintf('WARN  more than one method in the file; set MEA-NAP to the one you want\n');
end
method = methods{1};

% every node must carry the same field, or formatSpikeTimes errors
missingField = 0;
for i = 1:nNodes
    if ~isfield(spikeTimes{i}, method), missingField = missingField + 1; end
end
if missingField > 0
    nFail = nFail + 1;
    fprintf('FAIL  %d node(s) lack the field "%s"\n', missingField, method);
end

%% ---- fs and duration ----------------------------------------------------
sdr = S.spikeDetectionResult;
if isfield(sdr, 'params'), prm = sdr.params;
elseif isfield(sdr, 'Params'), prm = sdr.Params;
else
    nFail = nFail + 1;
    fprintf('FAIL  spikeDetectionResult has neither .params nor .Params\n');
    ok = false; return
end

if ~isfield(prm, 'fs')
    nFail = nFail + 1; fprintf('FAIL  params.fs missing\n'); fs = NaN;
else
    fs = double(prm.fs);
end
if ~isfield(prm, 'duration')
    nFail = nFail + 1; fprintf('FAIL  params.duration missing\n'); dur = NaN;
else
    dur = double(prm.duration);
end

fprintf('params.fs         : %g Hz\n', fs);
if isfield(prm, 'acquisitionFs')
    fprintf('acquisition fs    : %g Hz (the chip''s real rate)\n', double(prm.acquisitionFs));
end
fprintf('params.duration   : %g s\n', dur);

% MEA-NAP floors the duration when it bins
if isfinite(dur) && dur > floor(dur)
    nWarn = nWarn + 1;
    fprintf(['WARN  duration is not a whole number. MEA-NAP floors it to %d s when\n' ...
             '      building the spike matrix, so spikes after that point have no bin.\n'], ...
             floor(dur));
end

%% ---- spike matrix memory ------------------------------------------------
if isfinite(fs) && isfinite(dur)
    nSamples = floor(dur) * fs;
    gb = nSamples * nNodes * 8 / 2^30;
    fprintf('Spike matrix      : %d samples x %d nodes = %.2f GB when made full\n', ...
        round(nSamples), nNodes, gb);
    if gb > 2
        nWarn = nWarn + 1;
        fprintf(['WARN  MEA-NAP calls full() on this matrix in several places (raster,\n' ...
                 '      effective rank, NMF). %.1f GB is likely to run out of memory.\n' ...
                 '      Re-convert with a lower AnalysisFs (1000 Hz gives %.2f GB).\n'], ...
                 gb, nSamples / fs * 1000 * nNodes * 8 / 2^30);
    end
end

%% ---- spike times --------------------------------------------------------
nSpk = zeros(nNodes, 1);
badSort = 0; badNeg = 0; badLate = 0; badType = 0; badShape = 0; nonFinite = 0;
tMax = -Inf;

for i = 1:nNodes
    t = spikeTimes{i}.(method);
    if ~isa(t, 'double'), badType = badType + 1; end
    if ~isvector(t) && ~isempty(t), badShape = badShape + 1; end
    t = double(t(:));
    nSpk(i) = numel(t);
    if isempty(t), continue; end
    if any(~isfinite(t)), nonFinite = nonFinite + 1; end
    if any(diff(t) < 0), badSort = badSort + 1; end
    if any(t < 0), badNeg = badNeg + 1; end
    if isfinite(dur) && any(t >= floor(dur)), badLate = badLate + 1; end
    tMax = max(tMax, max(t));
end

fprintf('\nSpikes            : %d total\n', sum(nSpk));
fprintf('Per node          : min %d, median %d, max %d\n', ...
    min(nSpk), round(median(nSpk)), max(nSpk));
if isfinite(dur) && dur > 0
    fprintf('Firing rate       : median %.2f Hz, max %.2f Hz\n', ...
        median(nSpk)/dur, max(nSpk)/dur);
end
fprintf('Last spike at     : %.3f s\n', tMax);

if badType > 0
    nWarn = nWarn + 1;
    fprintf('WARN  %d node(s) store spike times as something other than double\n', badType);
end
if badShape > 0
    nFail = nFail + 1;
    fprintf('FAIL  %d node(s) store a non-vector in the method field\n', badShape);
end
if nonFinite > 0
    nFail = nFail + 1;
    fprintf('FAIL  %d node(s) contain NaN or Inf spike times\n', nonFinite);
end
if badSort > 0
    nWarn = nWarn + 1;
    fprintf('WARN  %d node(s) have unsorted spike times\n', badSort);
end
if badNeg > 0
    nFail = nFail + 1;
    fprintf('FAIL  %d node(s) have negative spike times\n', badNeg);
end
if badLate > 0
    nFail = nFail + 1;
    fprintf(['FAIL  %d node(s) have spikes at or beyond floor(duration) = %d s.\n' ...
             '      These fall outside the spike matrix MEA-NAP builds.\n'], ...
             badLate, floor(dur));
end

% a units sanity check: seconds, not samples or milliseconds
if isfinite(dur) && tMax > dur * 1.5
    nFail = nFail + 1;
    fprintf(['FAIL  the last spike (%.1f) is well beyond the duration (%.1f).\n' ...
             '      Spike times are probably in samples or ms, not seconds.\n'], tMax, dur);
end
if sum(nSpk) > 0 && tMax < dur * 0.1
    nWarn = nWarn + 1;
    fprintf('WARN  the last spike is very early relative to the duration. Check the units.\n');
end

%% ---- channels -----------------------------------------------------------
ch = S.channels;
fprintf('\nchannels          : %d entries, class %s\n', numel(ch), class(ch));
if numel(ch) ~= nNodes
    nFail = nFail + 1;
    fprintf('FAIL  channels has %d entries but there are %d nodes\n', numel(ch), nNodes);
end
if numel(unique(ch)) ~= numel(ch)
    nWarn = nWarn + 1;
    fprintf('WARN  channel labels are not unique (%d distinct of %d)\n', ...
        numel(unique(ch)), numel(ch));
end

%% ---- coordinates --------------------------------------------------------
if isfield(S, 'coords')
    c = S.coords;
    fprintf('coords            : %d x %d\n', size(c, 1), size(c, 2));
    if size(c, 1) ~= nNodes
        nFail = nFail + 1;
        fprintf('FAIL  coords has %d rows but there are %d nodes\n', size(c, 1), nNodes);
    end
    if any(~isfinite(c(:)))
        nFail = nFail + 1;
        fprintf('FAIL  coords contains NaN or Inf\n');
    else
        fprintf('  x range %.0f to %.0f, y range %.0f to %.0f\n', ...
            min(c(:,1)), max(c(:,1)), min(c(:,2)), max(c(:,2)));
        [~, iu] = unique(c, 'rows');
        nDup = nNodes - numel(iu);
        if nDup > 0
            fprintf('  %d node(s) share a position with another node\n', nDup);
        end
    end
else
    fprintf('coords            : not saved (MEA-NAP takes coordinates from\n');
    fprintf('                    getCoordsFromLayout, so this is informational only)\n');
end

%% ---- verdict ------------------------------------------------------------
fprintf('\n----------------------------------------\n');
if nFail == 0 && nWarn == 0
    fprintf('PASS - no problems found.\n');
elseif nFail == 0
    fprintf('PASS with %d warning(s).\n', nWarn);
else
    fprintf('%d FAILURE(S) and %d warning(s) - fix before running MEA-NAP.\n', nFail, nWarn);
end

fprintf(['\nMEA-NAP settings for this file:\n' ...
         '  Spike data folder         : %s\n' ...
         '  Spike method for analysis : %s\n' ...
         '  Detect spikes             : 0 (off)\n' ...
         '  Start analysis at step    : 2\n' ...
         '  Channel layout            : Custom (%d nodes)\n' ...
         '  Sampling frequency        : %g\n'], ...
         fileparts(matFile), method, nNodes, fs);

ok = (nFail == 0);
end
