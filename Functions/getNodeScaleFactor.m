function nodeScaleF = getNodeScaleFactor(xc, yc, defaultScaleF, Params)
%GETNODESCALEFACTOR Size electrode markers from their actual spatial spacing.
% The old diameter remains an upper cap. Coincident sorted units do not
% contribute zero distances, and Params.nodeScaleFactor can override sizing.

validateattributes(defaultScaleF, {'numeric'}, {'scalar', 'finite', 'positive'});
if nargin >= 4 && isfield(Params, 'nodeScaleFactor') && ~isempty(Params.nodeScaleFactor)
    validateattributes(Params.nodeScaleFactor, {'numeric'}, {'scalar', 'finite', 'positive'});
    nodeScaleF = Params.nodeScaleFactor;
    return
end

x = double(xc(:));
y = double(yc(:));
nodeScaleF = defaultScaleF;
if numel(x) ~= numel(y) || numel(x) < 2 || any(~isfinite([x; y]))
    return
end

nearest = inf(numel(x), 1);
chunkSize = 1000;
for first = 1:chunkSize:numel(x)
    last = min(first + chunkSize - 1, numel(x));
    distances = hypot(x(first:last) - x.', y(first:last) - y.');
    % Ignore self distances and all coincident units before taking minima.
    distances(distances == 0) = Inf;
    nearest(first:last) = min(distances, [], 2);
end
valid = nearest(isfinite(nearest) & nearest > 0);
if ~isempty(valid)
    nodeScaleF = min(defaultScaleF, 0.8 * median(valid));
end
end
