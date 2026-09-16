function varargout = axionStimEventsTool(action, varargin)
%AXIONSTIMEVENTSTOOL Shared logic for the 'axionStimEvents' stim detection method.
%
% Stimulation times are read straight from the Axion file's StimulationEvents
% (never from the voltage trace). A CSV with columns
%   (1) raw file name, (2) well, (3) stimulated electrode (channel id col*10+row)
% lists the stimulated electrodes of each well, one row per electrode. A well
% can have several rows (e.g. two electrodes stimulated alternately). Each
% StimulationEvent records the well and electrode(s) it was delivered to, so
% every listed electrode receives only the events delivered to it. The CSV is
% chosen on the General tab of the GUI (Params.axionStimCSV); the .raw files are
% looked for next to that CSV and in the MEA data folder. The same CSV/Params
% drive both the batch pipeline and the interactive stim detection app.
%
% Actions:
%   axionStimEventsTool('selectCSV', editFieldHandle, app)
%       File-picker callback for the 'Stim .raw CSV' button; writes the chosen
%       path into the edit field.
%   stimInfo = axionStimEventsTool('build', recName, channelNames, coords, Params, app)
%       Build stimInfo (one struct per channel) for one recording / well using
%       Params.axionStimCSV.
%   axionStimEventsTool('warnUnmatched', app)
%       Warn about CSV rows that never matched a built recording.
%   axionStimEventsTool('reset')
%       Forget the loaded CSV and cached events.

    persistent loadedCSVPath csvRows eventCache rowMatched

    if isempty(eventCache)
        eventCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end

    varargout = {};

    switch action
        case 'selectCSV'
            field = varargin{1};
            app = varargin{2};
            [csvName, csvPath] = uigetfile({'*.csv', 'CSV files (*.csv)'}, ...
                'Select stim .raw CSV (raw file name, well, electrode; one row per stimulated electrode)');
            if ~isequal(csvName, 0)
                field.Value = fullfile(csvPath, csvName);
            end
            if ~isempty(app) && isvalid(app)
                figure(app.UIFigure)   % return focus to the GUI after the dialog
            end

        case 'build'
            [recName, channelNames, coords, Params, app] = varargin{:};
            csvPath = '';
            if isfield(Params, 'axionStimCSV')
                csvPath = Params.axionStimCSV;
            end
            if isempty(csvPath)
                error('axionStimEvents:notConfigured', ...
                    ['No stim .raw CSV selected. Upload it with the "Stim .raw CSV" ' ...
                     'button on the General tab before using the axionStimEvents method.']);
            end
            if ~strcmp(csvPath, char(loadedCSVPath))
                csvRows = readAxionStimCSV(csvPath);
                loadedCSVPath = csvPath;
                eventCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
                rowMatched = false(numel(csvRows), 1);
            end
            if isempty(rowMatched) || numel(rowMatched) ~= numel(csvRows)
                rowMatched = false(numel(csvRows), 1);
            end
            rawFolders = candidateRawFolders(csvPath, Params);
            [stimInfo, eventCache, rowMatched] = buildForRecording( ...
                recName, channelNames, coords, Params, csvRows, rawFolders, ...
                eventCache, rowMatched, app);
            varargout = {stimInfo};

        case 'warnUnmatched'
            if ~isempty(rowMatched) && any(~rowMatched)
                warnUnmatchedRows(csvRows(~rowMatched), varargin{1});
            end

        case 'reset'
            loadedCSVPath = '';
            csvRows = [];
            eventCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
            rowMatched = [];

        otherwise
            error('axionStimEvents:badAction', 'Unknown axionStimEventsTool action "%s".', action);
    end
end


function folders = candidateRawFolders(csvPath, Params)
% Folders to search for the .raw files: next to the CSV, then the MEA data
% folder (where the Axion conversion co-locates .raw and .mat files).
    folders = {};
    csvFolder = fileparts(csvPath);
    if ~isempty(csvFolder)
        folders{end+1} = csvFolder;
    end
    if isfield(Params, 'rawData') && ~isempty(Params.rawData) && ischar(Params.rawData)
        folders{end+1} = Params.rawData;
    end
    folders = unique(folders, 'stable');
    if isempty(folders)
        folders = {pwd};
    end
end


function csvRows = readAxionStimCSV(csvFullPath)
% Read the CSV into a struct array with fields rawName / well / electrode, one
% element per stimulated electrode. Tolerates an optional header row, ignores
% blank rows, requires a numeric electrode (channel id), skips identical
% duplicate rows and flags malformed rows.

    if ~isfile(csvFullPath)
        error('axionStimEvents:csvNotFound', 'Stim .raw CSV not found: %s', csvFullPath);
    end

    raw = readcell(csvFullPath);
    if size(raw, 2) < 3
        error('axionStimEvents:badCSV', ...
            ['CSV "%s" must have at least 3 columns (raw file name, well, ' ...
             'electrode); found %d.'], csvFullPath, size(raw, 2));
    end

    csvRows = struct('rawName', {}, 'well', {}, 'electrode', {});
    seenKeys = {};

    for r = 1:size(raw, 1)
        if isBlankCell(raw{r, 1}) && isBlankCell(raw{r, 2}) && isBlankCell(raw{r, 3})
            continue   % skip fully empty rows
        end

        rawName = stripKnownExt(toChar(raw{r, 1}));
        well    = normalizeWell(toChar(raw{r, 2}));
        [elec, elecOk] = toElectrodeId(raw{r, 3});

        if ~elecOk
            if r == 1
                continue   % non-numeric electrode in row 1: treat as header
            end
            warning('axionStimEvents:malformedRow', ...
                'Skipping malformed CSV row %d (electrode "%s" is not a numeric channel id).', ...
                r, toChar(raw{r, 3}));
            continue
        end

        if isempty(rawName) || isempty(well)
            warning('axionStimEvents:malformedRow', ...
                'Skipping malformed CSV row %d (empty raw file name or well).', r);
            continue
        end

        % A well may list several electrodes; only an exact repeat is dropped.
        key = lower(sprintf('%s_%s_%d', rawName, well, elec));
        if any(strcmp(seenKeys, key))
            warning('axionStimEvents:duplicateRow', ...
                'Ignoring duplicate CSV row %d for %s (well %s, electrode %d).', ...
                r, rawName, well, elec);
            continue
        end

        csvRows(end+1) = struct('rawName', rawName, 'well', well, 'electrode', elec); %#ok<AGROW>
        seenKeys{end+1} = key; %#ok<AGROW>
    end

    if isempty(csvRows)
        error('axionStimEvents:emptyCSV', ...
            'No valid data rows found in CSV "%s".', csvFullPath);
    end
end


function [stimInfo, eventCache, rowMatched] = buildForRecording( ...
        recName, channelNames, coords, Params, csvRows, rawFolders, eventCache, rowMatched, app)
% Assemble stimInfo for one recording (one well). The well is identified by
% matching <rawName>_<well> from the CSV against the recording name. Every CSV
% row for that well names one stimulated electrode, which receives the
% StimulationEvents delivered to that electrode in that well.

    recKey = stripKnownExt(recName);

    matchIdx = [];
    for k = 1:numel(csvRows)
        if strcmpi([csvRows(k).rawName '_' csvRows(k).well], recKey)
            matchIdx(end+1) = k; %#ok<AGROW>
        end
    end

    if isempty(matchIdx)
        warning('axionStimEvents:noMatch', ...
            'No CSV row matches recording %s; assigning no stimulation to this well.', recName);
        stimInfo = buildStimInfo(channelNames, coords, [], {}, Params);
        return
    end

    rowMatched(matchIdx) = true;
    rows = csvRows(matchIdx);
    rawName = rows(1).rawName;
    well = rows(1).well;
    electrodes = [rows.electrode];

    missing = electrodes(~ismember(electrodes, channelNames));
    if ~isempty(missing)
        error('axionStimEvents:electrodeNotFound', ...
            ['Electrode(s) %s (CSV) are not channels in recording %s. ' ...
             'Available channels: %s.'], mat2str(missing), recName, mat2str(channelNames(:)'));
    end

    [events, eventCache] = getEvents(rawFolders, rawName, eventCache, app);
    stimTimes = cell(1, numel(electrodes));
    if isempty(events.times)
        warning('axionStimEvents:noEvents', ...
            'No StimulationEvents found in raw file "%s"; well %s will have no stimulation times.', ...
            rawName, recName);
        stimInfo = buildStimInfo(channelNames, coords, electrodes, stimTimes, Params);
        return
    end

    if isempty(events.pairTime)
        % The file does not say which electrode each event used. One listed
        % electrode can still take every event; several cannot be told apart.
        if numel(electrodes) > 1
            error('axionStimEvents:noElectrodeInfo', ...
                ['The stimulation events in raw file "%s" do not record which electrode ' ...
                 'they were delivered to, so they cannot be split between electrodes %s ' ...
                 'listed for well %s. List a single electrode for this well.'], ...
                rawName, mat2str(electrodes), well);
        end
        warning('axionStimEvents:noElectrodeInfo', ...
            ['The stimulation events in raw file "%s" do not record their electrode; ' ...
             'assigning all %d of them to electrode %d of well %s.'], ...
            rawName, numel(events.times), electrodes, well);
        stimTimes{1} = events.times;
        stimInfo = buildStimInfo(channelNames, coords, electrodes, stimTimes, Params);
        return
    end

    [wellRow, wellCol, wellOk] = parseWell(well);
    if ~wellOk
        error('axionStimEvents:badWell', ...
            'Cannot read well "%s" for raw file "%s"; expected a well name such as A1.', ...
            well, rawName);
    end

    inWell = events.pairWellRow == wellRow & events.pairWellCol == wellCol;
    wellChannels = unique(events.pairChannel(inWell))';
    if ~any(inWell)
        warning('axionStimEvents:noEventsInWell', ...
            ['Raw file "%s" has no stimulation events in well %s; stimulated wells in ' ...
             'this file: %s.'], rawName, well, strjoin(stimulatedWellNames(events), ', '));
    end

    for k = 1:numel(electrodes)
        stimTimes{k} = unique(events.pairTime(inWell & events.pairChannel == electrodes(k)));
        if isempty(stimTimes{k}) && any(inWell)
            warning('axionStimEvents:noEventsForElectrode', ...
                ['No stimulation events were delivered to electrode %d in well %s of raw file ' ...
                 '"%s"; electrodes stimulated in this well: %s.'], ...
                electrodes(k), well, rawName, mat2str(wellChannels));
        end
    end

    unlisted = setdiff(wellChannels, electrodes);
    if ~isempty(unlisted)
        warning('axionStimEvents:unlistedElectrodes', ...
            ['Electrode(s) %s in well %s of raw file "%s" were stimulated but are not in the ' ...
             'CSV; their stimulation times are ignored.'], mat2str(unlisted), well, rawName);
    end

    stimInfo = buildStimInfo(channelNames, coords, electrodes, stimTimes, Params);
end


function [events, eventCache] = getEvents(rawFolders, rawName, eventCache, app)
% Return the StimulationEvents of a raw file, read via AxisFile and cached per
% raw file (voltage is never loaded), as a struct with column vectors:
%   times            all event times (seconds, sorted)
%   pairTime, pairWellRow, pairWellCol, pairChannel
%                    one entry per (event, electrode) the event was delivered
%                    to; pairChannel is the channel id ElectrodeColumn*10 +
%                    ElectrodeRow, as written by rawConvertFunc
%   unlabelledTimes  times of events that name no electrode

    cacheKey = lower(rawName);
    if isKey(eventCache, cacheKey)
        events = eventCache(cacheKey);
        return
    end

    matchFile = '';
    for fi = 1:numel(rawFolders)
        rawList = dir(fullfile(rawFolders{fi}, '*.raw'));
        for f = 1:numel(rawList)
            if strcmpi(stripKnownExt(rawList(f).name), rawName)
                matchFile = fullfile(rawFolders{fi}, rawList(f).name);
                break
            end
        end
        if ~isempty(matchFile)
            break
        end
    end
    if isempty(matchFile)
        error('axionStimEvents:rawNotFound', ...
            'Raw file "%s.raw" referenced by the CSV was not found in: %s.', ...
            rawName, strjoin(rawFolders, '; '));
    end

    statusUpdate(app, sprintf('axionStimEvents: reading stimulation events from %s', rawName));
    fileData = AxisFile(matchFile);
    stimEvents = fileData.StimulationEvents;

    numEvents = numel(stimEvents);
    times = zeros(numEvents, 1);
    pairTime = cell(numEvents, 1);
    pairWellRow = cell(numEvents, 1);
    pairWellCol = cell(numEvents, 1);
    pairChannel = cell(numEvents, 1);
    labelled = true(numEvents, 1);

    for k = 1:numEvents
        times(k) = double(stimEvents(k).EventTime);
        % Electrodes is a ChannelMapping array, or a cell of them when the
        % stimulation block drives more than one channel group.
        mappings = stimEvents(k).Electrodes;
        if iscell(mappings)
            mappings = [mappings{:}];
        end
        if isempty(mappings)
            labelled(k) = false;
            continue
        end
        numMappings = numel(mappings);
        pairTime{k}    = repmat(times(k), numMappings, 1);
        pairWellRow{k} = double([mappings.WellRow])';
        pairWellCol{k} = double([mappings.WellColumn])';
        pairChannel{k} = double([mappings.ElectrodeColumn])' * 10 + double([mappings.ElectrodeRow])';
    end

    events = struct();
    events.times = sort(times);
    events.pairTime    = vertcat(zeros(0, 1), pairTime{:});
    events.pairWellRow = vertcat(zeros(0, 1), pairWellRow{:});
    events.pairWellCol = vertcat(zeros(0, 1), pairWellCol{:});
    events.pairChannel = vertcat(zeros(0, 1), pairChannel{:});
    events.unlabelledTimes = sort(times(~labelled));

    % Warned once per file here (events are cached); when no event names an
    % electrode, buildForRecording decides what to do instead.
    if any(labelled) && ~all(labelled)
        warning('axionStimEvents:unlabelledEvents', ...
            ['%d of %d stimulation events in raw file "%s" do not record their electrode ' ...
             'and are not assigned to any well.'], sum(~labelled), numEvents, rawName);
    end

    eventCache(cacheKey) = events;
end


function stimInfo = buildStimInfo(channelNames, coords, stimElectrodes, stimTimes, Params)
% Build the stimInfo cell (one struct per channel) in the exact format the
% existing methods produce. Electrode stimElectrodes(k) receives the times in
% stimTimes{k}; every other channel gets none. The blanking fields are
% populated so downstream artifact removal ignores [stimTime, stimTime +
% postStimWindowDur] (the GUI "post stim ignore duration").
%
% stimOrder ranks the stimulated electrodes of the well by their first
% stimulation time (1 = stimulated first; electrodes first stimulated at the
% same time share a rank); it is 0 for channels that were not stimulated.

    stimDur = Params.stimDuration;
    numChannels = length(channelNames);
    stimInfo = cell(numChannels, 1);

    firstTimes = inf(numel(stimElectrodes), 1);
    for k = 1:numel(stimElectrodes)
        if ~isempty(stimTimes{k})
            firstTimes(k) = min(stimTimes{k});
        end
    end
    stimOrder = zeros(numel(stimElectrodes), 1);
    fired = isfinite(firstTimes);
    if any(fired)
        [~, ~, stimOrder(fired)] = unique(firstTimes(fired));
    end

    for channel_idx = 1:numChannels
        k = find(stimElectrodes == channelNames(channel_idx), 1);
        if isempty(k)
            elecStimTimes = [];
            elecStimOrder = 0;
        else
            elecStimTimes = sort(stimTimes{k}(:));
            elecStimOrder = stimOrder(k);
        end

        stimStruct = struct();
        stimStruct.elecStimTimes = elecStimTimes;
        stimStruct.elecStimDur   = repmat(stimDur, length(elecStimTimes), 1);
        stimStruct.channelName   = channelNames(channel_idx);
        stimStruct.coords        = coords(channel_idx, :);
        stimStruct.stimOrder     = elecStimOrder;

        % Each stimulation time starts a blank; the blank duration is left at 0
        % so the ignored window equals [stimTime, stimTime + postStimWindowDur].
        stimStruct.blankStarts        = elecStimTimes;
        stimStruct.blankEnds          = elecStimTimes;
        stimStruct.nonStimBlankStarts = 0;
        stimStruct.nonStimBlankEnds   = 0;
        stimStruct.blankDurations     = zeros(length(elecStimTimes), 1);

        stimInfo{channel_idx} = stimStruct;
    end
end


function names = stimulatedWellNames(events)
    wells = unique([events.pairWellRow events.pairWellCol], 'rows');
    names = arrayfun(@(i) sprintf('%c%d', 'A' + wells(i, 1) - 1, wells(i, 2)), ...
        1:size(wells, 1), 'UniformOutput', false);
end


function warnUnmatchedRows(unmatchedRows, app)
    names = arrayfun(@(r) sprintf('%s_%s (electrode %d)', r.rawName, r.well, r.electrode), ...
        unmatchedRows, 'UniformOutput', false);
    msg = sprintf('axionStimEvents: %d CSV row(s) did not match any processed recording: %s', ...
        numel(unmatchedRows), strjoin(names, ', '));
    warning('axionStimEvents:unmatchedRows', '%s', msg);
    statusUpdate(app, msg);
end


function statusUpdate(app, msg)
    if ~isempty(app) && isvalid(app) && isprop(app, 'MEANAPStatusTextArea')
        app.MEANAPStatusTextArea.Value = [app.MEANAPStatusTextArea.Value; msg];
        drawnow
    end
end


function s = stripKnownExt(s)
% Remove any directory part and a trailing .raw/.mat extension only, preserving
% dots inside the name (e.g. "AK_HCNT24.4_DIV60").
    s = strtrim(s);
    sepIdx = find(s == '/' | s == '\', 1, 'last');
    if ~isempty(sepIdx)
        s = s(sepIdx+1:end);
    end
    if numel(s) >= 4 && (strcmpi(s(end-3:end), '.raw') || strcmpi(s(end-3:end), '.mat'))
        s = s(1:end-4);
    end
end


function w = normalizeWell(w)
    w = strtrim(w);
    w = regexprep(w, '^_+', '');   % drop leading underscores if present
end


function [wellRow, wellCol, ok] = parseWell(well)
% 'A1' -> row 1, column 1, matching the _A1 suffix written by rawConvertFunc and
% the WellRow / WellColumn of the Axion channel mappings.
    tok = regexp(well, '^([A-Za-z])0*(\d+)$', 'tokens', 'once');
    ok = ~isempty(tok);
    if ok
        wellRow = double(upper(tok{1})) - double('A') + 1;
        wellCol = str2double(tok{2});
    else
        wellRow = NaN;
        wellCol = NaN;
    end
end


function s = toChar(v)
    if ischar(v)
        s = strtrim(v);
    elseif isstring(v)
        s = strtrim(char(v));
    elseif isnumeric(v)
        if isempty(v) || (isscalar(v) && isnan(v))
            s = '';
        else
            s = strtrim(num2str(v));
        end
    elseif isa(v, 'missing')
        s = '';
    else
        s = strtrim(char(string(v)));
    end
end


function [id, ok] = toElectrodeId(v)
    ok = false;
    id = NaN;
    if isnumeric(v)
        if isscalar(v) && isfinite(v)
            id = double(v);
            ok = true;
        end
    else
        n = str2double(toChar(v));
        if ~isnan(n)
            id = n;
            ok = true;
        end
    end
    if ok
        id = round(id);
    end
end


function tf = isBlankCell(v)
    if isa(v, 'missing')
        tf = true;
    elseif isnumeric(v)
        tf = isempty(v) || all(isnan(v(:)));
    elseif ischar(v) || isstring(v)
        tf = strlength(strtrim(string(v))) == 0;
    else
        tf = false;
    end
end
