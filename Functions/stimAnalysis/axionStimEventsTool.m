function varargout = axionStimEventsTool(action, varargin)
%AXIONSTIMEVENTSTOOL Shared logic for the 'axionStimEvents' stim detection method.
%
% Stimulation times are read straight from the Axion file's StimulationEvents
% (never from the voltage trace), and so are the stimulated electrodes: each
% StimulationEvent records the well and electrode(s) it was delivered to. No
% CSV is needed. The recording name <rawName>_<well> (e.g. "Plate1_DIV14_A1",
% the name Axion file conversion gives each per-well .mat) says which raw file
% and well a recording comes from. Every electrode of that well that received
% stimulation events is a stimulated electrode, with the times of its events.
%
% The .raw files are looked for in the MEA data folder (Params.rawData, a
% folder or a cell array of folders), where file conversion leaves them next to
% the per-well .mat files, and then in the optional extraRawFolder.
%
% Events that do not record their electrode (e.g. LED/optical stimulation,
% events without channel tags) cannot be assigned. If no event in a raw file
% records its electrode, a warning is given once per raw file and its wells get
% no stimulation. If only some do, the others are reported and ignored, and so
% are the events AxisFile could not link to their tags
% (UnlinkedStimulationEvents). Wells without stimulation events get no
% stimulation. The same logic drives both the batch pipeline and the
% interactive stim detection app.
%
% Actions:
%   stimInfo = axionStimEventsTool('build', recName, channelNames, coords, Params, app)
%   stimInfo = axionStimEventsTool('build', recName, channelNames, coords, Params, app, extraRawFolder)
%       Build stimInfo (one struct per channel) for one recording / well.
%       extraRawFolder is searched for the .raw after Params.rawData.
%   axionStimEventsTool('reset')
%       Forget the cached events.

    persistent eventCache

    if isempty(eventCache)
        eventCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end

    varargout = {};

    switch action
        case 'build'
            [recName, channelNames, coords, Params, app] = varargin{1:5};
            extraRawFolder = '';
            if numel(varargin) >= 6
                extraRawFolder = varargin{6};
            end
            rawFolders = candidateRawFolders(Params, extraRawFolder);
            [stimInfo, eventCache] = buildForRecording( ...
                recName, channelNames, coords, Params, rawFolders, eventCache, app);
            varargout = {stimInfo};

        case 'reset'
            eventCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

        otherwise
            error('axionStimEvents:badAction', 'Unknown axionStimEventsTool action "%s".', action);
    end
end


function folders = candidateRawFolders(Params, extraRawFolder)
% Folders to search for the .raw files: the MEA data folder(s) in
% Params.rawData (where the Axion conversion co-locates .raw and .mat files),
% then extraRawFolder. Empty entries are dropped and repeats removed.
    folders = {};
    if isfield(Params, 'rawData') && ~isempty(Params.rawData)
        rawData = Params.rawData;
        if isstring(rawData)
            rawData = cellstr(rawData);
        elseif ~iscell(rawData)
            rawData = {rawData};
        end
        folders = cellfun(@char, rawData(:)', 'UniformOutput', false);
    end
    if ~isempty(extraRawFolder)
        folders{end+1} = char(extraRawFolder);
    end
    folders = folders(~cellfun(@isempty, folders));
    folders = unique(folders, 'stable');
    if isempty(folders)
        folders = {pwd};
    end
end


function [stimInfo, eventCache] = buildForRecording( ...
        recName, channelNames, coords, Params, rawFolders, eventCache, app)
% Assemble stimInfo for one recording (one well). The raw file and well come
% from the recording name <rawName>_<well>; every electrode of that well named
% by a StimulationEvent receives the events delivered to it.

    recKey = stripKnownExt(recName);
    tok = regexp(recKey, '^(.+)_([A-Za-z]0*\d+)$', 'tokens', 'once');
    if isempty(tok)
        error('axionStimEvents:badRecordingName', ...
            ['Cannot tell the well of recording "%s". With the axionStimEvents method ' ...
             'the recording name must end in _<well> (for example _A1), as produced by ' ...
             'Axion file conversion.'], recName);
    end
    rawName = tok{1};
    well = tok{2};

    [wellRow, wellCol, wellOk] = parseWell(well);
    if ~wellOk
        error('axionStimEvents:badWell', ...
            'Cannot read well "%s" of recording "%s"; expected a well name such as A1.', ...
            well, recName);
    end

    [events, eventCache] = getEvents(rawFolders, rawName, eventCache, app);

    if isempty(events.times)
        statusUpdate(app, sprintf( ...
            'axionStimEvents: raw file %s has no stimulation events; %s: no stimulation', ...
            rawName, well));
        stimInfo = buildStimInfo(channelNames, coords, [], {}, Params);
        return
    end

    if isempty(events.pairTime)
        % No event names an electrode, so none can be assigned. Warned once per
        % raw file (the flag is kept with the cached events).
        if ~events.warnedNoElectrode
            warning('axionStimEvents:noElectrodeInfo', ...
                ['Raw file "%s" has %d stimulation events, but none record which electrode ' ...
                 'they were delivered to (for example LED stimulation, or events without ' ...
                 'channel tags). Well %s, and every other well of this file, is assigned ' ...
                 'no stimulation.'], rawName, numel(events.times), well);
            events.warnedNoElectrode = true;
            eventCache(lower(rawName)) = events;
        end
        statusUpdate(app, sprintf( ...
            'axionStimEvents: %s: no stimulation (events in %s do not record their electrode)', ...
            well, rawName));
        stimInfo = buildStimInfo(channelNames, coords, [], {}, Params);
        return
    end

    inWell = events.pairWellRow == wellRow & events.pairWellCol == wellCol;
    if ~any(inWell)
        % Unstimulated (e.g. control) wells are normal, so no warning.
        statusUpdate(app, sprintf( ...
            'axionStimEvents: %s: no stimulation (stimulated wells in %s: %s)', ...
            well, rawName, strjoin(stimulatedWellNames(events), ', ')));
        stimInfo = buildStimInfo(channelNames, coords, [], {}, Params);
        return
    end

    electrodes = unique(events.pairChannel(inWell))';
    stimTimes = cell(1, numel(electrodes));
    for k = 1:numel(electrodes)
        stimTimes{k} = unique(events.pairTime(inWell & events.pairChannel == electrodes(k)));
    end

    notRecorded = ~ismember(electrodes, channelNames);
    if any(notRecorded)
        warning('axionStimEvents:electrodeNotInRecording', ...
            ['Electrode(s) %s of well %s were stimulated according to raw file "%s" but ' ...
             'are not channels of recording %s; their stimulation times are ignored.'], ...
            mat2str(electrodes(notRecorded)), well, rawName, recName);
        electrodes = electrodes(~notRecorded);
        stimTimes = stimTimes(~notRecorded);
    end

    if isempty(electrodes)
        statusUpdate(app, sprintf('axionStimEvents: %s: no stimulation', well));
    else
        statusUpdate(app, sprintf('axionStimEvents: %s: stimulated electrodes %s (%s events)', ...
            well, strjoin(arrayfun(@num2str, electrodes, 'UniformOutput', false), ', '), ...
            joinWithAnd(cellfun(@numel, stimTimes))));
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
%   warnedNoElectrode
%                    set by buildForRecording once it has warned that no
%                    event in this file names an electrode

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
            ['Raw file "%s.raw" was not found in: %s. The axionStimEvents method reads ' ...
             'stimulation events from the .raw, so keep it in the MEA data folder next ' ...
             'to the converted .mat files.'], rawName, strjoin(rawFolders, '; '));
    end

    statusUpdate(app, sprintf('axionStimEvents: reading stimulation events from %s', rawName));
    fileData = AxisFile(matchFile);
    if ~isempty(fileData.UnlinkedStimulationEvents)
        warning('axionStimEvents:unlinkedEvents', ...
            ['%d stimulation events in raw file "%s" lack the tags that say where they ' ...
             'were delivered and are ignored.'], ...
            numel(fileData.UnlinkedStimulationEvents), rawName);
    end
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
    events.warnedNoElectrode = false;

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


function s = joinWithAnd(values)
% [120 120 90] -> '120, 120 and 90'
    parts = arrayfun(@num2str, values, 'UniformOutput', false);
    if numel(parts) > 1
        s = [strjoin(parts(1:end-1), ', ') ' and ' parts{end}];
    else
        s = strjoin(parts, '');
    end
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
