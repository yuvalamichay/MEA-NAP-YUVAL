function peakZscores = computePeakZscoreForEachChannel(spikeTimes, allStimTimes, psth_window_s, baseline_window_s, psth_bin_width_s, psth_gaussian_width_ms, numChannels, numStimEvents, spikeMethod)
%computePeakZscoreForEachChannel Computes the peak z-score for each channel.
%   This function iterates over all channels and computes the peak z-score of
%   the post-stimulus time histogram (PSTH) for each one, using the provided
%   spike times and stimulus times.

peakZscores = zeros(numChannels, 1);

for chIdx = 1:numChannels
    
    if isempty(spikeTimes{chIdx}) || isempty(spikeTimes{chIdx}.(spikeMethod))
        peakZscores(chIdx) = 0;
        continue;
    end
    
    channelSpikeTimes = spikeTimes{chIdx}.(spikeMethod);
    
    peakZscores(chIdx) = computePeakZscore(channelSpikeTimes, allStimTimes, ...
        psth_window_s, baseline_window_s, psth_bin_width_s, psth_gaussian_width_ms, numStimEvents);
end

end
