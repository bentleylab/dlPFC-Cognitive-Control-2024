function [epoched_data,time,trials] = epochData(data,trials)

% StimulusTime is hard-coded to assume it was 10 kHz at recording and
% downsampled with factor of 10

%trials([trials.BlockType] == "S") = [];
%trials([trials.BlockType] == "R") = [];

% removeIndices = arrayfun(@(x) ~isempty(x.RT) && x.RT < 300, trials);
% trials(removeIndices) = [];

backwards = 2000; 
forwards = 2500; 
time = -backwards:forwards; % Set epoch length

epoched_data = NaN(length(time),numel(trials),min(size(data)));

for elecE = 1:min(size(data))
    % Get Go and Nogo trial epochs
    for iTrials = 1:length(trials)
            epoched_data(:,iTrials,elecE) = data(round(trials(iTrials).StimulusTime/10)-backwards...
                :round(trials(iTrials).StimulusTime/10)+forwards,elecE);    
    end
end
end