function [normalized_data,idxs] = FindArtifacts(data,trials,chan)


for elec = 1:min(size(data))
    for triali = 1:length(trials)
            epoched_data(:,triali,elec) = data(round(trials(triali).StimulusTime/10):round(trials(triali).StimulusTime/10)+1000,elec);
    end
end



normalized_data = (epoched_data - mean(epoched_data,1))./std(epoched_data,0,1);

idxs = find(any(abs(normalized_data(:,:,chan)) >= 5.5));

% for triali = 1:length(idxs)
%     figure(idxs(triali))
%     plot(normalized_data(:,idxs(triali),chan))
%     ylim([-5 5])
% end


end
    