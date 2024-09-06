threshs = [53.6 41.7 46.5 44.9 53.7 47.7 49.7 43.6 50.7 35.2 52.7 57.3 35.9]; %correct
% threshs = [53.6 41.7 46.5 44.9 47.7 49.7 43.6 50.7 35.2 52.7 57.3 35.9]; errors

% Change equality operator for posterior channels
for ii = 1:numel(Coordinates)
    antChan{ii} = Coordinates(ii).Electrodes(:,2) < threshs(ii);
end

%% Error power data structs were named differently 
for ii = 1:numel(PDPowers)
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,:,:,antChan{ii});
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,:,:,antChan{ii});
end

% for ii = 1:numel(PDPowers)
%     PDPowers(ii).goPowers = PDPowers(ii).goPowers(:,:,:,antChan{ii});
%     PDPowers(ii).ngPowers = PDPowers(ii).ngPowers(:,:,:,antChan{ii});
% end


%%
for ii = 1:numel(antChan1)
    antChan1{ii} = double(antChan1{ii});
end
% labels = [];
% for ii = 1:numel(antChan1)
%     labels = [labels; antChan1{ii}];
% end

%% Remove subjects with no electrodes remaining for rostral group
PDPowers([3 4 6]) = [];

PDPowers([3 4 7]) = [];
PDPhases([3 4 7]) = [];
PDMotorS([3 4 7]) = [];
NGErrors([3 4 6]) = [];
%% Remove FEF/pre-motor area channels
PDPowers(3).goPower(:,:,:,end) = [];
PDPowers(3).ngPower(:,:,:,end) = [];

PDPowers(4).goPower(:,:,:,end) = [];
PDPowers(4).ngPower(:,:,:,end) = [];

PDPowers(7).goPower(:,:,:,end-1:end) = [];
PDPowers(7).ngPower(:,:,:,end-1:end) = [];

PDPowers(9).goPower(:,:,:,end-1:end) = [];
PDPowers(9).ngPower(:,:,:,end-1:end) = [];

%% Include Only FEF channels
PDPowers(3).goPower = PDPowers(3).goPower(:,:,:,end);
PDPowers(3).ngPower = PDPowers(3).ngPower(:,:,:,end);

PDPowers(4).goPower = PDPowers(4).goPower(:,:,:,end);
PDPowers(4).ngPower = PDPowers(4).ngPower(:,:,:,end);

PDPowers(7).goPower = PDPowers(7).goPower(:,:,:,end-1:end);
PDPowers(7).ngPower = PDPowers(7).ngPower(:,:,:,end-1:end);

PDPowers(9).goPower = PDPowers(9).goPower(:,:,:,end-1:end);
PDPowers(9).ngPower = PDPowers(9).ngPower(:,:,:,end-1:end);

%% Remove subjects with no electrodes remaining for rostral group

PDPowers([8 13]) = [];

%% Remove FEF channels ERRORS
PDPowers(3).ngPower(:,:,:,end) = [];

PDPowers(4).ngPower(:,:,:,end) = [];

PDPowers(6).ngPower(:,:,:,end-1:end) = [];

PDPowers(8).ngPower(:,:,:,end-1:end) = [];

%% Remove FEF Coordinates
Coordinates(3).yDist(end) = [];

Coordinates(4).yDist(end) = [];

Coordinates(7).yDist(end-1:end) = [];

Coordinates(9).yDist(end-1:end) = [];

