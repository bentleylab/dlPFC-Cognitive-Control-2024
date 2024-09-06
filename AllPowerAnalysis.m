%%
baseDir = '/Users/anaskhan/Library/CloudStorage/Box-Box/Lab work/Anas Khan/DLPFC/'; % Set to where you stored the data
subjects = {
'DLPFC023',
'DLPFC025',
'DLPFC026',
'DLPFC030',
'DLPFC031',
'DLPFC032',
'DLPFC033',
'DLPFC034',
'DLPFC041',
'DLPFC044',
'DLPFC045',
'DLPFC052',
'DLPFC061',
'DLPFC074'};
files = {
'lfp2',
'task3',
'LFP3',
'task2',
'task1',
'Task_2',
'task2',
'task2',
'task_1',
'baseline',
'baseline',
'Baseline',
'baseline',
'baseline'};

%%
subnum = 1;
subjectID = subjects{subnum};

fname = files{subnum};
% Load and preprocess ECOG data
[INFO,ECOG_nostim,trials] = PreProcessECOG(subjectID,fname,baseDir);

%% Remove trials with artifacts by visual inspection
% Artifact identification
[~,idxs] = FindArtifacts(ECOG_nostim,trials,1);

% trials(idxs) = [];
%%
% Epoch Data
[ECOG_epoched,time,trials] = epochData(ECOG_nostim,trials);
% Adjacent bipolar rereferencing
ECOG_bipolar = ECOG_epoched(:,:,2:end) - ECOG_epoched(:,:,1:end-1);
%% Remove desired channels based on location in brain (e.g., not dlPFC)
ECOG_bipolar(:,:,1) = [];
%% ECOG
% Get powers
num_trials = size(ECOG_bipolar,2);
[CWT,frex,trimmedT] = MyCWT(ECOG_bipolar,num_trials,time);
powers = CWT.*conj(CWT);

% Normalize to pre-stim baseline for all trials

% Concatenate all trials to calculate mean and std for each freq and each electrode
BL_duration = length(find(trimmedT >= -500 & trimmedT <= -200));
BL_window = trimmedT >= -500 & trimmedT <= -200;
allBLs = powers(:,BL_window,~ismember(1:numel(trials),idxs),:);
allBLs = reshape(allBLs,[size(powers,1) BL_duration*size(allBLs,3) size(powers,4)]);

% Initialize means and stds matrices
meanBL = NaN(size(powers,1),size(powers,4));
sdBL = NaN(size(powers,1),size(powers,4));

% Loop over frequencies and electrodes
for elec = 1:size(powers,4)
    for fi = 1:size(powers,1)
        meanBL(fi,elec) = mean(allBLs(fi,:,elec));
        sdBL(fi,elec) = std(allBLs(fi,:,elec));
    end
end

clear allBLs
% Use means and stds to z-score raw data

% Initialize z-scored matrices
z_powers = NaN(size(powers,1),length(trimmedT),num_trials,size(powers,4));

for elec = 1:size(powers,4)
    for fi = 1:size(powers,1)
        z_powers(fi,:,:,elec) = (powers(fi,:,:,elec) - meanBL(fi,elec))./sdBL(fi,elec);
    end
end

%% This block of code is for the task-switching control analysis

% % Initialize an empty array to store indices of Go-NoGo trials in the GNG block
% % goNoGoIndices = [];
% 
% % Loop through the trials structure array
% for i = 2:length(trials)
%     % Check if the current trial and the previous trial are both Go trials in the GNG block
%     if strcmp(trials(i).Condition, 'NOGO') && ...
%        strcmp(trials(i).BlockType, 'GNG') && ...
%        strcmp(trials(i-1).Condition, 'GO') && ...
%        strcmp(trials(i-1).BlockType, 'GNG') && ...
%        trials(i).Repetition == trials(i-1).Repetition && ...
%        trials(i).ACC == 1 && ...
%        ~any(ismember(idxs,i))
%         % Add the index to the goGoIndices array
%         goNoGoIndices(end+1) = i;
%     end
% end
% 
% % Initialize an empty array to store indices of NoGo-NoGo trials in the GNG block
% noGoNoGoIndices = [];
% 
% % Loop through the trials structure array
% for i = 2:length(trials)
%     % Check if the current trial and the previous trial are both NoGo trials in the GNG block
%     if strcmp(trials(i).Condition, 'NOGO') && ...
%        strcmp(trials(i).BlockType, 'GNG') && ...
%        strcmp(trials(i-1).Condition, 'NOGO') && ...
%        strcmp(trials(i-1).BlockType, 'GNG') && ...
%        trials(i).Repetition == trials(i-1).Repetition && ...
%        trials(i).ACC == 1 && ...
%        ~any(ismember(idxs,i))
%         % Add the index to the noGoNoGoIndices array
%         noGoNoGoIndices(end+1) = i;
%     end
% end

% gNG_powers = z_powers(:,:,goNoGoIndices,:);
% ngNG_powers = z_powers(:,:,noGoNoGoIndices,:);
% 
% PDPowers(subnum).Subject = subjectID;
% PDPowers(subnum).Time = trimmedT;
% PDPowers(subnum).Frex = frex;
% PDPowers(subnum).gNGPowers = gNG_powers;
% PDPowers(subnum).ngNGPowers = ngNG_powers;
% PDPowers(subnum).Trials = trials;
% save('TaskSwitchingControlNG.mat','PDPowers','-append');

% AllPowers(subnum).Subject = subjectID;
% AllPowers(subnum).Time = trimmedT;
% AllPowers(subnum).Powers = z_powers;
% AllPowers(subnum).Trials = trials;
% AllPowers(subnum).Frex = frex;
% save('MVPAData.mat','AllPowers','-append');

clearvars -except files subjects baseDir Coordinates channums PDPowers AllPowers
%% This block of code gets different versions of errors
% incorrectTrials = [trials.ACC] == 0;
% etrials = trials(incorrectTrials);
% ez_powers = z_powers(:,:,incorrectTrials,:);
% 
% eNGtrials = [etrials.BlockType] == "GNG" & [etrials.Condition] == "NOGO";
% eNG_powers = ez_powers(:,:,eNGtrials,:);
% 
% eGOtrials = [etrials.BlockType] == "GNG" & [etrials.Condition] == "GO";
% eGO_powers = ez_powers(:,:,eGOtrials,:);
% 
% subji = length(PDPowers)+1;
% PDPowers(subji).Subject = subjectID;
% PDPowers(subji).Time = trimmedT;
% PDPowers(subji).Frex = frex;
% PDPowers(subji).goPowers = eGO_powers;
% PDPowers(subji).ngPowers = eNG_powers;
% 
% clearvars -except PDPowers

%% LCD vs HCD block trials (old names were motor an cognitive)

correctTrials = [trials.ACC] == 1;
trials = trials(correctTrials);
z_powers = z_powers(:,:,correctTrials,:);
powers = powers(:,:,correctTrials,:);

mGOtrials = [trials.BlockType] ~= "GNG";
mGO_powers = z_powers(:,:,mGOtrials,:);

GOtrials = [trials.BlockType] == "GNG" & [trials.Condition] == "GO";
GO_powers = z_powers(:,:,GOtrials,:);

NGtrials = [trials.BlockType] == "GNG" & [trials.Condition] == "NOGO";
NG_powers = z_powers(:,:,NGtrials,:);

cALLtrials = [trials.BlockType] == "GNG";

motor_BL = squeeze(mean(powers(:,126:626,mGOtrials,:),2));
cognitive_BL = squeeze(mean(powers(:,126:626,cALLtrials,:),2));

% motorRTs = [trials(mGOtrials).RT];
% cognitiveRTs = [trials(GOtrials).RT];
%

PDPowers(subnum).Subject = subjectID;
PDPowers(subnum).Time = trimmedT;
PDPowers(subnum).Frex = frex;
PDPowers(subnum).goPowers = GO_powers;
PDPowers(subnum).ngPowers = NG_powers;
PDPowers(subnum).motorBLPowers = motor_BL;
PDPowers(subnum).cogBLPowers = cognitive_BL;

% PDPowers(subnum).mgoPowers = mGO_powers;
% PDPowers(subnum).goRTs = cognitiveRTs;
% PDPowers(subnum).mgoRTs = motorRTs;

clearvars -except PDPowers
