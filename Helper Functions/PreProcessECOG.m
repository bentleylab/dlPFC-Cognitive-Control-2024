function [INFO,ECOG_data,trials,FileInfo] = PreProcessECOG(subjectID,fname,baseDir)

% Load INFO file and data files
INFO = LoadInfoFile; % Load INFO file from Box Drive
subjectPath = [baseDir subjectID '/'];
taskfile = dir([subjectPath fname '*taskdata.mat']);
S = load([subjectPath taskfile.name],'trials','restperiods');
trials = S.trials;
load([subjectPath fname '.mat'],'Data','FileInfo')

% Preprocessing

FileNames = {INFO(strcmp({INFO.ID}, subjectID)).File.FileName}; % Collect all names of recordings


ns_idx = find(contains(FileNames,fname)); % Find recording of interest

% The below chunk of code identifies all good ECOG contacts
ECOG_bad_contacts = INFO(strcmp({INFO.ID}, subjectID)).File(ns_idx).BadECOGContacts;
ECOG_EIDs = INFO(strcmp({INFO.ID}, subjectID)).File(ns_idx).ECOGElectrodeIDs;
ECOG_good_IDs = ECOG_EIDs(~ismember(ECOG_EIDs, ECOG_EIDs(ECOG_bad_contacts)));
ECOG_good_idxs = ismember(FileInfo.chanlabels, ECOG_good_IDs);

% Extract data from good channels
Data = double(Data)';
ECOG_mat = Data(:,ECOG_good_idxs);

Fs = FileInfo.srate; % Extract sampling rate

% Bandpass Filter
[b,a] = butter(2, [0.5 500]/(Fs/2), 'bandpass'); % Setting butterworth filter; 2nd order 0.5 to 500 hz
ECOG_mat = filtfilt(b,a,ECOG_mat); % zero-phase filter

% Notch filter for line noise 60 Hz and harmonics
d1 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5, ...
               'DesignMethod','butter','SampleRate',Fs);
d2 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',119.5,'HalfPowerFrequency2',120.5, ...
               'DesignMethod','butter','SampleRate',Fs);
d3 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',179.5,'HalfPowerFrequency2',180.5, ...
               'DesignMethod','butter','SampleRate',Fs);
d4 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',239.5,'HalfPowerFrequency2',240.5, ...
               'DesignMethod','butter','SampleRate',Fs);

% zero-phase filters
ECOG_mat = filtfilt(d1,ECOG_mat);
ECOG_mat = filtfilt(d2,ECOG_mat);
ECOG_mat = filtfilt(d3,ECOG_mat);
ECOG_mat = filtfilt(d4,ECOG_mat);

% hard code downsampling by factor of 10 (usually gets data from 10 kHz to
% 1 kHz, but some subjects have higher than 10 kHz Fs)
ECOG_data = ECOG_mat(1:10:end,:);
end