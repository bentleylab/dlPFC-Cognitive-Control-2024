%% Broadband Go vs No-go 
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,3),size(PDPowers(subject).ngPower,3),true);
        bootGO(permi,:,:,:) = squeeze(mean(PDPowers(subject).goPower(:,:,idxs,:),3));
    end
    
        subjectBoots(subject).Boots = squeeze(mean(bootGO)); % average over 1000 iterations     

        clear bootGO
end

trueDiffs = NaN(length(PDPowers(1).Frex),length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,:,subject) = mean(squeeze(mean(PDPowers(subject).ngPower,3)) - subjectBoots(subject).Boots,3);
end

trueDiff = mean(trueDiffs,3);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1),size(trueDiff,2));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),size(trueDiff,2),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,3),size(PDPowers(subject).ngPower,3),true);

        % Combine trials
        allTrials = cat(3,PDPowers(subject).goPower(:,:,idxs,:),PDPowers(subject).ngPower);
    
        % Making mapping vector
        mapping = [ones(size(PDPowers(subject).ngPower,3),1); ones(size(PDPowers(subject).ngPower,3),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,:,subject) = mean(squeeze(mean( allTrials(:,:,new_map==2,:) ,3)) - squeeze(mean( allTrials(:,:,new_map==1,:) ,3)),3);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:,:) = mean(diff,3);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs));
std_h0 = squeeze(std(perm_diffs));

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:,:));
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

% Go
GoPower = [];

for ii = 1:numel(PDPowers)
    GoPower(:,:,ii) = mean(squeeze(mean(PDPowers(ii).goPower,3)),3);
end

% Go EMG
GoEMG = [];

for ii = 1:numel(PDEMG)
    GoEMG(:,ii) = mean(PDEMG(ii).correctGO,2);
end
GOEMG = mean(GoEMG,2);


% NoGo
NGPower = [];

for ii = 1:numel(PDPowers)
    NGPower(:,:,ii) = mean(squeeze(mean(PDPowers(ii).ngPower,3)),3);
end

% NoGo EMG
NGEMG = [];

for ii = 1:numel(PDEMG)
    NGEMG(:,ii) = mean(PDEMG(ii).correctNG,2);
end
NGEMG = mean(NGEMG,2);

[~,fidx] = arrayfun(@(x) min(abs(x-PDPowers(1).Frex)), [2 4 10 40 194]);
frexticks = PDPowers(1).Frex(fidx);
clear fidx

fig = figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
contourf(PDPowers(1).Time,PDPowers(1).Frex,mean(GoPower,3),50,'LineColor','none');
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 10 40 194])
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-0.5 0.5])
set(cb,'yTick',[-0.5 0 0.5])
cb.Label.String = ['\bf' 'Power (z)' '\rm'];
yyaxis right
plot(PDEMG(1).Time,GOEMG,'--k','LineWidth',2);
ax = gca;
ax.YAxis(2).Visible = 'off';
yl = ylim;


fig = figure(2);

left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
contourf(PDPowers(1).Time,PDPowers(1).Frex,mean(NGPower,3),50,'LineColor','none');
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 10 40 194])
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-0.5 0.5])
set(cb,'yTick',[-0.5 0 0.5])
cb.Label.String = ['\bf' 'Power (z)' '\rm'];
yyaxis right
plot(PDEMG(1).Time,NGEMG,'--k','LineWidth',2);
ax = gca;
ax.YAxis(2).Visible = 'off';
ylim(yl)

fig = figure(3);

left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
contourf(PDPowers(1).Time,PDPowers(1).Frex,trueDiff,50,'LineColor','none');
hold on
contour(PDPowers(1).Time,PDPowers(1).Frex,logical(zmap),1,'linecolor','k','linewidth',1.0);
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 10 40 194])
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-0.2 0.2])
set(cb,'yTick',[-0.2 0 0.2])
cb.Label.String = ['\bf' 'Power (z)' '\rm'];
yyaxis right
ax = gca;
ax.YAxis(2).Visible = 'off';
ylim(yl)


%% Theta Cognitive Control 
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
end

% Get theta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = squeeze(mean(PDPowers(ii).ngPower(1:17,:,:,:)));
    PDPowers(ii).goPower = squeeze(mean(PDPowers(ii).goPower(1:17,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(PDPowers(subject).ngPower,2),true);
        bootGO(permi,:,:) = squeeze(mean(PDPowers(subject).goPower(:,idxs,:),2));
    end
    if numel(size(bootGO)) == 3
        subjectBoots(subject).Boots = squeeze(mean(bootGO)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(bootGO))'; % average over 1000 iterations
    end

    clear bootGO
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(squeeze(mean(PDPowers(subject).ngPower,2)) - subjectBoots(subject).Boots,2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(PDPowers(subject).ngPower,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).goPower(:,idxs,:),PDPowers(subject).ngPower);
    
        % Making mapping vector
        mapping = [ones(size(PDPowers(subject).ngPower,2),1); ones(size(PDPowers(subject).ngPower,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    thetaPowerNG(:,ii) = mean(mean(PDPowers(ii).ngPower,2),3);
    thetaPowerGO(:,ii) = mean(mean(PDPowers(ii).goPower,2),3);
end

plot(t,thetaPowerGO(:,end),'b','LineWidth',2)
hold on
plot(t,thetaPowerNG(:,end),'r','LineWidth',2)

thetaPowerGO = thetaPowerGO(:,[2:3,5:end]);
thetaPowerNG = thetaPowerNG(:,[2:3,5:end]);

mthetaPowerGO = mean(thetaPowerGO,2)';
mthetaPowerNG = mean(thetaPowerNG,2)';

semthetaPowerGO = std(thetaPowerGO,0,2)/sqrt(size(thetaPowerGO,2));
semthetaPowerNG = std(thetaPowerNG,0,2)/sqrt(size(thetaPowerNG,2));

semthetaPowerGO = semthetaPowerGO';
semthetaPowerNG = semthetaPowerNG';

%
t = PDPowers(1).Time;
sigClusters = t(zmap==1);


%% Control for task switching theta No-go trials only

% Get theta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngNGPowers = squeeze(mean(PDPowers(ii).ngNGPowers(1:17,:,:,:)));
    PDPowers(ii).gNGPowers = squeeze(mean(PDPowers(ii).gNGPowers(1:17,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).gNGPowers,2),size(PDPowers(subject).ngNGPowers,2),true);
        bootGO(permi,:,:) = squeeze(mean(PDPowers(subject).gNGPowers(:,idxs,:),2));
    end
    if numel(size(bootGO)) == 3
        subjectBoots(subject).Boots = squeeze(mean(bootGO)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(bootGO))'; % average over 1000 iterations
    end

    clear bootGO
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(squeeze(mean(PDPowers(subject).ngNGPowers,2)) - subjectBoots(subject).Boots,2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).gNGPowers,2),size(PDPowers(subject).ngNGPowers,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).gNGPowers(:,idxs,:),PDPowers(subject).ngNGPowers);
    
        % Making mapping vector
        mapping = [ones(size(PDPowers(subject).ngNGPowers,2),1); ones(size(PDPowers(subject).ngNGPowers,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    thetaPowerNG(:,ii) = mean(mean(PDPowers(ii).ngNGPowers,2),3);
    thetaPowerGO(:,ii) = mean(mean(PDPowers(ii).gNGPowers,2),3);
end


mthetaPowerGO = mean(thetaPowerGO,2);
mthetaPowerNG = mean(thetaPowerNG,2);

semthetaPowerGO = std(thetaPowerGO,0,2)/sqrt(size(thetaPowerGO,2));
semthetaPowerNG = std(thetaPowerNG,0,2)/sqrt(size(thetaPowerNG,2));

semthetaPowerGO = semthetaPowerGO';
semthetaPowerNG = semthetaPowerNG';

[hl,~] = boundedline(PDPowers(1).Time,mthetaPowerNG,semthetaPowerNG,'b',PDPowers(1).Time,mthetaPowerGO,semthetaPowerGO,'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-202 1000])
ylim([-0.2 0.7])
xlabel('Time (ms)')
ylabel('LF Power (z)')

%
t = PDPowers(1).Time;
sigClusters = t(zmap==1);
%% Beta Cognitive Control 
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
end

% Downnsample to 40 Hz
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,1:25:end,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,1:25:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(1:25:end);
end

% Get beta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = squeeze(mean(PDPowers(ii).ngPower(27:41,:,:,:)));
    PDPowers(ii).goPower = squeeze(mean(PDPowers(ii).goPower(27:41,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(PDPowers(subject).ngPower,2),true);
        bootGO(permi,:,:) = squeeze(mean(PDPowers(subject).goPower(:,idxs,:),2));
    end
    if numel(size(bootGO)) == 3
        subjectBoots(subject).Boots = squeeze(mean(bootGO)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(bootGO))'; % average over 1000 iterations
    end

    clear bootGO
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(squeeze(mean(PDPowers(subject).ngPower,2)) - subjectBoots(subject).Boots,2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(PDPowers(subject).ngPower,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).goPower(:,idxs,:),PDPowers(subject).ngPower);
    
        % Making mapping vector
        mapping = [ones(size(PDPowers(subject).ngPower,2),1); ones(size(PDPowers(subject).ngPower,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    betaPowerNG(:,ii) = mean(mean(PDPowers(ii).ngPower,2),3);
    betaPowerGO(:,ii) = mean(mean(PDPowers(ii).goPower,2),3);
end

plot(t,betaPowerGO(:,end),'b','LineWidth',2)
hold on
plot(t,betaPowerNG(:,end),'r','LineWidth',2)


mbetaPowerGO = mean(betaPowerGO,2)';
mbetaPowerNG = mean(betaPowerNG,2)';

sembetaPowerGO = std(betaPowerGO,0,2)/sqrt(size(betaPowerGO,2));
sembetaPowerNG = std(betaPowerNG,0,2)/sqrt(size(betaPowerNG,2));

sembetaPowerGO = sembetaPowerGO';
sembetaPowerNG = sembetaPowerNG';

%% Theta Cognitive Control 
for ii = 1:numel(PDPowers)
    PDPowers(ii).goPowers = PDPowers(ii).goPowers(:,426:1626,:,:);
    PDPowers(ii).mgoPowers = PDPowers(ii).mgoPowers(:,426:1626,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
end

% Get theta
for ii = 1:numel(PDPowers)
    PDPowers(ii).goPowers = squeeze(mean(PDPowers(ii).goPowers(1:17,:,:,:)));
    PDPowers(ii).mgoPowers = squeeze(mean(PDPowers(ii).mgoPowers(1:17,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPowers,2),size(PDPowers(subject).mgoPowers,2),true);
        bootGO(permi,:,:) = squeeze(mean(PDPowers(subject).goPowers(:,idxs,:),2));
    end
    if numel(size(bootGO)) == 3
        subjectBoots(subject).Boots = squeeze(mean(bootGO)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(bootGO))'; % average over 1000 iterations
    end

    clear bootGO
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(subjectBoots(subject).Boots - squeeze(mean(PDPowers(subject).mgoPowers,2)),2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPowers,2),size(PDPowers(subject).mgoPowers,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).mgoPowers,PDPowers(subject).goPowers(:,idxs,:));
    
        % Making mapping vector
        mapping = [ones(size(PDPowers(subject).mgoPowers,2),1); ones(size(PDPowers(subject).mgoPowers,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    thetaPowerMG(:,ii) = mean(mean(PDPowers(ii).mgoPowers,2),3);
    thetaPowerGO(:,ii) = mean(mean(PDPowers(ii).goPowers,2),3);
end

mthetaPowerGO = mean(thetaPowerGO,2)';
mthetaPowerMG = mean(thetaPowerMG,2)';

semthetaPowerGO = std(thetaPowerGO,0,2)/sqrt(size(thetaPowerGO,2));
semthetaPowerMG = std(thetaPowerMG,0,2)/sqrt(size(thetaPowerMG,2));

semthetaPowerGO = semthetaPowerGO';
semthetaPowerMG = semthetaPowerMG';

%
t = PDPowers(1).Time;

plot(t,mthetaPowerGO,'r','LineWidth',2)
patch([t fliplr(t)],[mthetaPowerGO-semthetaPowerGO fliplr(mthetaPowerGO+semthetaPowerGO)],'r','FaceAlpha',0.35,'EdgeColor','none')
hold on
plot(t,mthetaPowerMG,'b','LineWidth',2)
patch([t fliplr(t)],[mthetaPowerMG-semthetaPowerMG fliplr(mthetaPowerMG+semthetaPowerMG)],'b','FaceAlpha',0.35,'EdgeColor','none')

ylim([-0.25 0.35])
xlabel('Time (ms)')
ylabel('Theta Power (z)')
meanRT = mean([TaskMetrics.RT]);
xline(meanRT,'--k','LineWidth',1.2)
legend('C-Go','','M-Go','','')
%% Beta correct vs error (No-go)
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,301:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,301:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(301:end);

end

% Downsample to 40 Hz
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,1:25:end,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,1:25:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(1:25:end);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,1:25:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,1:25:end,:,:);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(1:25:end);
end


% Get beta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = squeeze(mean(PDPowers(ii).ngPower(27:41,:,:,:)));
    PDPowers(ii).goPower = squeeze(mean(PDPowers(ii).goPower(27:41,:,:,:)));
    ErrorPowers(ii).ngPowers = squeeze(mean(ErrorPowers(ii).ngPowers(27:41,:,:,:)));
    ErrorPowers(ii).goPowers = squeeze(mean(ErrorPowers(ii).goPowers(27:41,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Correct trials = to # of Error trials
        idxs = randsample(size(PDPowers(subject).ngPower,2),size(ErrorPowers(subject).ngPowers,2),true);
        boot(permi,:,:) = squeeze(mean(PDPowers(subject).ngPower(:,idxs,:),2));
    end
    if numel(size(boot)) == 3
        subjectBoots(subject).Boots = squeeze(mean(boot)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(boot))'; % average over 1000 iterations
    end

    clear boot
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(squeeze(mean(ErrorPowers(subject).ngPowers,2)) - subjectBoots(subject).Boots,2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).ngPower,2),size(ErrorPowers(subject).ngPowers,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).ngPower(:,idxs,:),ErrorPowers(subject).ngPowers);
    
        % Making mapping vector
        mapping = [ones(size(ErrorPowers(subject).ngPowers,2),1); ones(size(ErrorPowers(subject).ngPowers,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    betaPowerCorr(:,ii) = mean(mean(subjectBoots(ii).Boots,2),3);
    betaPowerErr(:,ii) = mean(mean(ErrorPowers(ii).ngPowers,2),3);
end

t = PDPowers(1).Time;

[hl,~] = boundedline(t,mean(betaPowerErr,2),std(betaPowerErr,[],2)./sqrt(size(betaPowerErr,2)),'r',t,mean(betaPowerCorr,2),std(betaPowerCorr,[],2)./sqrt(size(betaPowerCorr,2)),'b');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-202 1000])
ylim([-0.3 0.3])
yticks(-0.3:0.1:0.3);
yticklabels(-0.3:0.1:0.3)
xlabel('Time (ms)')
ylabel('Beta Power (z)')
%
sigClusters = t(zmap==1);


plot(t,mbetaPowerGO,'b','LineWidth',2)
patch([t fliplr(t)],[mbetaPowerGO-sembetaPowerGO fliplr(mbetaPowerGO+sembetaPowerGO)],'b','FaceAlpha',0.35,'EdgeColor','none')
hold on
plot(t,mbetaPowerNG,'r','LineWidth',2)
patch([t fliplr(t)],[mbetaPowerNG-sembetaPowerNG fliplr(mbetaPowerNG+sembetaPowerNG)],'r','FaceAlpha',0.35,'EdgeColor','none')

plot([sigClusters(1) sigClusters(end)],[-0.35 -0.35],'-k','LineWidth',2)
plot(mean(sigClusters),-0.13,'*k','Linewidth', 1.5,'MarkerSize',6)
ylim([-0.4 0.4])
xlabel('Time (ms)')
ylabel('Beta Power (z)')
meanRT = mean([TaskMetrics.RT]);
xline(meanRT,'--k','LineWidth',1.2)
legend('Go','','NoGo','','')

%% Beta correct vs error (Go)

for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,301:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,301:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(301:end);

end

% Downsample to 40 Hz
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,1:25:end,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,1:25:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(1:25:end);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,1:25:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,1:25:end,:,:);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(1:25:end);
end

% Get beta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = squeeze(mean(PDPowers(ii).ngPower(27:41,:,:,:)));
    PDPowers(ii).goPower = squeeze(mean(PDPowers(ii).goPower(27:41,:,:,:)));
    ErrorPowers(ii).ngPowers = squeeze(mean(ErrorPowers(ii).ngPowers(27:41,:,:,:)));
    ErrorPowers(ii).goPowers = squeeze(mean(ErrorPowers(ii).goPowers(27:41,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Correct trials = to # of Error trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(ErrorPowers(subject).goPowers,2),true);
        boot(permi,:,:) = squeeze(mean(PDPowers(subject).goPower(:,idxs,:),2));
    end
    if numel(size(boot)) == 3
        subjectBoots(subject).Boots = squeeze(mean(boot)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(boot))'; % average over 1000 iterations
    end

    clear boot
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(squeeze(mean(ErrorPowers(subject).goPowers,2)) - subjectBoots(subject).Boots,2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(ErrorPowers(subject).goPowers,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).goPower(:,idxs,:),ErrorPowers(subject).goPowers);
    
        % Making mapping vector
        mapping = [ones(size(ErrorPowers(subject).goPowers,2),1); ones(size(ErrorPowers(subject).goPowers,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    betaPowerCorr(:,ii) = mean(mean(subjectBoots(ii).Boots,2),3);
    betaPowerErr(:,ii) = mean(mean(ErrorPowers(ii).goPowers,2),3);
end

t = PDPowers(1).Time;

[hl,~] = boundedline(t,mean(betaPowerErr,2),std(betaPowerErr,[],2)./sqrt(size(betaPowerErr,2)),'r',t,mean(betaPowerCorr,2),std(betaPowerCorr,[],2)./sqrt(size(betaPowerCorr,2)),'b');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-202 1000])
ylim([-0.3 0.3])
yticks(-0.3:0.1:0.3);
yticklabels(-0.3:0.1:0.3)
xlabel('Time (ms)')
ylabel('Beta Power (z)')
%
sigClusters = t(zmap==1);


plot(t,mbetaPowerGO,'b','LineWidth',2)
patch([t fliplr(t)],[mbetaPowerGO-sembetaPowerGO fliplr(mbetaPowerGO+sembetaPowerGO)],'b','FaceAlpha',0.35,'EdgeColor','none')
hold on
plot(t,mbetaPowerNG,'r','LineWidth',2)
patch([t fliplr(t)],[mbetaPowerNG-sembetaPowerNG fliplr(mbetaPowerNG+sembetaPowerNG)],'r','FaceAlpha',0.35,'EdgeColor','none')

plot([sigClusters(1) sigClusters(end)],[-0.35 -0.35],'-k','LineWidth',2)
plot(mean(sigClusters),-0.13,'*k','Linewidth', 1.5,'MarkerSize',6)
ylim([-0.4 0.4])
xlabel('Time (ms)')
ylabel('Beta Power (z)')
meanRT = mean([TaskMetrics.RT]);
xline(meanRT,'--k','LineWidth',1.2)
legend('Go','','NoGo','','')

%% Beta No-go error vs correct Go (motor response in both, rules different)

for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,301:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,301:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(301:end);

end

% Downsample to 40 Hz
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,1:25:end,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,1:25:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(1:25:end);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,1:25:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,1:25:end,:,:);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(1:25:end);
end

% Get beta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = squeeze(mean(PDPowers(ii).ngPower(27:40,:,:,:)));
    PDPowers(ii).goPower = squeeze(mean(PDPowers(ii).goPower(27:41,:,:,:)));
    ErrorPowers(ii).ngPowers = squeeze(mean(ErrorPowers(ii).ngPowers(27:40,:,:,:)));
    ErrorPowers(ii).goPowers = squeeze(mean(ErrorPowers(ii).goPowers(27:40,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Correct trials = to # of Error trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(ErrorPowers(subject).ngPowers,2),true);
        boot(permi,:,:) = squeeze(mean(PDPowers(subject).goPower(:,idxs,:),2));
    end
    if numel(size(boot)) == 3
        subjectBoots(subject).Boots = squeeze(mean(boot)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(boot))'; % average over 1000 iterations
    end

    clear boot
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(squeeze(mean(ErrorPowers(subject).ngPowers,2)) - subjectBoots(subject).Boots,2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).goPower,2),size(ErrorPowers(subject).ngPowers,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).goPower(:,idxs,:),ErrorPowers(subject).ngPowers);
    
        % Making mapping vector
        mapping = [ones(size(ErrorPowers(subject).ngPowers,2),1); ones(size(ErrorPowers(subject).ngPowers,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    betaPowerCorr(:,ii) = mean(mean(subjectBoots(ii).Boots,2),3);
    betaPowerErr(:,ii) = mean(mean(ErrorPowers(ii).ngPowers,2),3);
end

t = PDPowers(1).Time;

[hl,~] = boundedline(t,mean(betaPowerErr,2),std(betaPowerErr,[],2)./sqrt(size(betaPowerErr,2)),'r',t,mean(betaPowerCorr,2),std(betaPowerCorr,[],2)./sqrt(size(betaPowerCorr,2)),'b');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-202 1000])
ylim([-0.3 0.3])
yticks(-0.3:0.1:0.3);
yticklabels(-0.3:0.1:0.3)
xlabel('Time (ms)')
ylabel('Beta Power (z)')
hold on
plot([sigClusters(1) sigClusters(end)],[-0.29 -0.29],'-k','LineWidth',2)
plot(mean(sigClusters),-0.275,'*k','Linewidth', 1.5,'MarkerSize',6)


legend('','','Error NoGo','Correct Go','','')

%% Beta Go error vs correct No-go (motor response in neither, rules different)

for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,301:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,301:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(301:end);
end

% Downsample to 40 Hz
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,1:25:end,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,1:25:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(1:25:end);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,1:25:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,1:25:end,:,:);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(1:25:end);
end

% Get beta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = squeeze(mean(PDPowers(ii).ngPower(27:41,:,:,:)));
    PDPowers(ii).goPower = squeeze(mean(PDPowers(ii).goPower(27:41,:,:,:)));
    ErrorPowers(ii).ngPowers = squeeze(mean(ErrorPowers(ii).ngPowers(27:41,:,:,:)));
    ErrorPowers(ii).goPowers = squeeze(mean(ErrorPowers(ii).goPowers(27:41,:,:,:)));
end

flipidxs = [];
for subject = 1:numel(PDPowers) % loop over subjects
    if size(ErrorPowers(subject).goPowers,2) > size(PDPowers(subject).ngPower,2)
        flipidxs = [flipidxs;subject];
        for permi = 1:1000 % repeat subsampling procedure 1000 times
    
            % in each electrode, subsample set of error go trials = to # of correct nogo trials
            idxs = randsample(size(ErrorPowers(subject).goPowers,2),size(PDPowers(subject).ngPower,2),true);
            boot(permi,:,:) = squeeze(mean(ErrorPowers(subject).goPowers(:,idxs,:),2));
        end
        if numel(size(boot)) == 3
            subjectBoots(subject).Boots = squeeze(mean(boot)); % average over 1000 iterations     
        else
            subjectBoots(subject).Boots = squeeze(mean(boot))'; % average over 1000 iterations
        end
        clear boot
    else
        for permi = 1:1000 % repeat subsampling procedure 1000 times
    
            % in each electrode, subsample set of Correct nogo trials = to # of Error go trials
            idxs = randsample(size(PDPowers(subject).ngPower,2),size(ErrorPowers(subject).goPowers,2),true);
            boot(permi,:,:) = squeeze(mean(PDPowers(subject).ngPower(:,idxs,:),2));
        end
        if numel(size(boot)) == 3
            subjectBoots(subject).Boots = squeeze(mean(boot)); % average over 1000 iterations     
        else
            subjectBoots(subject).Boots = squeeze(mean(boot))'; % average over 1000 iterations
        end
        clear boot
    end
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    if ismember(subject,flipidxs)
        trueDiffs(:,subject) = mean(subjectBoots(subject).Boots - squeeze(mean(PDPowers(subject).ngPower,2)),2);
    else
        trueDiffs(:,subject) = mean(squeeze(mean(ErrorPowers(subject).goPowers,2)) - subjectBoots(subject).Boots,2);
    end
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)
        if ismember(subject,flipidxs)
            % Get indices for subsampled set of Go trials = to # of No-Go trials
            idxs = randsample(size(ErrorPowers(subject).goPowers,2),size(PDPowers(subject).ngPower,2),true);
    
            % Combine trials
            allTrials = cat(2,PDPowers(subject).ngPower,ErrorPowers(subject).goPowers(:,idxs,:));
        
            % Making mapping vector
            mapping = [ones(size(PDPowers(subject).ngPower,2),1); ones(size(PDPowers(subject).ngPower,2),1)+1];
            
            % Shuffle the mapping vector
            new_map = mapping(randperm(numel(mapping)));
        
            % Calculate mean difference
            diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
        else
            % Get indices for subsampled set of No-Go trials = to # of Go trials
            idxs = randsample(size(PDPowers(subject).ngPower,2),size(ErrorPowers(subject).goPowers,2),true);
    
            % Combine trials
            allTrials = cat(2,PDPowers(subject).ngPower(:,idxs,:),ErrorPowers(subject).goPowers);
        
            % Making mapping vector
            mapping = [ones(size(ErrorPowers(subject).goPowers,2),1); ones(size(ErrorPowers(subject).goPowers,2),1)+1];
            
            % Shuffle the mapping vector
            new_map = mapping(randperm(numel(mapping)));
        
            % Calculate mean difference
            diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
        end
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    if ismember(ii,flipidxs)
        betaPowerCorr(:,ii) = mean(mean(PDPowers(ii).ngPower,2),3);
        betaPowerErr(:,ii) = mean(mean(subjectBoots(ii).Boots,2),3);
    else
        betaPowerCorr(:,ii) = mean(mean(subjectBoots(ii).Boots,2),3);
        betaPowerErr(:,ii) = mean(mean(ErrorPowers(ii).goPowers,2),3);
    end
end

t = PDPowers(1).Time;

[hl,~] = boundedline(t,mean(betaPowerErr,2),std(betaPowerErr,[],2)./sqrt(size(betaPowerErr,2)),'r',t,mean(betaPowerCorr,2),std(betaPowerCorr,[],2)./sqrt(size(betaPowerCorr,2)),'b');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-202 1000])
ylim([-0.3 0.3])
yticks(-0.3:0.1:0.3);
yticklabels(-0.3:0.1:0.3)
xlabel('Time (ms)')
ylabel('Beta Power (z)')
hold on
plot([sigClusters(1) sigClusters(end)],[-0.29 -0.29],'-k','LineWidth',2)
plot(mean(sigClusters),-0.275,'*k','Linewidth', 1.5,'MarkerSize',6)


legend('','','Error Go','Correct NoGo','','')
%% Theta correct vs error (No-go)

for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = PDPowers(ii).ngPower(:,426:1626,:,:);
    PDPowers(ii).goPower = PDPowers(ii).goPower(:,426:1626,:,:);
    ErrorPowers(ii).ngPowers = ErrorPowers(ii).ngPowers(:,301:end,:,:);
    ErrorPowers(ii).goPowers = ErrorPowers(ii).goPowers(:,301:end,:,:);
    PDPowers(ii).Time = PDPowers(ii).Time(426:1626);
    ErrorPowers(ii).Time = ErrorPowers(ii).Time(301:end);
end

% Get theta
for ii = 1:numel(PDPowers)
    PDPowers(ii).ngPower = squeeze(mean(PDPowers(ii).ngPower(1:17,:,:,:)));
    PDPowers(ii).goPower = squeeze(mean(PDPowers(ii).goPower(1:17,:,:,:)));
    ErrorPowers(ii).ngPowers = squeeze(mean(ErrorPowers(ii).ngPowers(1:17,:,:,:)));
    ErrorPowers(ii).goPowers = squeeze(mean(ErrorPowers(ii).goPowers(1:17,:,:,:)));
end

for subject = 1:numel(PDPowers) % loop over subjects
    for permi = 1:1000 % repeat subsampling procedure 1000 times

        % in each electrode, subsample set of Correct trials = to # of Error trials
        idxs = randsample(size(PDPowers(subject).ngPower,2),size(ErrorPowers(subject).ngPowers,2),true);
        boot(permi,:,:) = squeeze(mean(PDPowers(subject).ngPower(:,idxs,:),2));
    end
    if numel(size(boot)) == 3
        subjectBoots(subject).Boots = squeeze(mean(boot)); % average over 1000 iterations     
    else
        subjectBoots(subject).Boots = squeeze(mean(boot))'; % average over 1000 iterations
    end

    clear boot
end

trueDiffs = NaN(length(PDPowers(1).Time),numel(PDPowers));

for subject = 1:numel(PDPowers)
    trueDiffs(:,subject) = mean(squeeze(mean(ErrorPowers(subject).ngPowers,2)) - subjectBoots(subject).Boots,2);
end

trueDiff = mean(trueDiffs,2);

clear trueDiffs subject idxs ii

% Permutations for testing H0 = 0 at every T-F point

npermutations = 1000;
alpha = 0.05;
perm_diffs = NaN(npermutations,size(trueDiff,1));

for permi = 1:npermutations
    
    diff = NaN(size(trueDiff,1),numel(PDPowers)); % initialize diff matrix

    for subject = 1:numel(PDPowers)

        % Get indices for subsampled set of Go trials = to # of No-Go trials
        idxs = randsample(size(PDPowers(subject).ngPower,2),size(ErrorPowers(subject).ngPowers,2),true);

        % Combine trials
        allTrials = cat(2,PDPowers(subject).ngPower(:,idxs,:),ErrorPowers(subject).ngPowers);
    
        % Making mapping vector
        mapping = [ones(size(ErrorPowers(subject).ngPowers,2),1); ones(size(ErrorPowers(subject).ngPowers,2),1)+1];
        
        % Shuffle the mapping vector
        new_map = mapping(randperm(numel(mapping)));
    
        % Calculate mean difference
        diff(:,subject) = mean(squeeze(mean( allTrials(:,new_map==2,:) ,2)) - squeeze(mean( allTrials(:,new_map==1,:) ,2)),2);
    end
    
    % Store values for this iteration
    perm_diffs(permi,:) = mean(diff,2);
end


% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_diffs))';
std_h0 = squeeze(std(perm_diffs))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_diffs(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end

for ii = 1:numel(PDPowers)
    thetaPowerCorr(:,ii) = mean(mean(subjectBoots(ii).Boots,2),3);
    thetaPowerErr(:,ii) = mean(mean(ErrorPowers(ii).ngPowers,2),3);
end

t = PDPowers(1).Time;

[hl,~] = boundedline(t,mean(thetaPowerErr,2),std(thetaPowerErr,[],2)./sqrt(size(thetaPowerErr,2)),'r',t,mean(thetaPowerCorr,2),std(thetaPowerCorr,[],2)./sqrt(size(thetaPowerCorr,2)),'b');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-202 1000])
ylim([-0.2 0.7])
xlabel('Time (ms)')
ylabel('Low Frequency Power (z)')

%
sigClusters = t(zmap==1);


plot(t,mbetaPowerGO,'b','LineWidth',2)
patch([t fliplr(t)],[mbetaPowerGO-sembetaPowerGO fliplr(mbetaPowerGO+sembetaPowerGO)],'b','FaceAlpha',0.35,'EdgeColor','none')
hold on
plot(t,mbetaPowerNG,'r','LineWidth',2)
patch([t fliplr(t)],[mbetaPowerNG-sembetaPowerNG fliplr(mbetaPowerNG+sembetaPowerNG)],'r','FaceAlpha',0.35,'EdgeColor','none')

plot([sigClusters(1) sigClusters(end)],[-0.35 -0.35],'-k','LineWidth',2)
plot(mean(sigClusters),-0.13,'*k','Linewidth', 1.5,'MarkerSize',6)
ylim([-0.4 0.4])
xlabel('Time (ms)')
ylabel('Beta Power (z)')
meanRT = mean([TaskMetrics.RT]);
xline(meanRT,'--k','LineWidth',1.2)
legend('Go','','NoGo','','')
