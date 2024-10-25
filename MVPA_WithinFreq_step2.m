
%% Managing perception-action representations through the interplay of alpha and theta band activity
% Paul Wendiggensen, Astrid Prochnow, Charlotte Pscherer, Alexander Münchau, Christian Frings, Christian Beste
% 
% D2 - MVPA on source-level alpha and theta (temporal), group-level
% Paul Wendiggensen & Astrid Prochnow, 2022
% based on code by Moritz Mückschel


% MVPA 2nd step Group Analysis
% Perform level 2 MVPA analysis, i.e. group level
%
%
% Initial version 2021-11 by Moritz Mückschel
%
% v2 2021-11 Moritz Mückschel
%
% v3 2021-12-09 Moritz Mückschel
% * change figure export from 'saveas' to 'print', support specification of
%   export resolution
% * Support plotting in different formats
% * Allow specifying figure font name
% 
% v4 2021-12-20 Astrid Prochnow
% * save time and AUC range of significant classification for reporting
%
% v6 2022-02-28 Moritz Mückschel
% * do not use hard coded alpha values, specify in conf instead
%
% v7 2022-04-11
%
% v7_1 2022-04-28
% * disabled unmasked temporal generalization plot
%

clc
clearvars
close all

%% CONFIG

% path toolbox
conf.path_mvpalight = '/Users/anaskhan/Desktop/Apps/MVPA-Light-master';
conf.path_fieldtrip = '/Users/anaskhan/Desktop/Apps/fieldtrip';

conf.comparison = {...
    'MVPA_step1_beta',...
    'MVPA_step1_theta',...
    };

% String indicating baseline settings
conf.baseline_string = '';

% Folder containing output of step1 analysis (i.e. level 1 individual
% subjects)
conf.step1outputfolder = ['/Users/anaskhan/Documents/Bentley Lab/Analysis/Data/DLPFC/PD','/D1_MVPA_step1'];

% Plotting folder
conf.plotfolder = ['/Users/anaskhan/Documents/Bentley Lab/Analysis/Data/DLPFC/PD','/D2_MVPA_step2/'];

% Number of permutations
conf.n_permutations = 1000;


% clustercritval is a parameter that needs to be selected by the user. In a
% Level 2 (group level) test, it represents the critical cutoff value for
% the statistic. Here, we selected Wilcoxon, so clustercritval corresponds
% to the cutoff value for the z-statistic which is obtained by a normal
% approximation
conf.clustercritval = 0.05;

% cfg.stat_alpha
conf.clusteralpha = 0.05;


% Figure formattype 
% Cell of chars or strings, speciying the figure formattype as used by the 
% 'print' command, e.g. '-dpng' for PNG or '-dsvg' for SVG. Multiple
% formattypes are supported. Example: 
% conf.figure_formattype = {'-dsvg','-dpng'};
conf.figure_formattype = {'-dsvg','-dpng'};

% Figure export resolution
% Resolution, specified as a character vector or a string containing -r and
% an integer value indicating the resolution in dots per inch. For example,
% '-r300' sets the output resolution to 300 dots per inch.
conf.figure_resolution = '-r600';

% Font name for figures
conf.figure_fontname = 'Arial';

%%

% Add fieldtrip and MVPA light toolbox
restoredefaultpath;
addpath(conf.path_fieldtrip);
ft_defaults;
addpath([conf.path_mvpalight '/startup'])
addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections'))
startup_MVPA_Light

% Create outputfolder
if ~exist(conf.plotfolder,'dir')
    mkdir(conf.plotfolder)
end

%% Loop comparisons
for cond = 1:numel(conf.comparison)
    % Load data
    data = load([conf.step1outputfolder filesep conf.comparison{cond} conf.baseline_string '.mat']);
         
    time_range = -200:25:1000;
    %% across time
    
    % calculate/plot the merge and average
    results = data.result_across_time;
    result_merge = mv_combine_results(results, 'merge');
    result_merge = mv_select_result(result_merge, 'auc'); % select AUC only
    % mv_plot_result(result_merge)
    
    % Instead of plotting all subjects together, we can also calculate the
    % grand average across subjects and plot this. To this end, we only need to
    % replace 'merge' by 'average' when calling mv_combine_results. Note that
    % the shaded area is now the standard deviation across subjects.
    result_average = mv_combine_results(results, 'average');
    result_average = mv_select_result(result_average, 'auc');
    % mv_plot_result(result_average);
    % ylim([0 1]);
    
    cfg_stat = [];
    cfg_stat.metric          = 'auc';
    cfg_stat.test            = 'permutation';
    cfg_stat.correctm        = 'cluster';  % correction method is cluster
    cfg_stat.n_permutations  = conf.n_permutations;
    
    % Clusterstatistic is the actual statistic used for the clustertest.
    % Normally the default value 'maxum' is used, we are setting it here
    % explicitly for clarity. Maxsum adds up the statistics calculated at each
    % time point (the latter are set below using cfg_stat.statistic)
    cfg_stat.clusterstatistic = 'maxsum';
    cfg_stat.alpha           = conf.clusteralpha; % use standard significance threshold of 5%
    
    % Level 2 stats design: we have to choose between within-subject and
    % between-subjects. Between-subjects is relevant when there is two
    % different experimental groups (eg patients vs controls) and we want to
    % investigate whether their MVPA results are significantly different. Here,
    % we have only one group and we want to see whether the AUC is
    % significantly different from a null value, hence the statistical design
    % is within-subject
    cfg_stat.design          = 'within';
    % cfg_stat.statistic defines how the difference between the AUC values and
    % the null is calculated at each time point (across subjects).
    % We can choose t-test or its nonparametric counterpart Wilcoxon test. We
    % choose Wilcoxon here.
    cfg_stat.statistic       = 'wilcoxon';
    % The null value for AUC (corresponding to a random classifier) is 0.5
    cfg_stat.null            = 0.5;
    
    % clustercritval is a parameter that needs to be selected by the user. In a
    % Level 2 (group level) test, it represents the critical cutoff value for
    % the statistic. Here, we selected Wilcoxon, so clustercritval corresponds
    % to the cutoff value for the z-statistic which is obtained by a normal
    % approximation
    cfg_stat.clustercritval = abs(norminv(conf.clustercritval/2));


    
    stat_level2 = mv_statistics(cfg_stat, results);
    
    % save time and AUC range and significant classification for reporting
    AcrossTime.comparison = conf.comparison(cond);
    AcrossTime.time_range = time_range;
    AcrossTime.AUCvalues = result_average.perf';
    AcrossTime.significant = stat_level2.mask;
    save([char(conf.plotfolder),char('/stats_across_time_'),char(conf.comparison(cond)),char('.mat')],'AcrossTime');
    
    % plot the grand average result again and indicate the cluster in bold
    h = mv_plot_result(result_average, time_range, 'mask', stat_level2.mask);
    set(gca, 'FontName', conf.figure_fontname)
    set(gcf,'color','w');
    set(gcf, 'Units', 'Centimeters', 'Position', [0, 0, 16, 16], 'PaperUnits', 'Centimeters', 'PaperSize', [16,16]);
    ylim([0 1]);
    legend show;
    
    if contains(conf.comparison(cond),'beta')
        title('Beta Power MVPA');
        xlabel('Time (ms)')
        ylabel('AUC')
        grid off
        legend off
    else
        title('LF Power MVPA');
        xlabel('Time (ms)')
        ylabel('AUC')
        grid off
        legend off
    end

    for ii = 1:numel(conf.figure_formattype)
        print(gcf,[conf.plotfolder filesep conf.comparison{cond} conf.baseline_string '_ACROSS_TIME'],conf.figure_formattype{ii},conf.figure_resolution);
    end
    if numel(conf.comparison) > 1
        close
    end
end