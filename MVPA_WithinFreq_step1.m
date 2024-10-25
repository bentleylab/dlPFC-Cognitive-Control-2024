% Path setup
PATH_MVPA   = '/Users/anaskhan/Desktop/Apps/MVPA-Light-master/startup';
PATH_fs = '/Users/anaskhan/Desktop/Apps/fieldtrip';
addpath(PATH_MVPA)
startup_MVPA_Light

load MVPAData.mat
conf.frequency      = {'theta','beta'};

for freq = 1:size(conf.frequency,2) 
    if freq == 1 % Low frequency (theta)
        frexidxs = 1:17;
    else
        frexidxs = 27:40;
    end
        for sbj = 1:13 %subject
              
              data = squeeze(mean(AllPowers(sbj).Powers(frexidxs,:,:,:)));
              data = permute(data,[2 3 1]);

              clabel = double([AllPowers(sbj).Trials.Condition] == "GO")';
              clabel = clabel + 1;
                % classification across time
                cfg             = [];
                cfg.metric      = {'acc','auc'};
                cfg.classifier  = 'svm';
                [~,result_across_time{sbj}] = mv_classify(cfg,data,clabel);
        end
        
        save(fullfile('/Users/anaskhan/Documents/Bentley Lab/Analysis/Data/DLPFC/PD','D1_MVPA_step1',['MVPA_step1_',conf.frequency{freq},'.mat']),'result_timeXtime','result_across_time');
end