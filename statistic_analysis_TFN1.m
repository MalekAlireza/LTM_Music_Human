function [stat] = statistic_analysis_TFN1(data, sel)
%%% sel = a vector to determine the sub : [ 1 2 4 5 8 9 10] 
%%% => removing sub 3 6 7  
for index_task = 1 : 2
    cfg = [];
    cfg.keepindividual = 'yes';
    grandavg{index_task} = ft_freqgrandaverage(cfg, data{sel,index_task});
end
% grandavg{1,1}.powspctrm = grandavg{1,1}.powspctrm - grandavg{1,2}.powspctrm;
% grandavg{1,2}.powspctrm = grandavg{1,2}.powspctrm - grandavg{1,2}.powspctrm;
cfg = [];
cfg.channel          = {'EEG'};
cfg.latency          = [-0.7 5];%% [-0.7 5] and [-0.0 10]
cfg.frequency        = [4 40];%% [3 40] and [3 30]
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.method    = 'analytic';
% cfg.correctm  = 'fdr';
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.1;
cfg.clusterstatistic = 'maxsize';
cfg.minnbchan        = 1;
cfg.tail             = 1;
% cfg.correcttail      = 'prob';
cfg.clustertail      = 1;
cfg.alpha            = 0.05; %% 0.05 and  0.1
cfg.numrandomization = 500;%% 200 and 1000
% specifies with which sensors other sensors can form clusters
% prepare_neighbours determines what sensors may form clusters DB
cfg_neighb = [];
cfg_neighb.method = 'template';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, grandavg{1});

subj = length(sel);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, grandavg{1}, grandavg{2});

end
