function [stat] = statistic_analysis_TFN4(grandavg, sel)
% % %%% sel = a vector to determine the sub : [ 1 2 4 5 8 9 10] 
% % %%% => removing sub 3 6 7  
% % for index_task = 1 : 2
% %     cfg = [];
% %     cfg.keepindividual = 'yes';
% %     grandavg{index_task} = ft_freqgrandaverage(cfg, data{sel,index_task});
% % end
% grandavg{1,1}.powspctrm = grandavg{1,1}.powspctrm - grandavg{1,2}.powspctrm;
% grandavg{1,2}.powspctrm = grandavg{1,2}.powspctrm - grandavg{1,2}.powspctrm;
cfg = [];
cfg.channel          = {'EEG'};
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.1;
cfg.clusterstatistic = 'maxsize';
cfg.minnbchan        = 3;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 3000;
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
cfg = [];
cfg.parameter = 'stat';
cfg.layout = 'easycapM1.mat';
cfg.alpha  = 0.05;
%cfg.zlim = [-4 4];
cfg.subplotsize = [1 1];
cfg.highlightseries         = {'on', 'on', 'on'};
cfg.highlightsymbolseries   =['*', 'x', '+'];
cfg.highlightcolorpos       = [1 0 1];
cfg.highlightcolorneg       = [1 1 0];
cfg.highlightsizeseries     = [16 16 12 12 12];
cfg.marker          =  'on';
cfg.colorbar        = 'no';%%East  no
cfg.colormap        = 'jet';
cfg.markersize      = 8;
cfg.comment         = 'no'; %%no  xlim
cfg.commentpos      = 'layout';
cfg.highlightsymbol    = 'o';
cfg.highlightsize      = 10;
cfg.highlightfontsize  = 10;
ft_clusterplot(cfg, stat);
% stat.stat = mean(grandavg{1,1}.powspctrm-grandavg{1,2}.powspctrm)./...
%     std(grandavg{1,1}.powspctrm-grandavg{1,2}.powspctrm,[],1);
figure;
stat.stat = stat.stat+0.75;
cfg.zlim = [-1.5 1.5];
ft_topoplotER(cfg, stat)
end

