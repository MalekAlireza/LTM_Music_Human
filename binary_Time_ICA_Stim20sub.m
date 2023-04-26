%%% This code is writen to finilized the time-lock analysis 
% binary_time_ICA_first_post
clear 
path(pathdef);
clc
close all
%% add path (EEGLAB - Fieldtrip - Chronux)
matlabrc
current_path = 'F:/EEGdata/ClassicalSongParadigm/';
addpath 'F:/EEGdata/ClassicalSongParadigm/fieldtrip-20191008'
ft_defaults
%% subjects
subjectsFolders =  {'code_2019_06_05_Sebastian',...
                    'code_2019_06_17_Florian',...
                    'code_2019_06_18_Augosto',...
                    'code_2019_07_02_Lukas',...
                    'code_2019_07_09_Sai',...
                    'code_2019_07_27_Mostafa',...
                    'code_2019_07_29_George',...
                    'code_2019_07_30_Farnam',...
                    'code_2019_07_30_Poorya',...
                    'code_2019_08_04_Behnam',...
                    'code_2019_08_04_Hadi',...
                    'code_2019_08_05_Ali',...
                    'code_2019_08_07_Babak',...
                    'code_2019_08_08_Mahdi',...
                    'code_2019_09_24_Julian',...
                    'code_2022_01_13_Mohsen',...
                    'code_2022_01_14_Foad',...
                    'code_2022_01_16_Petar',...
                    'code_2022_01_18_Ehsan',...
                    'code_2022_01_27_Arturo'};
subjectsNames ={'Sebastian',...
                'Florian',...
                'Augosto',...
              	'Lukas',...
              	'Sai',...
              	'Mostafa',...
               	'George',...
               	'Farnam',...
               	'Poorya',...
              	'Behnam',...
               	'Hadi',...
                'Ali',...
               	'Babak',...
             	'Mahdi',...
              	'Julian',...
                'Mohsen',...
                'Foad',...
                'Peter',...
                'Ehsan',...
                'Arturo'};
                
removingBADsubjects = [ ];
subjectsFolders(removingBADsubjects) = [];
subjectsNames(removingBADsubjects) = [];
%% Worked perfect
uselabel = [1,  7,  6, 	2,  5,      NaN;...   	%S1     OK
            1,  7,  6, 	5,  NaN,    NaN;...  	%S2     OK
            1,  7,  2, 	3,  NaN,    NaN;... 	%S3     OK
            1,  7,  6,	2,  5,      NaN;...     %S4     OK
            1,  7,  6, 	2,  5,      NaN;...     %S5     OK
            1,  7,  2,  3,  NaN,    NaN;... 	%S6     
            1,  7,  6, 	2,  NaN,    NaN;...     %S7     OK
            1,  7,  6,  2,  5,      NaN;...   	%S8     OK
            1,  7,  6,  2,  5,      NaN;...     %S9     OK
            1,  7,  6,  2,  5,      NaN;...   	%S10    OK
            1,  7,  6,  2,  NaN,    NaN;...     %S11    OK
            1,  7,  6,  2,  5,      3;  ...     %S12    OK
            1,  7,  6,  2,  5,      NaN;...     %S13    OK
            1,  7,  6,  2,  NaN,  	NaN;...  	%S14            
            1,  7,  6,  5,  NaN,    NaN; ...  	%S15    OK  
            1,  7,  6,  2,  5,      NaN;... 	%S16    OK
            1,  7,  6,  2,  3,      NaN;...     %S17    OK
            1,  7,NaN,  NaN,  NaN,  NaN;...  	%S18    OK
            1,  7,  6,  2,  3,      NaN;... 	%S19    OK
            1,  7,  NaN,NaN,NaN,    NaN];     	%S20    OK
%% temporal response calculation
for index_sub = 1 : length(subjectsFolders)
    load(fullfile(current_path,subjectsFolders{index_sub},strcat('S_',subjectsNames{index_sub},'_Binary_ICA_clean')));%% load tenStim_seg_AR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% baseline average
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg = [];
    cfg.demean        = 'yes';%'no' or 'yes', whether to apply baseline correction (default = 'no')
    cfg.lpfilter   = 'yes';                              % apply lowpass filter
    cfg.lpfreq     = 30;                                 % lowpass at 35 Hz.
  	cfg.baselinewindow = [-0.2 -0.001];
    data_eeg = ft_preprocessing(cfg, data_eeg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% extract first and post
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg = [];
    cfg.latency = [-0.2 1];
    data_eeg_first = ft_selectdata(cfg, data_eeg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% EEG channel selection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [eeg_channel] = ft_channelselection('eeg', data_eeg.label);
    cfg = [];
    cfg.channel = eeg_channel;
    [data_eeg_first] = ft_selectdata(cfg, data_eeg_first);
    clear data_eeg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% time_lock analysis 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for index_task = 1:2
        if index_task == 1
            cfg.trials = [];
            ind_con = [];
            for indlow = uselabel(index_sub,uselabel(index_sub,:)<4)
                cfg.trials  = cat(1,cfg.trials, find(data_eeg_first.trialinfo(:,end) == indlow));
                ind_con     = cat(1,ind_con,    find(data_eeg_first.trialinfo(:,end) == indlow));
            end
        elseif index_task == 2
            cfg.trials = [];
            ind_con = [];
            for indlow = uselabel(index_sub,uselabel(index_sub,:)>4)
                cfg.trials  = cat(1,cfg.trials, find(data_eeg_first.trialinfo(:,end) == indlow));
                ind_con     = cat(1,ind_con,    find(data_eeg_first.trialinfo(:,end) == indlow));
            end
        else
            error('There is only two tasks')
        end
     	ERP_first{index_sub,index_task} = ft_timelockanalysis(cfg, data_eeg_first);
    end  
    for index_task = 1:2
        cfg = [];
        cfg.baselinewindow = [-0.2 -0.001];
        ERP_first{index_sub,index_task} = ft_preprocessing(cfg, ERP_first{index_sub,index_task});
    end
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    ERP_diff_first{index_sub}   = ft_math(cfg, ERP_first{index_sub,1}, ERP_first{index_sub,2});
end
%% Statistical analysis for ERP
alpha = 0.05;
rem_sub = nchoosek(1:20,20);
% rem_sub(ii)=[];
% calculate grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
ERP_first_grand{1}         = ft_timelockgrandaverage(cfg,ERP_first{rem_sub,1});
ERP_first_grand{2}         = ft_timelockgrandaverage(cfg,ERP_first{rem_sub,2});

cfg_neighb        = [];
cfg_neighb.method = 'template';
neighbours        = ft_prepare_neighbours(cfg_neighb, ERP_first_grand{1});
cfg = [];
cfg.channel = {'EEG'};
cfg.latency = [0 1]; %% you need to change j = [0:timestep:1]; 
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = alpha; 
cfg.clusterstatistic = 'maxsize';
cfg.minnbchan = 4;
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = alpha;
cfg.resampling    = 'permutation';%permutation, bootstrap
cfg.numrandomization = 1000;

subj = length(rem_sub);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
[stat] = ft_timelockstatistics(cfg, ERP_first{rem_sub,1}, ERP_first{rem_sub,2});
%%%%% plot
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
ERP_first_grand_sub = ft_math(cfg, ERP_first_grand{1}, ERP_first_grand{2});


% define parameters for plotting
timestep = 0.05;      %(in seconds)
sampling_rate = 1000;
sample_count = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
% get relevant (significant) values
try
    pos_cluster_pvals = [stat.posclusters(:).prob];
catch
end
try
    neg_cluster_pvals = [stat.negclusters(:).prob];
catch
end
% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = alpha; end % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.

try
    pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
    pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
catch
    pos = zeros(size(stat.posclusterslabelmat));
end

try
    neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
    neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
catch
    neg = zeros(size(stat.negclusterslabelmat));
end
% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(ERP_first_grand_sub.label, stat.label);

%%% plot
for k = 1:20
    figure('units','normalized','outerposition',[0 0 1 1])
    cfg = [];
    cfg.xlim=[j(k) j(k+1)];
    
    pos_int = zeros(numel(ERP_first_grand_sub.label),1);
    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
    neg_int = zeros(numel(ERP_first_grand_sub.label),1);
    neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

    cfg.zlim      = [-1.5 1.5];
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos_int | neg_int);%%
    cfg.highlightcolor   = [1 1 0];
    cfg.comment      = 'no'; %%no  xlim
    cfg.commentpos   = 'layout';
    cfg.layout       = 'easycapM1.mat';
    cfg.colorbar     = 'no';%%East  no
    cfg.colormap     = 'jet';
    cfg.markersize   = 6;
    cfg.highlightsymbol    = 'o';
    cfg.highlightsize      = 10;
    cfg.highlightfontsize  = 10;
    ft_topoplotER(cfg, ERP_first_grand_sub);
%     print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/ERP/erpTopoBinary','_',num2str(100*j(k))))
%     close all
end
%sgtitle(num2str(rem_sub(ind,:)));
effectSizeERPpaired(ERP_first(1:20,2),ERP_first(1:20,1),[0.4 0.45],{'FCz','FC1','FC2','Cz','C1','C2'},20)
%% Look at the whole erp at one glance
cfg = [];
cfg.showlabels  = 'yes';
cfg.layout      = 'easycapM1.mat';
figure; ft_multiplotER(cfg, ERP_first_grand{1}, ERP_first_grand{2})
%% STD shade (average erp for all subjects)
chan = {'Fz','FCz','FC1','FC2','Cz','C1','C2','Pz'};
% chan = {'P3','P5','P7','CP3','CP5','TP7','C3','C5','T7'};
fam = zeros(20, size(ERP_first{1,1}.avg,2),length(chan));
unfam = zeros(20, size(ERP_first{1,2}.avg,2),length(chan));
count = 1;
for i = 1:20
    for j = 1:length(chan)
        fam(count,:,j) = ERP_first{i,2}.avg(strcmp(ERP_first{i,2}.label , chan{j}),:);
        unfam(count,:,j) = ERP_first{i,1}.avg(strcmp(ERP_first{i,1}.label , chan{j}),:);
    end
    count = count +1;
end
fam = mean(fam,3);
unfam = mean(unfam,3);
ymax = 16;
figure('units','normalized','outerposition',[0 0 1 1])
hold on
time = [0.4 0.45];
rectangle('Position',[time(1) -3.5 (time(2)-time(1)) ymax],'FaceColor',[0.5 0.5 0.5]);
hold on
stdshade(fam, 0.1, 'b', ERP_first{1,1}.time , 3);
hold on;
stdshade(unfam, 0.1, 'r', ERP_first{1,1}.time , 3);
xticks([-0.2:0.1:1]);
grid;
ylim([-3.5 9])
xlim([-0.2 1])
ylabel('Amplitude(\muV)')
xlabel('Time(s)')
legend('Familiarity SEM','Familiarity','Unfamiliarity SEM','Unfamiliarity')
set(gca,'FontSize',25)
print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/ERP/erpAllSub'))
% print('-dtiff','-r500',strcat('F:\EEGdata\ClassicalSongParadigm\Results20Sub\SuppERP/SUPPerpAllSubLeftposterior'))
%% ERP for each subject individually
chan = {'FCz','FC1','FC2','Cz','C1','C2'};
time = [0.4 0.45];
% Scaling of the vertical axis for the plots below
ymax = 16;
figure;
ss=[];
sss=[];
for isub = [1:20]
    fam = zeros(length(chan),size(ERP_first{1,1}.avg,2));
    unfam = zeros(length(chan),size(ERP_first{1,2}.avg,2));
    for inChan = 1 : length(chan)
        fam(inChan,:) = ERP_first{isub,2}.avg(strcmp(ERP_first{isub,2}.label , chan{inChan}),:);
        unfam(inChan,:) = ERP_first{isub,1}.avg(strcmp(ERP_first{isub,1}.label , chan{inChan}),:);
    end
    figure;
    % plot the lines in front of the rectangle
    plot(ERP_first{isub,1}.time,mean(fam,1),'b','LineWidth',2);
    hold on;
    plot(ERP_first{isub,2}.time,mean(unfam,1),'r','LineWidth',2);
    title(strcat('Subject ',num2str(isub)))
    ylabel('Amplitude (\muV)')
    xlabel('Time(s)')
    set(gca,'FontSize',15)
    grid;
    ylim([-16 23])
    xlim([-.2 1])
    print('-dtiff','-r1000',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/ERP/erpSub','_',num2str(isub)))
    close all
end