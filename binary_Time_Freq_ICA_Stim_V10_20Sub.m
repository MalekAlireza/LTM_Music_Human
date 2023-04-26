%%% This code is writen to finilized the time frequency analysis 
clear
path(pathdef);
clc
close all
%% add path (EEGLAB - Fieldtrip - Chronux)
matlabrc
current_path = 'F:/EEGdata/ClassicalSongParadigm/';
addpath 'F:/EEGdata/ClassicalSongParadigm/fieldtrip-20191008/'
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
%% time frequency response with Hanning
total{length(subjectsFolders),2} = []; 
for index_sub = 1:length(subjectsFolders)
    load(fullfile(current_path,subjectsFolders{index_sub},strcat('S_',subjectsNames{index_sub},'_Binary_ICA_clean')));%% load tenStim_seg_AR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% baseline average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg = [];
    cfg.demean        = 'yes';%'no' or 'yes', whether to apply baseline correction (default = 'no')
    cfg.baselinewindow = [-0.2 -0.001];
    data_eeg = ft_preprocessing(cfg, data_eeg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% EEG channel selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [eeg_channel] = ft_channelselection('eeg', data_eeg.label);
    cfg = [];
    cfg.channel = eeg_channel;
    [data_eeg] = ft_selectdata(cfg, data_eeg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% frequency analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for index_task = 1:2
        cfg = [];
        cfg.output = 'pow';
        cfg.method = 'mtmconvol';
        cfg.taper  = 'hanning';
        cfg.keeptrials  = 'yes';
        if index_task == 1
            cfg.trials  = find(data_eeg.trialinfo(:,end) == 1 | data_eeg.trialinfo(:,end) == 2 | data_eeg.trialinfo(:,end) == 3);
            ind_con     = find(data_eeg.trialinfo(:,end) == 1 | data_eeg.trialinfo(:,end) == 2 | data_eeg.trialinfo(:,end) == 3);
        elseif index_task == 2
            cfg.trials  = find(data_eeg.trialinfo(:,end) == 6 | data_eeg.trialinfo(:,end) == 7 | data_eeg.trialinfo(:,end) == 5);
            ind_con     = find(data_eeg.trialinfo(:,end) == 6 | data_eeg.trialinfo(:,end) == 7 | data_eeg.trialinfo(:,end) == 5);
        else
            error('There is only two tasks')
        end
        cfg.foi          = 3:0.5:40;
        cfg.t_ftimwin    = 0.5 * ones(1,length(cfg.foi));  % 7 cycles per time window
        cfg.toi          = -1:0.05:10;
     	total{index_sub,index_task} = ft_freqanalysis(cfg, data_eeg);
    end
end
clear data_eeg 
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
            1,  7,  NaN,NaN,NaN,    NaN;...  	%S18    OK
            1,  7,  6,  2,  3,      NaN;... 	%S19    OK
            1,  7,  NaN,NaN,NaN,    NaN];     	%S20    OK
%%
numTrials1 = zeros(20,2);
for index_sub = 1:length(subjectsFolders)
    load(fullfile(current_path,subjectsFolders{index_sub},strcat('S_',subjectsNames{index_sub},'_Binary_ICA_clean')));%% load tenStim_seg_AR
    a = [];
    for indlow = uselabel(index_sub,uselabel(index_sub,:)<4)
        a  = cat(1,a, find(data_eeg.trialinfo(:,end) == indlow));
    end
    numTrials1(index_sub,1) = length(a);
    a = [];
    for indhigh = uselabel(index_sub,uselabel(index_sub,:)>4)
        a  = cat(1,a, find(data_eeg.trialinfo(:,end) == indhigh));
    end
    numTrials1(index_sub,2) = length(a);
end
%% averaging and normalizing via bootstrap to minimized the effect of unbalanced trials for each condition per subject
nonphaselock_ave  = total;
for index_sub = 1 : length(subjectsFolders)
    disp(index_sub)
    numTrials = zeros(1,2);
    a = [];
    for indlow = uselabel(index_sub,uselabel(index_sub,:)<4)
        a  = cat(1,a, find(nonphaselock_ave{index_sub,1}.trialinfo(:,end) == indlow));
    end
    numTrials(1,1) = length(a);
    a = [];
    for indhigh = uselabel(index_sub,uselabel(index_sub,:)>4)
        a  = cat(1,a, find(nonphaselock_ave{index_sub,2}.trialinfo(:,end) == indhigh));
    end
    numTrials(1,2) = length(a);
    if numTrials(1,1) > numTrials(1,2)
        UnFlag = 1;
        faFlag = 0;
    elseif numTrials(1,1) < numTrials(1,2)
        UnFlag = 0;
        faFlag = 1;
    else
        UnFlag = 0;
        faFlag = 0;
    end
    for index_task = 1:2
        if index_task == 1
            if UnFlag == 1
                n = numTrials(1,1);
                m = numTrials(1,2);
                subpow = [];
                for surr = 1 : 100
                    cfg = [];
                    seltrials = [];
                    for indlow = uselabel(index_sub,uselabel(index_sub,:)<4)
                        seltrials  = cat(1,seltrials, find(nonphaselock_ave{index_sub,index_task}.trialinfo(:,end) == indlow));
                    end
                    cfg.trials = seltrials(randperm(n,m));
                    cfg.avgoverrpt  = 'yes';
                    aveTrial = ft_selectdata(cfg, nonphaselock_ave{index_sub,index_task});
                    subpow = cat(4,subpow,aveTrial.powspctrm); %Ch * Freq * Time
                end
                aveTrial.powspctrm = mean(subpow,4);
            else
                cfg = [];
                cfg.trials = [];
                for indlow = uselabel(index_sub,uselabel(index_sub,:)<4)
                    cfg.trials  = cat(1,cfg.trials, find(nonphaselock_ave{index_sub,index_task}.trialinfo(:,end) == indlow));
                end
                cfg.avgoverrpt  = 'yes';
                aveTrial = ft_selectdata(cfg, nonphaselock_ave{index_sub,index_task});
            end
        elseif index_task == 2
            if faFlag == 1
                n = numTrials(1,2);
                m = numTrials(1,1);
                subpow = [];
                for surr = 1 : 100
                    cfg = [];
                    seltrials = [];
                    for indhigh = uselabel(index_sub,uselabel(index_sub,:)>4)
                        seltrials  = cat(1,seltrials, find(nonphaselock_ave{index_sub,index_task}.trialinfo(:,end) == indhigh));
                    end
                    cfg.trials = seltrials(randperm(n,m));
                    cfg.avgoverrpt  = 'yes';
                    aveTrial = ft_selectdata(cfg, nonphaselock_ave{index_sub,index_task});
                    subpow = cat(4,subpow,aveTrial.powspctrm); %Ch * Freq * Time
                end
                aveTrial.powspctrm = mean(subpow,4);
            else
                cfg = [];
                cfg.trials = [];
                for indhigh = uselabel(index_sub,uselabel(index_sub,:)>4)
                    cfg.trials  = cat(1,cfg.trials, find(nonphaselock_ave{index_sub,index_task}.trialinfo(:,end) == indhigh));
                end
                cfg.avgoverrpt  = 'yes';
                aveTrial = ft_selectdata(cfg, nonphaselock_ave{index_sub,index_task});
            end
        else
            error('There is only two tasks')
        end
        nonphaselock_ave{index_sub,index_task} = aveTrial;
    end
end
save(fullfile(current_path,'totalv10_20SubFinal.mat'),'nonphaselock_ave','-v7.3');
%% baseline correction on the average
nonphaselock1  = nonphaselock_ave;
for index_sub = 1 : length(subjectsFolders)
    disp(index_sub)
    for index_task = 1:2
        cfg = [];
        cfg.baseline     = [-0.95 -0.25];%[-0.85 -0.45];  [-0.95 -0.25]
        cfg.baselinetype  = 'db'; %, 'relative', 'relchange', 'normchange', 'db' or 'zscore' (default = 'absolute')
        nonphaselock1{index_sub,index_task} = ft_freqbaseline(cfg, nonphaselock_ave{index_sub,index_task});
    end
end
%% statistic between two conditions
sel = 1:20;
[stat_nonphase] = statistic_analysis_TFN1(nonphaselock1, sel);
cfg                 = [];
cfg.parameter       = 'stat';
cfg.layout          = 'easycapM1.mat';
cfg.marker          =  'on';
cfg.maskparameter   = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle       = 'outline';% opacity  outline
cfg.maskalpha       = 0.05; %% 0.05 and  0.1
cfg.colorbar    	= 'yes';
cfg.colormap        = redblue(50);
cfg.colormap        = 'jet';
figure; ft_multiplotTFR(cfg, stat_nonphase);
%% careful with print
print('-dtiff','-r300',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TFR/TRF_Pzdiff'))
% print('-dtiff','-r250',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/SuppTFR/TRF_F4diff'))
%% Frequency Analysis of Alpha band (Fig.2A)
nonphaselock1_freq = nonphaselock1;
for index_sub = 1 : length(subjectsFolders)
    for index_task = 1 : 2
        cfg = [];
        cfg.latency     = [0.7 5];
        cfg.avgovertime = 'yes';
        cfg.frequency     = [8 12];
        cfg.avgoverfreq = 'yes';
        nonphaselock1_freq{index_sub,index_task} = ft_selectdata(cfg, nonphaselock1_freq{index_sub,index_task});
    end
end
[stat_nonphase_freq0] = statistic_analysis_TFN2(nonphaselock1_freq, 1:20);
cfg = [];
cfg.parameter               = 'stat';
cfg.layout                  = 'easycapM1.mat';
cfg.alpha                   = 0.05;
%cfg.zlim = [-4 4];
cfg.subplotsize             = [1 1];
cfg.highlightseries         = {'on', 'on', 'on'};
cfg.highlightsymbolseries   =['*', 'x', '+'];
cfg.highlightcolorpos       = [1 0 1];
cfg.highlightcolorneg       = [1 1 0];
cfg.highlightsizeseries     = [16 16 12 12 12];
cfg.marker                  =  'on';
cfg.colorbar                = 'no';%%East  no
cfg.colormap                = 'jet';
cfg.markersize              = 8;
cfg.comment                 = 'no'; %%no  xlim
cfg.commentpos              = 'layout';
cfg.highlightsymbol         = 'o';
cfg.highlightsize           = 10;
cfg.highlightfontsize       = 10;
ft_clusterplot(cfg, stat_nonphase_freq0);
figure;
stat_nonphase_freq0.stat = stat_nonphase_freq0.stat-1;
cfg.contournum              = 10;
ft_topoplotER(cfg, stat_nonphase_freq0)
% print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoAlpha'))
%% Frequency Analysis of low Beta band (Fig.2A)
nonphaselock1_freq = nonphaselock1;
for index_sub = 1 : length(subjectsFolders)
    for index_task = 1 : 2
        cfg = [];
        cfg.latency     = [0.7 5];
        cfg.avgovertime = 'yes';
        cfg.frequency     = [12 15.9];
        cfg.avgoverfreq = 'yes';
        nonphaselock1_freq{index_sub,index_task} = ft_selectdata(cfg, nonphaselock1_freq{index_sub,index_task});
    end
end
[stat_nonphase_freq] = statistic_analysis_TFN3(nonphaselock1_freq, 1:20);
cfg                         = [];
cfg.parameter               = 'stat';
cfg.layout                  = 'easycapM1.mat';
cfg.alpha                   = 0.05;
%cfg.zlim = [-4 4];
cfg.subplotsize             = [1 1];
cfg.highlightseries         = {'on', 'on', 'on'};
cfg.highlightsymbolseries   =['*', 'x', '+'];
cfg.highlightcolorpos       = [1 0 1];
cfg.highlightcolorneg       = [1 1 0];
cfg.highlightsizeseries     = [16 16 12 12 12];
cfg.marker                  =  'on';
cfg.colorbar                = 'no';%%East  no
cfg.colormap                = 'jet';
cfg.markersize              = 8;
cfg.comment                 = 'no'; %%no  xlim
cfg.commentpos              = 'layout';
cfg.highlightsymbol         = 'o';
cfg.highlightsize           = 10;
cfg.highlightfontsize       = 10;
ft_clusterplot(cfg, stat_nonphase_freq);
figure;
stat_nonphase_freq.stat = stat_nonphase_freq.stat-1.5;
cfg.zlim                    = [-1 1];
cfg.contournum              = 10;
ft_topoplotER(cfg, stat_nonphase_freq)
% print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoBeta'))
%% To find out that Unfamiliar has no effect but familiar has effect Alpha band(Fig.2A)
%%%%%%%%%%%%% Alpha %%%%%%%%%%%%%%
clear grandavg
for index_task = 1 : 2
    cfg = [];
    cfg.keepindividual = 'yes';
    grandavg{index_task} = ft_freqgrandaverage(cfg, nonphaselock1{1:20,index_task});
end
for index_task = 1 : 2
    cfg = [];
    cfg.latency         = [0.7 5];
    cfg.avgovertime     = 'yes';
    cfg.frequency       = [8 12];
    cfg.avgoverfreq     = 'yes';
    %cfg.trials         = 'all';
    %cfg.avgoverrpt     = 'yes';
    cfg.nanmean         = 'yes';
    grandavgS{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
    %grandavg{index_task}.powspctrm = squeeze(grandavg{index_task}.powspctrm)';
    cfg.latency         = [-0.95 -0.25];
    grandavgB{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
end
clear grandavgA
grandavgA{1} = grandavgS{1};grandavgA{1}.time = 0;
grandavgA{2} = grandavgB{1};grandavgA{2}.time = 0;
grandavgA{2}.powspctrm = grandavgB{1}.powspctrm;
[stat1] = statistic_analysis_TFN4(grandavgA,1:20);
print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoAlphaUnf'))

clear grandavgA
grandavgA{1} = grandavgS{2};grandavgA{1}.time = 0;
grandavgA{2} = grandavgB{2};grandavgA{2}.time = 0;
% grandavgA{2}.powspctrm = (grandavgB{1}.powspctrm+grandavgB{2}.powspctrm)/2;
grandavgA{2}.powspctrm = grandavgB{2}.powspctrm;
[stat2] = statistic_analysis_TFN4(grandavgA,1:20);
print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoAlphaFam'))

for index_sub = 1 : length(subjectsFolders)
for index_task = 1 : 2
    cfg = [];
    cfg.latency     = [0.7 5];
    cfg.avgovertime = 'yes';
    cfg.frequency     = [8 12];
    cfg.avgoverfreq = 'yes';
    %cfg.trials   = 'all';
    %cfg.avgoverrpt = 'yes';
    cfg.nanmean     = 'yes';
    grandavgSS{index_sub,index_task} = ft_selectdata(cfg, nonphaselock1{index_sub,index_task});
    %grandavg{index_task}.powspctrm = squeeze(grandavg{index_task}.powspctrm)';
    cfg.latency     = [-0.95 -0.25];
    grandavgBB{index_sub,index_task} = ft_selectdata(cfg, nonphaselock1{index_sub,index_task});
end
% grandavgBB{index_sub,1}.powspctrm = (grandavgBB{index_sub,1}.powspctrm+grandavgBB{index_sub,2}.powspctrm)/2;
end
effectSizeTFRpaired(grandavgSS(1:20,1),grandavgBB(1:20,1),[-1 10],[8 12],...
    {'C1','C3','CP1','CP3'},20)
effectSizeTFRpaired(grandavgSS(1:20,2),grandavgBB(1:20,2),[-1 10],[8 12],...
    {'C1','C3','CP1','CP3'},20)
%% To find out that Unfamiliar has no effect but familiar has effect LowBeta band(Fig.2A)
%%%%%%%%%%%%% Beta %%%%%%%%%%%%%%
clear grandavg
for index_task = 1 : 2
    cfg = [];
    cfg.keepindividual = 'yes';
    grandavg{index_task} = ft_freqgrandaverage(cfg, nonphaselock1{1:20,index_task});
end
for index_task = 1 : 2
    cfg = [];
    cfg.latency     = [0.7 5];
    cfg.avgovertime = 'yes';
    cfg.frequency     = [12 15.99];
    cfg.avgoverfreq = 'yes';
    %cfg.trials   = 'all';
    %cfg.avgoverrpt = 'yes';
    cfg.nanmean     = 'yes';
    grandavgS{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
    %grandavg{index_task}.powspctrm = squeeze(grandavg{index_task}.powspctrm)';
    cfg.latency     = [-0.95 -0.25];
    grandavgB{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
end
clear grandavgA
grandavgA{1} = grandavgS{1};grandavgA{1}.time = 0;
grandavgA{2} = grandavgB{1};grandavgA{2}.time = 0;
% grandavgA{2}.powspctrm = (grandavgB{1}.powspctrm+grandavgB{2}.powspctrm)/2;
grandavgA{2}.powspctrm = grandavgB{1}.powspctrm;
[stat1] = statistic_analysis_TFN5(grandavgA,1:20);
print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoBetaUnf'))

clear grandavgA
grandavgA{1} = grandavgS{2};grandavgA{1}.time = 0;
grandavgA{2} = grandavgB{2};grandavgA{2}.time = 0;
% grandavgA{2}.powspctrm = (grandavgB{1}.powspctrm+grandavgB{2}.powspctrm)/2;
grandavgA{2}.powspctrm = grandavgB{2}.powspctrm;
[stat2] = statistic_analysis_TFN5(grandavgA,1:20);
print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoBetaFam'))

for index_sub = 1 : length(subjectsFolders)
for index_task = 1 : 2
    cfg = [];
    cfg.latency     = [0.7 5];
    cfg.avgovertime = 'yes';
    cfg.frequency     = [12 15.99];
    cfg.avgoverfreq = 'yes';
    %cfg.trials   = 'all';
    %cfg.avgoverrpt = 'yes';
    cfg.nanmean     = 'yes';
    grandavgSS{index_sub,index_task} = ft_selectdata(cfg, nonphaselock1{index_sub,index_task});
    %grandavg{index_task}.powspctrm = squeeze(grandavg{index_task}.powspctrm)';
    cfg.latency     = [-0.95 -0.25];
    grandavgBB{index_sub,index_task} = ft_selectdata(cfg, nonphaselock1{index_sub,index_task});
end
grandavgBB{index_sub,1}.powspctrm = (grandavgBB{index_sub,1}.powspctrm+grandavgBB{index_sub,2}.powspctrm)/2;
grandavgBB{index_sub,2}.powspctrm = grandavgBB{index_sub,1}.powspctrm;
end
effectSizeTFRpaired(grandavgSS(1:20,1),grandavgBB(1:20,2),[-1 10],[0 16],...
    {'C3','CP3',},20)
effectSizeTFRpaired(grandavgSS(1:20,2),grandavgBB(1:20,2),[-1 10],[0 16],...
    {'C3','CP3'},20)
%% variation of frequencies (Alpha and low-Beta) over time (Fig.2b)
cfg = [];
cfg.parameter       = 'stat';
cfg.frequency       = [8 12];%% [8 12] [12 15.99]
cfg.avgoverfreq     = 'yes';
stat_nonphase_freq  = ft_selectdata(cfg, stat_nonphase);
% for the multiple plots also
Xlim = -1:1:5; %%  0:1:10; -1:1:5; 0:0.5:5;
for ind_Xlim = 1 : length(Xlim)-1
    cfg = [];
    cfg.parameter           = 'stat';
    cfg.xlim                = [Xlim(ind_Xlim) Xlim(ind_Xlim+1)];
    cfg.zlim                = [-1.5 1.5];
    cfg.marker              = 'on';
    cfg.colorbar            = 'no';%%East  no
    cfg.colormap            = 'jet';
    cfg.comment             = 'no'; %%no  xlim
    cfg.commentpos          = 'layout';
    cfg.layout              = 'easycapM1.mat';
    cfg.markersize          = 8;
    cfg.highlightsymbol     = 'o';
    cfg.highlightsize       = 10;
    cfg.highlightfontsize   = 10;
    %figure('units','normalized','outerposition',[0 0 1 1])
    ft_topoplotTFR(cfg,stat_nonphase_freq);
    print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoVarAlpha','_',num2str(ind_Xlim)))
%     print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/SuppTopo/topoVarBeta','_',num2str(ind_Xlim)))
    close all
end
%% ft_multiplotTFR for each condition and topoplot variation for each condition (Fig2 and Fig3)
clear grandavg
for index_task = 1 : 2
    cfg = [];
    cfg.keepindividual = 'no';
    grandavg{index_task} = ft_freqgrandaverage(cfg, nonphaselock1{1:20,index_task});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot condition 1
cfg = [];
cfg.parameter   = 'powspctrm';
cfg.layout      = 'easycapM1.mat';
cfg.channel     = {'EEG'};
cfg.xlim    	= [-0.7 5];%% [-0.7 5] and [-0.0 10]
cfg.ylim        = [3 40];%% [3 40] and [3 30]
cfg.colorbar	= 'yes';
cfg.colormap    = 'jet';
figure; ft_multiplotTFR(cfg, grandavg{1});
%%%%%% careful with print
print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TFR/TRF_PzUnf'))
% print('-dtiff','-r250',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/SuppTFR/TRF_F4Fam'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% variation of frequencies (Alpha) over time (Fig.2b)
cfg = [];
cfg.parameter   = 'powspctrm';
cfg.frequency	= [8 12];
cfg.avgoverfreq = 'yes';
stat_nonphase_freq = ft_selectdata(cfg, grandavg{2});%% grandavg{1}:unfamiliar grandavg{2}:Familiar
% for the multiple plots also
Xlim = -1:1:5; %%  0:1:10;  0:0.5:5;  Xlim = 5:1:10;
for ind_Xlim = 1 : length(Xlim)-1
	cfg = [];
    cfg.parameter       = 'powspctrm';
    cfg.xlim            = [Xlim(ind_Xlim) Xlim(ind_Xlim+1)];
    cfg.zlim            = [-1 1];
    cfg.marker          = 'on';
    cfg.colorbar        = 'no';%%East  no
    cfg.colormap        = 'jet';
    cfg.comment         = 'no'; %%no  xlim
    cfg.commentpos      = 'layout';
    cfg.layout          = 'easycapM1.mat';
    cfg.markersize          = 8;
    cfg.highlightsymbol     = 'o';
    cfg.highlightsize       = 10;
    cfg.highlightfontsize   = 10;
    %figure('units','normalized','outerposition',[0 0 1 1])
    ft_topoplotTFR(cfg,stat_nonphase_freq);
    print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoVarAlphaFam','_',num2str(ind_Xlim)))
    %print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/SuppTopo/topoVarAlphaUnf','_',num2str(ind_Xlim+5)))
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot condition 2
cfg = [];
cfg.parameter   = 'powspctrm';
cfg.layout      = 'easycapM1.mat';
cfg.channel   	= {'EEG'};
cfg.xlim     	= [-0.7 5];
cfg.ylim     	= [3 40];
cfg.colorbar  	= 'yes';
cfg.colormap    = 'jet';
figure; ft_multiplotTFR(cfg, grandavg{2});
%%%%%% careful with print
print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TFR/TRF_PzFam'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% variation of frequencies (Low-Beta) over time (Fig.2b)
cfg = [];
cfg.parameter   = 'powspctrm';
cfg.frequency	= [12 15.99];
cfg.avgoverfreq = 'yes';
stat_nonphase_freq = ft_selectdata(cfg, grandavg{1});
% for the multiple plots also
Xlim = -1:1:5; %%  0:1:10;  0:0.5:5;  Xlim = 5:1:10;
for ind_Xlim = 1 : length(Xlim)-1
	cfg = [];
    cfg.parameter       = 'powspctrm';
    cfg.xlim            = [Xlim(ind_Xlim) Xlim(ind_Xlim+1)];
    cfg.zlim            = [-1.25 1.25];
    cfg.marker          = 'on';
    cfg.colorbar        = 'no';%%East  no
    cfg.colormap        = 'jet';
    cfg.comment         = 'no'; %%no  xlim
    cfg.commentpos      = 'layout';
    cfg.layout          = 'easycapM1.mat';
    cfg.markersize          = 8;
    cfg.highlightsymbol     = 'o';
    cfg.highlightsize       = 10;
    cfg.highlightfontsize   = 10;
    %figure('units','normalized','outerposition',[0 0 1 1])
    ft_topoplotTFR(cfg,stat_nonphase_freq);
    print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TopoFreq/topoVarBetaUnf','_',num2str(ind_Xlim)))
    %print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/SuppTopo/topoVarBetaFam','_',num2str(ind_Xlim+5)))
    close all
end
%% Cohen's d calcuation
%%%% Alpha and beta for unfamiliar versus familiar
nonphaselock1  = nonphaselock_ave;
for index_sub = 1 : length(subjectsFolders)
    disp(index_sub)
    for index_task = 1:2
        cfg = [];
        cfg.baseline     = [-0.95 -0.25];%[-0.85 -0.45];
    	cfg.baselinetype  = 'db'; %, 'relative', 'relchange', 'normchange', 'db' or 'zscore' (default = 'absolute')
        nonphaselock1{index_sub,index_task} = ft_freqbaseline(cfg, nonphaselock_ave{index_sub,index_task});
    end
end
%%%% familiar vs unfamiliar
effectSizeTFRpaired(nonphaselock1(1:20,1),nonphaselock1(1:20,2),[0.8 4.8],[8 12],...
    {'Fz','FCz','F1','FC1','Pz','CP1'},20)
effectSizeTFRpaired(nonphaselock1(1:20,1),nonphaselock1(1:20,2),[0.8 4.8],[12 15.99],...
    {'FCz','Fz','FC1','F1','FC3','F3','AF3','FC5','F5'},20)
%% Frequency Analysis other freqs
clear stat_nonphase_freqOther
otherfreq = [4 8;16 22;22 32;32 50];
for indfreq = 1 : size(otherfreq,1)
    nonphaselock1_freq = nonphaselock1;
    for index_sub = 1 : length(subjectsFolders)
        for index_task = 1 : 2
            cfg = [];
            cfg.latency     = [0.5 5];
            cfg.avgovertime = 'yes';
            cfg.frequency     = otherfreq(indfreq,:);
            cfg.avgoverfreq = 'yes';
            nonphaselock1_freq{index_sub,index_task} = ft_selectdata(cfg, nonphaselock1_freq{index_sub,index_task});
        end
    end
    [stat_nonphase_freqOther{indfreq}] = statistic_analysis_TFN2(nonphaselock1_freq, 1:20);
end
%%
%%%% familiar vs unfamiliar
effectSizeTFRpaired(nonphaselock1(:,1),nonphaselock1(:,2),[0.25 5],[7.5 11.5],...
    {'F1','Fz','FCz','FC1','Pz','CP1'},15)
effectSizeTFRpaired(nonphaselock1(:,1),nonphaselock1(:,2),[0.25 5],[12 16],...
    {'FCz','Fz','FC1','F1','FC3','F3','AF3','FC5','F5','C5','F7','FT7','T7'},15)
%%%% familiar vs baseline
nonphaselockF = nonphaselock1;
for indSub = 1 : 15
    for indTask = 1 : 2
        if indTask == 1
            cfg = [];
            cfg.latency   = [-0.85 -0.25];
            cfg.avgovertime = 'yes';
            cfg.nanmean     = 'yes';
            [nonphaselockF{indSub,indTask}] = ft_selectdata(cfg, nonphaselockF{indSub,2});
            nonphaselockF{indSub,indTask}.time = 2.624999999999906;
        else
            cfg = [];
            cfg.latency   = [0.25 5];
            cfg.avgovertime = 'yes';
            cfg.nanmean     = 'yes';
            [nonphaselockF{indSub,indTask}] = ft_selectdata(cfg, nonphaselockF{indSub,indTask});
        end
    end
end 
effectSizeTFRpaired(nonphaselockF(:,1),nonphaselockF(:,2),[0.25 5],[7.5 11.5],...
    {'AF3', 'Fz','F1', 'FC1','F5','FT7','F7'},15)
effectSizeTFRpaired(nonphaselockF(:,1),nonphaselockF(:,2),[0.25 5],[12 16],...
    {'FCz', 'Fz', 'Cz', 'FC1', 'F1', 'C1', 'FC3', 'F3', 'C3', 'AF3', 'FC5', 'F5', 'CP1', 'CP3'},15)
%%%% unfamiliar vs baseline
nonphaselockU = nonphaselock1;
for indSub = 1 : 15
    for indTask = 2 : -1:1
        if indTask == 2
            cfg = [];
            cfg.latency   = [-1 -0.25];
            cfg.avgovertime = 'yes';
            cfg.nanmean     = 'yes';
            [nonphaselockU{indSub,indTask}] = ft_selectdata(cfg, nonphaselockU{indSub,1});
            nonphaselockU{indSub,indTask}.time = 2.624999999999906;
        else
            cfg = [];
            cfg.latency   = [0.25 5];
            cfg.avgovertime = 'yes';
            cfg.nanmean     = 'yes';
            [nonphaselockU{indSub,indTask}] = ft_selectdata(cfg, nonphaselockU{indSub,indTask});
        end
    end
end 
effectSizeTFRpaired(nonphaselockU(:,1),nonphaselockU(:,2),[0.25 5],[7.5 11.5],...
    {'F1','Fz','FCz','FC1','Pz','CP1'},15)
effectSizeTFRpaired(nonphaselockU(:,1),nonphaselockU(:,2),[0.25 5],[12 16],...
    {'FCz','Fz','FC1','F1','FC3','F3','AF3','FC5','F5','C5','F7','FT7','T7'},15)
%% Fieldtrip to EEGLAB for 3D pictures
ChanLoc = readlocs('Standard-10-20-Cap81.locs');
clear ChanLocnew
for indloc = 1 : length(grandavg{1,1}.label)
    for indtest = 1 : 81
        if strcmp(ChanLoc(indtest).labels, grandavg{1,1}.label{indloc})
            ChanLocnew(indloc) = ChanLoc(indtest);
        end
    end
end
%%% To find out that Unfamiliar has no effect but familiar has effect
for index_task = 1 : 2
    cfg = [];
    cfg.keepindividual = 'yes';
    grandavg{index_task} = ft_freqgrandaverage(cfg, nonphaselock1{:,index_task});
end
for index_task = 1 : 2
    cfg = [];
    cfg.keepindividual = 'yes';
    grandavg{index_task} = ft_freqgrandaverage(cfg, nonphaselock1{:,index_task});
end
clear EEGbeta EEGalpha restBeta restAlpha
for index_task = 1 : 2
    cfg = [];
    cfg.latency     = [0 5];
    cfg.avgovertime = 'yes';
    cfg.frequency     = [12 16];
    cfg.avgoverfreq = 'yes';
    cfg.trials   = 'all';
    cfg.avgoverrpt = 'yes';
    cfg.nanmean     = 'yes';
    EEGbeta{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
    EEGbeta{index_task}.powspctrm = squeeze(EEGbeta{index_task}.powspctrm)';
    cfg.latency     = [-1 0];
    restBeta{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
    restBeta{index_task}.powspctrm = squeeze(restBeta{index_task}.powspctrm)';
    cfg.latency     = [0 5];
    cfg.frequency     = [8 12];
    EEGalpha{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
    EEGalpha{index_task}.powspctrm = squeeze(EEGalpha{index_task}.powspctrm)';
    cfg.latency     = [-1 0];
    restAlpha{index_task} = ft_selectdata(cfg, grandavg{:,index_task});
    restAlpha{index_task}.powspctrm = squeeze(restAlpha{index_task}.powspctrm)';
end
restBeta{1}.powspctrm = (restBeta{1}.powspctrm+restBeta{2}.powspctrm)./2;
DataField = EEGbeta{1};
clear DataF 
DataF. trial{1} = DataField.powspctrm;
DataF.fsample = 1000;
DataF.time{1} = DataField.time;
[EEG_Unf] = fieldtrip2eeglab(DataF);
EEG_Unf.chanlocs = ChanLocnew;
EEGOUT = pop_headplot( EEG_Unf, 1);% [], 'Alpha band [7.5-11.5 Hz]');  Beta band [8-12 Hz]


%% pooling regression Fig3
ttime = [-1 0; 0 1;1 2;2 3;3 4;4 5;5 6;6 7;7 8;8 9;9 10];
cchan = {{'Pz';'CP1'},{'F1';'Fz'},{'AF3';'F3';'F5';'FC3';'FC5'},{'AF4';'F4';'F6';'FC4';'FC6'},... %% alpha
         {'Pz';'CP1'},{'F1';'Fz'},{'AF3';'F3';'F5';'FC3';'FC5'},{'AF4';'F4';'F6';'FC4';'FC6'}};%% low beta
ffreq = [8 12; 12 15.99];
datareg = zeros(2,11,8,19); % two cond * 11 time point * 8 regions * 19 subjects
count = 1;
for index_sub = 1:20
    disp(index_sub)
    for index_task = 1:2
        for indchan = 1:4
            for indtime = 1:11
                cfg = [];
                cfg.nanmean     = 'yes';
                cfg.avgoverfreq = 'yes';
                cfg.avgovertime = 'yes';
                cfg.avgoverchan = 'yes';
                cfg.channel     = cchan{indchan}; 
                cfg.frequency   = ffreq(1,:);
                cfg.latency     = ttime(indtime,:);
                X = ft_selectdata(cfg, nonphaselock1{index_sub,index_task});
                datareg(index_task,indtime,indchan,count) = X.powspctrm;
            end
        end
        for indchan = 4:8
            for indtime = 1:11
                cfg = [];
                cfg.nanmean     = 'yes';
                cfg.avgoverfreq = 'yes';
                cfg.avgovertime = 'yes';
                cfg.avgoverchan = 'yes';
                cfg.channel     = cchan{indchan}; 
                cfg.frequency   = ffreq(2,:);
                cfg.latency     = ttime(indtime,:);
                X = ft_selectdata(cfg, nonphaselock1{index_sub,index_task});
                datareg(index_task,indtime,indchan,index_sub) = X.powspctrm;
            end
        end
    end
    count = count + 1;
end
acolor1 = {[0.75 1 0.75];[1 0.75 0.75];[0.75 0.75 1]};
acolor2 = {[0 1 0];[1 0 0];[0 0 1]};
%% Continue pooling regression Fig3
ttitle = {{'Reg_Alpha_Pz'};...
            {'Reg_Alpha_Fz'};...
            {'Reg_Alpha_F5'};...
            {'Reg_Alpha_F6'};...
            {'Reg_Beta_Pz'};...
            {'Reg_Beta_Fz'};...
            {'Reg_Beta_F5'};...
            {'Reg_Beta_F6'}};
for creg = 1:1:8
    figure('units','normalized','outerposition',[0 0 1 1])
    errorbar(-0.5:1:4.5,...
        squeeze(mean(datareg(1,1:6,creg,:),4)),...
        squeeze(0.5*(std(datareg(1,1:6,creg,:),[],4))),...
        '-*','MarkerSize',15,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1],'LineWidth', 5,'Color',[0,0,1])
    hold on
    errorbar(-0.5:1:4.5,...
        squeeze(mean(datareg(2,1:6,creg,:),4)),...
        squeeze(0.5*(std(datareg(2,1:6,creg,:),[],4))),...
        '-s','MarkerSize',15,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0],'LineWidth', 5,'Color',[1,0,0])
    hold on
    ylim([-2, 1.5]);
    yticks(-2:0.5:1.5);
    xlim([-.75, 4.75]);
    xticks(-0.5:1:4.5);
    xticklabels({'-0.5','0.5', '1.5', '2.5','3.5', '4.5'});
    grid; 
    xlabel('Time (s)')
    ylabel('Average of power')
    legend('Unfamiliar','Familiar','Location','Northeast','NumColumns',2,'FontSize',30)
    %title(ttitle{creg})
    set(gca,'FontSize',30)
    print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/TFR/',ttitle{creg}{1}))
end
%% pooling regression per subject
for index_sub = 1 : length(subjectsFolders)
    figure;
    plot(squeeze(datareg(1,:,6,index_sub)))
    hold on
    plot(squeeze(datareg(2,:,6,index_sub)))
end
%% pooling after 5 seconds
ttitle = {{'Reg_Alpha_Pz'};...
            {'Reg_Alpha_Fz'};...
            {'Reg_Alpha_F5'};...
            {'Reg_Alpha_F6'};...
            {'Reg_Beta_Pz'};...
            {'Reg_Beta_Fz'};...
            {'Reg_Beta_F5'};...
            {'Reg_Beta_F6'}};
for creg = 1:1:8
    figure('units','normalized','outerposition',[0 0 1 1])
    errorbar(-0.5:1:9.5,...
        squeeze(mean(datareg(1,1:11,creg,:),4)),...
        squeeze(0.5*(std(datareg(1,1:11,creg,:),[],4))),...
        '-*','MarkerSize',15,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1],'LineWidth', 5,'Color',[0,0,1])
    hold on
    errorbar(-0.5:1:9.5,...
        squeeze(mean(datareg(2,1:11,creg,:),4)),...
        squeeze(0.5*(std(datareg(2,1:11,creg,:),[],4))),...
        '-s','MarkerSize',15,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0],'LineWidth', 5,'Color',[1,0,0])
    hold on
    ylim([-2, 1.5]);
    yticks(-2:0.5:1.5);
    xlim([-0.75, 9.75]);
    xticks(-0.5:1:9.5);
    xticklabels({'-0.5','0.5', '1.5', '2.5','3.5', '4.5','5.5','6.5', '7.5', '8.5','9.5'});
    grid; 
    xlabel('Time (s)')
    ylabel('Average of power')
    legend('Unfamiliar','Familiar','Location','Northeast','NumColumns',2,'FontSize',30)
    %title(ttitle{creg})
    set(gca,'FontSize',30)
    print('-dtiff','-r500',strcat('F:/EEGdata/ClassicalSongParadigm/Results20Sub/SuppTFR/',ttitle{creg}{1}))
end
%% Figure 1
indsub = [1:1:20];
indcond = [1 2];
for subind = indsub
    for condind = indcond
f4 = figure('units','normalized','outerposition',[0 0 1 1]); 
y = 1:length(nonphaselock_ave{subind,condind}.label);
x = nonphaselock_ave{subind,condind}.freq;
z = nonphaselock_ave{subind,condind}.time(1:122);% 5 seconds
v = nonphaselock_ave{subind,condind}.powspctrm(:,:,1:122);
[x,y,z] = meshgrid(x,y,z);
slice(x,y,z,v,...
    [nonphaselock_ave{subind,condind}.freq(1) nonphaselock_ave{subind,condind}.freq(end)],...
    [1 length(nonphaselock_ave{subind,condind}.label)],...
    [nonphaselock_ave{subind,condind}.time(1) nonphaselock_ave{subind,condind}.time(122)])
xlabel('Frequency (Hz)')
xlim([nonphaselock_ave{subind,condind}.freq(1) nonphaselock_ave{subind,condind}.freq(end)])
ylabel('Channel')
ylim([1 51])
zlabel('Time (s)')
zlim([nonphaselock_ave{subind,condind}.time(1) nonphaselock_ave{subind,condind}.time(122)])
tickLocations = [-0.9,1,3,5]; % change to whatever you want
tickLabels    = {'-0.9','1','3','5'}; % change to whatever you want
set(gca,'zTick',tickLocations,'zTickLabel',tickLabels)
tickLocations = [8,18,28,38]; % change to whatever you want
tickLabels    = {'C3','P8','F8','TP7'}; % change to whatever you want
set(gca,'yTick',tickLocations,'yTickLabel',tickLabels)
colormap('jet')
ax = gca; 
ax.FontSize = 30;
print('-dtiff','-r150',strcat('fig_',num2str(subind),'_',num2str(condind)))
close all
    end
end
%%%%%%%%%%%%%%%%%%%%%
indsub = [1:1:20];
indcond = [1 2];
for subind = indsub
    for condind = indcond
f4 = figure('units','normalized','outerposition',[0 0 1 1]); 
y = 1:length(nonphaselock1{subind,condind}.label);
x = nonphaselock1{subind,condind}.freq;
z = nonphaselock1{subind,condind}.time(1:122);% 5 seconds
v = nonphaselock1{subind,condind}.powspctrm(:,:,1:122);
[x,y,z] = meshgrid(x,y,z);
slice(x,y,z,v,...
    [nonphaselock1{subind,condind}.freq(1) nonphaselock1{subind,condind}.freq(end)],...
    [1 length(nonphaselock1{subind,condind}.label)],...
    [nonphaselock1{subind,condind}.time(1) nonphaselock1{subind,condind}.time(122)])
xlabel('Frequency (Hz)')
xlim([nonphaselock1{subind,condind}.freq(1) nonphaselock1{subind,condind}.freq(end)])
ylabel('Channel')
ylim([1 51])
zlabel('Time (s)')
zlim([nonphaselock1{subind,condind}.time(1) nonphaselock1{subind,condind}.time(122)])
tickLocations = [-0.9,1,3,5]; % change to whatever you want
tickLabels    = {'-0.9','1','3','5'}; % change to whatever you want
set(gca,'zTick',tickLocations,'zTickLabel',tickLabels)
tickLocations = [8,18,28,38]; % change to whatever you want
tickLabels    = {'C3','P8','F8','TP7'}; % change to whatever you want
set(gca,'yTick',tickLocations,'yTickLabel',tickLabels)
colormap('jet')
caxis([-3 3]);
ax = gca; 
ax.FontSize = 30;
print('-dtiff','-r150',strcat('figNorm_',num2str(subind),'_',num2str(condind)))
close all
    end
end
%%%%%%%%%%%
freqAll = [4 8;8 12; 12 16;16 22;22 32;32 50];
for indf = 1 : 6
    nonphaselock1_freq = nonphaselock1;
    for index_sub = 3
        for index_task = 1 : 1
            cfg = [];
            cfg.latency     = [-0.75 -0.25];
            cfg.avgovertime = 'yes';
            cfg.frequency     = [freqAll(indf,1) freqAll(indf,2)];
            cfg.avgoverfreq = 'yes';
            nonphaselock2 = ft_selectdata(cfg, nonphaselock1_freq{index_sub,index_task});
            nonphaselock2.time=3;
            cfg = [];
            cfg.xlim = [1 5];
            cfg.ylim = [freqAll(indf,1) freqAll(indf,2)];
%             cfg.zlim = [-2 +2];
            cfg.layout = 'easycapM1.mat';
            cfg.comment = 'no';
            f4 = figure('units','normalized','outerposition',[0 0 1 1]); 
            ft_topoplotTFR(cfg,nonphaselock2); 
            colormap('jet')
            print('-dtiff','-r150',strcat('figtopoBase_',num2str(index_sub),'_',num2str(indf),'_',num2str(index_task)))
            close all
        end
    end
end