%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% binary_ICA_clean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add path Fieldtrip
clear 
path(pathdef);
clc
close all

matlabrc
current_path = 'J:/Alireza Malekmohammadi/PhD_Projects/01_MusicFamiliarity/';
addpath 'F:/EEGdata/ClassicalSongParadigm/fieldtrip-20191008'
ft_defaults
%% add Subject folders and names 
subjectsRawFolders =  {'sbj_2019_06_05_Sebastian',...
                    'sbj_2019_06_17_Florian',...
                    'sbj_2019_06_18_Augosto',...
                    'sbj_2019_07_02_Lukas',...
                    'sbj_2019_07_09_Sai',...
                    'sbj_2019_07_27_Mostafa',...
                    'sbj_2019_07_29_George',...
                    'sbj_2019_07_30_Farnam',...
                    'sbj_2019_07_30_Poorya',...
                    'sbj_2019_08_04_Behnam',...
                    'sbj_2019_08_04_Hadi',...
                    'sbj_2019_08_05_Ali',...
                    'sbj_2019_08_07_Babak',...
                    'sbj_2019_08_08_Mahdi',...
                    'sbj_2019_09_24_Julian',...
                    '2022.01.13.Mohsen',...
                    '2022.01.14.Foad',...
                    '2022.01.16.Petar',...
                    '2022.01.18.Ehsan',...
                    '2022.01.27.Arturo'};
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
subjectsNames =    {'sebastian1',...
                    'FlorianS1',...
                    'augosto',...
                    'Lukas1',...
                    'Sai1',...
                    'Mostafa1',...
                    'George1',...%%%%attention  George1_1
                    'Farnam1',...
                    'Poorya1',...
                    'Behnam1',...
                    'Hadi1',...
                    'Ali1',...
                    'Babak1',...
                    'Mahdi1',...
                    'Julian1',...
                    'mohsen',...
                    'Foad',...
                    'petar',...
                    'ehsan',...
                    'arturo'};
sbjNames =    {'Sebastian',...
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
removingBADsubjects = [];
subjectsRawFolders(removingBADsubjects) = [];
subjectsFolders(removingBADsubjects) = [];
subjectsNames(removingBADsubjects) = [];
sbjNames(removingBADsubjects) = [];
%% load data, trails,
index_sub = 16;
cfg = [];
cfg.dataset = fullfile(current_path,subjectsRawFolders{index_sub},strcat(subjectsNames{index_sub},'.vhdr'));
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = 'S  1';
cfg.trialdef.prestim = 2.5;
cfg.trialdef.poststim = 12.5;
cfg = ft_definetrial(cfg);
cfg.continuous = 'yes';
cfg.lpfilter      = 'no'; %''no' or 'yes'  lowpass filter (default = 'no')
cfg.hpfilter      = 'no'; %'no' or 'yes'  highpass filter (default = 'no')
cfg.bpfilter      = 'yes'; %'no' or 'yes'  bandpass filter (default = 'no')
cfg.bsfilter      = 'no';%'no' or 'yes'  bandstop filter (default = 'no')
cfg.dftfilter     = 'yes'; %'no' or 'yes'  line noise removal using discrete fourier transform (default = 'no')
cfg.medianfilter  = 'no'; %'no' or 'yes'  jump preserving median filter (default = 'no')
cfg.lpfreq        = []; %lowpass  frequency in Hz
cfg.hpfreq        = []; %highpass frequency in Hz
cfg.bpfreq        = [0.5 100]; %bandpass frequency range, specified as [lowFreq highFreq] in Hz
cfg.bsfreq        = [];%bandstop frequency range, specified as [low high] in Hz (or as Nx2 matrix for notch filter)
cfg.dftfreq       = [50 100 ]; %line noise frequencies in Hz for DFT filter (default = [50 100 150])
cfg.lpfilttype    = []; %digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
cfg.hpfilttype    = [];%digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
cfg.bpfilttype    = 'but';%digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
cfg.bsfilttype    = [];%digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
cfg.lpfiltdir     = [];%filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
cfg.hpfiltdir     = [];%filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
cfg.bpfiltdir     = 'twopass';%filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
cfg.bsfiltdir     = [];%filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
cfg.plotfiltresp  = 'yes';%'no' or 'yes', plot filter responses (firws, default = 'no')
cfg.usefftfilt    = 'no';%'no' or 'yes', use fftfilt instead of filter (firws, default = 'no')
cfg.demean        = 'yes';%'no' or 'yes', whether to apply baseline correction (default = 'no')
cfg.detrend       = 'no';%'no' or 'yes', remove linear trend from the data (done per trial) (default = 'no')
cfg.polyremoval   = 'no';%'no' or 'yes', remove higher order trend from the data (done per trial) (default = 'no')
cfg.polyorder     = 0; %polynome order for poly trend removal (default = 2; note that all lower-order trends will also be removed when using cfg.polyremoval)
cfg.baselinewindow = [-1 -0.01];
if index_sub < 4
cfg.reref         = 'yes';
cf.channel        = 'all';
cfg.refchannel    = {'M1', 'TP9'};
cfg.implicitref   = 'M1';
end

data = ft_preprocessing(cfg);

if index_sub < 4
cfg = [];
cfg.channel = {'all', '-M1', '-TP9'};

[data] = ft_selectdata(cfg, data);
end
%% trials
data.trialinfo = xlsread(fullfile('F:/EEGdata/ClassicalSongParadigm/',subjectsFolders{index_sub},strcat(sbjNames{index_sub},'_docB')));
%% chnage label 
data.label{find(strcmp(data.label, 'EOGb'), 1)}='Fpz';
data.label{find(strcmp(data.label, 'EOGL'), 1)}='AF7';
data.label{find(strcmp(data.label, 'EOGR'), 1)}='AF8';
%% remove EOG channels
cfg = [];
cfg.channel = {'all', '-AF7', '-Fpz', '-AF8'};

[data] = ft_selectdata(cfg, data);
%% extract the time
cfg = [];
cfg.latency = [-1 11.499];
data = ft_selectdata(cfg, data);
%% ICA ECG, EOG removal
cfg = [];
cfg.method = 'sobi';
cfg.channel = 'all';
cfg.trials  = [1:85];

[comp] = ft_componentanalysis(cfg, data);
%% ploting components
cfg = [];
cfg.component = 26:50;
cfg.layout = 'easycapM1.mat';
cfg.comment = 'no';
cfg.colormap = 'jet';

ft_topoplotIC(cfg, comp);
%% show the component
cfg =[];
cfg.layout = 'easycapM1.mat';
cfg.viewmode = 'component';
cfg.channelcolormap = 'jet';

ft_databrowser(cfg, comp);
%% show raw data
cfg =[];
cfg.viewmode = 'vertical';
cfg.layout = 'easycapM1.mat';

ft_databrowser(cfg, data)
%% reject component
cfg =[];
cfg.component = [1 2 3 7 9 13 14 17 21 23 33];

[data1] = ft_rejectcomponent(cfg, comp, data);
%% show the data again
cfg =[];
cfg.viewmode = 'vertical';
cfg.layout = 'easycapM1.mat';

ft_databrowser(cfg, data1)
%% remove EOG channels
% cfg = [];
% cfg.channel = {'all', '-AF7', '-Fpz', '-AF8'};
% 
% [data1] = ft_selectdata(cfg, data1);
%% Trail view
cfg          = [];
cfg.method   = 'trial';
cfg.alim     = 100;
cfg.eegscale = 1;

ft_rejectvisual(cfg,data1);
%% reject other trial/channels
cfg          = [];
cfg.method   = 'summary';

data_eeg        = ft_rejectvisual(cfg,data1);
%% save data
save(fullfile(subjectsFolders{index_sub},strcat('S','_',sbjNames{index_sub},'_Binary_ICA_clean.mat')),'data_eeg')