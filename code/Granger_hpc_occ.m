clearvars
tic
%
global ft_default
ft_default.spmversion = 'spm12';
ft_defaults
ft_hastoolbox('brewermap', 1);
region = 'anterior'; baseline = 0; resp_lock = 0;
nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];
hpc_nsubs = [1,3,4,7,9,10,12,13];

millis = linspace(-1,1,1001)*1000;
iszurich = logical([zeros(1,13),ones(1,3)]);

% create stp and patient_data structs
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);
p=1;
for subI = 1:numel(nsubs)
    if ~ismember(subI,hpc_nsubs)
        continue
    end
    fprintf(['Getting clean trials for Patient ',num2str(nsubs(subI)), '\n'])
    foldn = sprintf('Information/Patient%d+',nsubs(subI));
    cd(foldn)
    cd(sprintf('hpc_%s',patient_data(subI).hpc_axis))
    hpc_trials{p} = load('clean_trials_bipolar.mat');
    cd ../../../
    hpc_cleandata{p} = preproc_filtering(subI,nsubs,stp,patient_data,hpc_trials{p}.clean.trl,85);
    p=p+1;
end

region = 'occ'; baseline = 0; resp_lock = 0;
nsubs = [3,8,9,16,25,31,36,37];
iszurich = logical([zeros(1,8)]);

% create stp and patient_data structs
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);
for subI = 1:numel(nsubs)
    foldn = sprintf('Information/Patient%d+',nsubs(subI));
    cd(foldn)
    cd(sprintf('%s',patient_data(subI).region))
    occ_trials{subI} = load('clean_trials_bipolar.mat');
    cd ../../../
    occ_cleandata{subI} = preproc_filtering(subI,nsubs,stp,patient_data,occ_trials{subI}.clean.trl,85);
end

%%
for subI = 1:numel(nsubs)
    [keeptrials,ihpc,iocc] = intersect(hpc_trials{subI}.clean.trl,occ_trials{subI}.clean.trl);
    
    cfg=[];
    cfg.trials = ihpc;
    hpc_cleandata_sub{subI} = ft_selectdata(cfg,hpc_cleandata{subI});
    cfg=[];
    %cfg.avgoverchan = 'yes';
    cfg.channel = hpc_cleandata_sub{subI}.label{1};
    hpc_cleandata_avg{subI} = ft_selectdata(cfg,hpc_cleandata_sub{subI});
    
    cfg=[];
    cfg.trials = iocc;
    occ_cleandata_sub{subI} = ft_selectdata(cfg,occ_cleandata{subI});
    cfg=[];
    %         cfg.avgoverchan = 'yes';
    cfg.channel = occ_cleandata_sub{subI}.label{1};
    occ_cleandata_avg{subI} = ft_selectdata(cfg,occ_cleandata_sub{subI});
    
    cleandata = ft_appenddata([], hpc_cleandata_avg{subI},occ_cleandata_avg{subI});
    
    
%     % dft filter
%     cfg = [];
%     cfg.dftfilter ='yes';
%     cleandata = ft_preprocessing(cfg,cleandata);
    
    % resample
    cfg            = [];
    cfg.resamplefs = 250;
    cfg.detrend    = 'yes';
    cleandata = ft_resampledata(cfg,cleandata);
    cfg=[];
    cfg.latency = [-1 0];
    cleandata=ft_selectdata(cfg,cleandata);
        timepoint = size(cleandata.trial{1},2);

    info=hpc_trials{subI}.clean.info(ihpc,:);
    
    median_ent = median(info(:,1));
    small_idx=info(:,1) < median_ent;
    large_idx=~small_idx;
    
    cfg=[];
    cfg.trials = small_idx;
    small = ft_selectdata(cfg,cleandata);
    cfg=[];
    cfg.trials = large_idx;
    large = ft_selectdata(cfg,cleandata);
        % lat
        cfg=[];
        cfg.latency = [-1 -0.5];
        small=ft_selectdata(cfg,small);
        large=ft_selectdata(cfg,large);
    
    
        cfg            = [];
        cfg.output     = 'fourier';
        cfg.method     = 'mtmfft';
        cfg.taper      = 'dpss';
        cfg.foilim     = [0 80];
        cfg.pad        = 'nextpow2';
        cfg.tapsmofrq  = 10;
        small_ph       = ft_freqanalysis(cfg, small);
        large_ph       = ft_freqanalysis(cfg, large);
        cfg            = [];
        cfg.method = 'granger';
        granger_large{subI} = ft_connectivityanalysis(cfg, large_ph);
        granger_small{subI} = ft_connectivityanalysis(cfg, small_ph);

        for freqI = 1:size(granger_large{subI}.freq,2)
            large_hpc_to_occ{subI}.powspctrm(1,freqI) = log(granger_large{subI}.grangerspctrm(1,2,freqI));
            large_occ_to_hpc{subI}.powspctrm(1,freqI) = log(granger_large{subI}.grangerspctrm(2,1,freqI));
            small_hpc_to_occ{subI}.powspctrm(1,freqI) = log(granger_small{subI}.grangerspctrm(1,2,freqI));
            small_occ_to_hpc{subI}.powspctrm(1,freqI) = log(granger_small{subI}.grangerspctrm(2,1,freqI));
        end

    large_hpc_to_occ{subI}.freq = granger_large{subI}.freq;
    large_hpc_to_occ{subI}.dimord = 'chan_freq';
    large_hpc_to_occ{subI}.label = {'hpc-occ'};
    large_occ_to_hpc{subI}.freq = granger_large{subI}.freq;
    large_occ_to_hpc{subI}.dimord = 'chan_freq';
    large_occ_to_hpc{subI}.label = {'occ-hpc'};
    small_hpc_to_occ{subI}.freq = granger_large{subI}.freq;
    small_hpc_to_occ{subI}.dimord = 'chan_freq';
    small_hpc_to_occ{subI}.label = {'hpc-occ'};
    small_occ_to_hpc{subI}.freq = granger_large{subI}.freq;
    small_occ_to_hpc{subI}.dimord = 'chan_freq';
    small_occ_to_hpc{subI}.label = {'occ-hpc'};
    clear small large small_idx large_idx info keeptrials ihpc iocc
end

%%
cfg=[];
cfg.parameter='powspctrm';
cfg.keepindividual='yes';
avg_large_hpc_to_occ = ft_freqgrandaverage(cfg,large_hpc_to_occ{:});
avg_small_hpc_to_occ = ft_freqgrandaverage(cfg,small_hpc_to_occ{:});
avg_large_occ_to_hpc = ft_freqgrandaverage(cfg,large_occ_to_hpc{:});
avg_small_occ_to_hpc = ft_freqgrandaverage(cfg,small_occ_to_hpc{:});

cfg = [];
cfg.parameter = 'powspctrm';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';%'diff';
cfg.correctm         = 'cluster';
cfg.tail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 'all';
cfg.clustertail      = 0;
cfg.clusteralpha     = 0.05;
cfg.minnbchan        = 0;
cfg.correcttail      = 'alpha';
cfg.neighbours       = [];
cfg.ivar = 2;
cfg.uvar = 1;
% design matrices
clear design;
design(1,:) = [1:size(hpc_trials,2), 1:size(hpc_trials,2)];
design(2,:) = [ones(1,size(hpc_trials,2)), ones(1,size(hpc_trials,2)) * 2];
cfg.design = design;

gc_hpc2occ_stats = ft_freqstatistics(cfg,avg_large_hpc_to_occ,avg_small_hpc_to_occ);
gc_occ2hpc_stats = ft_freqstatistics(cfg,avg_large_occ_to_hpc,avg_small_occ_to_hpc);

cfg=[];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
% cfg.maskstyle = 'box';
figure();ft_singleplotER(cfg,gc_occ2hpc_stats);ylabel('GC T value');xlabel('Frequency (Hz)');
figure();ft_singleplotER(cfg,gc_hpc2occ_stats);ylabel('GC T value');xlabel('Frequency (Hz)');
figure();ft_singleplotER(cfg,directionality);ylabel('GC T value');xlabel('Frequency (Hz)');


figure('position',[10 10 500 400]);% 
bounds = nanstd(squeeze(avg_large_hpc_to_occ.powspctrm))/sqrt(numel(nsubs));
boundedline(avg_large_hpc_to_occ.freq, nanmean(squeeze(avg_large_hpc_to_occ.powspctrm)), bounds, 'alpha', 'transparency', 0.1,'cmap',[0.8500 0.3250 0.0980]);
hold on;
bounds = nanstd(squeeze(avg_small_hpc_to_occ.powspctrm))/sqrt(numel(nsubs));
boundedline(avg_large_hpc_to_occ.freq, nanmean(squeeze(avg_small_hpc_to_occ.powspctrm)), bounds, 'alpha', 'transparency', 0.1,'cmap',[0 0.4470 0.7410]);
legend({'','High entropy','','Low entropy'});ylabel('log(GC)');xlabel('Frequency')
xlim([1 80])

figure('position',[10 10 500 400]);
bounds = nanstd(squeeze(avg_large_occ_to_hpc.powspctrm))/sqrt(numel(nsubs));
boundedline(avg_large_hpc_to_occ.freq, nanmean(squeeze(avg_large_occ_to_hpc.powspctrm)), bounds, 'alpha', 'transparency', 0.1,'cmap',[0.8500 0.3250 0.0980]);
hold on;
bounds = nanstd(squeeze(avg_small_occ_to_hpc.powspctrm))/sqrt(numel(nsubs));
boundedline(avg_large_hpc_to_occ.freq, nanmean(squeeze(avg_small_occ_to_hpc.powspctrm)), bounds, 'alpha', 'transparency', 0.1,'cmap',[0 0.4470 0.7410]);
legend({'','High entropy','','Low entropy'});ylabel('log(GC)');xlabel('Frequency')
xlim([1 80]); hold off;
