%% load fus
clearvars
global ft_default
ft_default.spmversion = 'spm12';
ft_defaults
% get ts
region = 'fus'; baseline = 0; resp_lock = 0;
nsubs = [3,8,13,15,31,32,36];
iszurich = logical([zeros(1,7)]);
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);
for subI = 1:numel(nsubs)
    fprintf(['Getting clean trials for Patient ',num2str(nsubs(subI)), '\n'])
    if stp.zurich(subI) == 1
        foldn = sprintf('Information/P%dz',nsubs(subI));
    else
        foldn = sprintf('Information/Patient%d+',nsubs(subI));
    end
    cd(foldn)
    cd(sprintf('%s',patient_data(subI).region))
    load('clean_trials_bipolar.mat','clean')
    %this now has clean_trials:
    sub(subI).clean_trials = clean.info;
    cd ../../../
    stp.poststim = 1.4;
    stp.prestim = 1.8;
    timest={-stp.prestim:0.002:stp.poststim};
    
    cleandata = preproc(subI,nsubs,stp,patient_data,clean.trl);
    %timest = cleandata.time(1,1);
    cfg         = [];
    cfg.time    = repmat({timest{1}(1:end-1)},1,size(cleandata.trial,2));
    cleandata   = ft_resampledata(cfg,cleandata);
    
    cfg             = [];
    cfg.demean      = 'yes';
    fus_demeaned{subI} = ft_preprocessing(cfg,cleandata);
    cfg             = [];
    cfg.avgoverchan      = 'yes';
    fus_demeaned{subI} = ft_selectdata(cfg,fus_demeaned{subI});
    
    fus_demeaned{subI}.preds = clean.info;
    clear cleandata
    
end
fus = sub;

%% load hpc ripple and clean data
region = 'anterior'; baseline = 0; resp_lock = 0;
nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];
hpc_nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];
millis = linspace(-1,1,1001)*1000;
iszurich = logical([zeros(1,13),ones(1,3)]);

% create stp and patient_data structs
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);
hpfilt = 200;
ripdur = 25;
fname = sprintf('HPCRipples/HPCAnterior_vaz_hpf%d_%dms_%dsubjs_Jan22_reclean.mat',hpfilt, ripdur, numel(nsubs));
load(fname)
hpc = sub;

for subI = 1:numel(nsubs)
    if nsubs(subI) == 22 
        continue
    end
    k=1;
    fprintf(['Getting clean trials for Patient ',num2str(nsubs(subI)), '\n'])
    if stp.zurich(subI) == 1
        foldn = sprintf('Information/P%dz',nsubs(subI));
    else
        foldn = sprintf('Information/Patient%d+',nsubs(subI));
    end
    cd(foldn)
    cd(sprintf('hpc_%s',patient_data(subI).hpc_axis))
    load('clean_trials_bipolar.mat','clean')
    hpc(subI).clean_trials = clean.info;
    cd ../../../
end
%% subset hpc patients for those with fus contacts and get corresponding trls
hpc_fus=[1,3,5,6,10,11,12];
if sum(~ismember(hpc_nsubs(hpc_fus),nsubs)) ~= 0
    error('wrong patients subsetting!')
end
hpc_fus_sub = hpc(hpc_fus);

% find hpc ripple trials and corresponding fus trials
for subI = 1:size(hpc_fus_sub,2)
    x = 1;paddings=cell(1,10);
    k=0;invalid_rips(subI,2)=0;j=0;
    % find same ripple in different channels
    hpc_fus_sub(subI).unique_rips(:,13) = 1:size(hpc_fus_sub(subI).unique_rips,1);
    rips_sorted_by_chan_start{subI} = sortrows(hpc_fus_sub(subI).unique_rips,[2,3]);
    dif=diff(rips_sorted_by_chan_start{subI}(:,2:3));
    same_rip =find(abs(dif(:,2))<10 & dif(:,1)==0);
    rips_sorted_by_chan_start{subI}(same_rip,12) = 1;
    % return to original order to remove rip
    rips_sorted_by_chan_start{subI} = sortrows(rips_sorted_by_chan_start{subI},13);
    hpc_fus_sub(subI).same_rip=sum(rips_sorted_by_chan_start{subI}(:,12));
    % remove rip
    same_rip_rmv = find(hpc_fus_sub(subI).unique_rips(:,12));
    hpc_fus_sub(subI).unique_rips(same_rip_rmv,:) = [];
    
    % get hpc trls based on rip time or no rip
    idx_hpc_pre_stim_rip_trls =  hpc_fus_sub(subI).unique_rips(:,3) < 501;
    hpc_pre_stim_rip_trls = hpc_fus_sub(subI).unique_rips(idx_hpc_pre_stim_rip_trls ,11);
    hpc_pre_stim_rip_trls_unq = unique(hpc_pre_stim_rip_trls);
    hpc_post_stim_rip_trls_unq = unique(hpc_fus_sub(subI).unique_rips(~idx_hpc_pre_stim_rip_trls ,11));
    
    idx_hpc_no_rip_trls = ~ismember(hpc_fus_sub(subI).clean_trials(:,5),unique(hpc_fus_sub(subI).unique_rips(:,11)));
    hpc_no_rip_trls = hpc_fus_sub(subI).clean_trials(idx_hpc_no_rip_trls,5);
    
    % find corresponding fus trl for no rip trials
    fus_trls = fus(subI).clean_trials(:,5);
    no_rip_combined_trl = intersect(hpc_no_rip_trls, fus_trls);
    [~,fus_idx] = ismember(no_rip_combined_trl,fus_trls);
    for trlI = 1:numel(no_rip_combined_trl)
        fus_no_rip_demeaned{1,subI}.trial{1,trlI} = fus_demeaned{1, subI}.trial{:, fus_idx(trlI)}(:,501:end);
        fus_no_rip_baseline{1,subI}.trial{1,trlI} = fus_demeaned{1, subI}.trial{:, fus_idx(trlI)}(:,1:500);
        fus_no_rip_demeaned{1,subI}.time{1,trlI} = fus_demeaned{1, subI}.time{1,1}(501:end);
        fus_no_rip_baseline{1,subI}.time{1,trlI} = fus_demeaned{1, subI}.time{1,1}(1:500);
        no_infos(trlI,:) = fus_demeaned{1, subI}.preds(fus_idx(trlI),:);
    end
    fus_no_rip_demeaned{1,subI}.fsample = 500;
    fus_no_rip_demeaned{1,subI}.label = {'fus'};
    fus_no_rip_baseline{1,subI}.fsample = 500;
    fus_no_rip_baseline{1,subI}.label = {'fus'};
    
    % find corresponding fus trl for pre rip trials
    pre_rip_combined_trl = intersect(hpc_pre_stim_rip_trls_unq, fus_trls);
    [~,fus_idx_pre] = ismember(pre_rip_combined_trl,fus_trls);
    for trlI = 1:numel(pre_rip_combined_trl)
        fus_pre_rip_demeaned{1,subI}.trial{1,trlI} = fus_demeaned{1, subI}.trial{:, fus_idx_pre(trlI)}(:,501:end);
        fus_pre_rip_baseline{1,subI}.trial{1,trlI} = fus_demeaned{1, subI}.trial{:, fus_idx_pre(trlI)}(:,1:500);
        fus_pre_rip_demeaned{1,subI}.time{1,trlI} = fus_demeaned{1, subI}.time{1,1}(501:end);
        fus_pre_rip_baseline{1,subI}.time{1,trlI} = fus_demeaned{1, subI}.time{1,1}(1:500);
        pre_infos(trlI,:) = fus_demeaned{1, subI}.preds(fus_idx_pre(trlI),:);
    end
    fus_pre_rip_demeaned{1,subI}.fsample = 500;
    fus_pre_rip_demeaned{1,subI}.label ={'fus'};
    fus_pre_rip_baseline{1,subI}.fsample = 500;
    fus_pre_rip_baseline{1,subI}.label = {'fus'};
    
    l_foi = 35:2.5:160;%2.5:2.5:32.5;
    l_tw = 0.4*ones(length(l_foi),1)';
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.keeptrials   = 'yes';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'dpss';%'hanning'; %
    cfg.foi          = l_foi;
    cfg.pad          = 'maxperlen';
    cfg.tapsmofrq    = 10; % specifies half the spectral concentration
    cfg.t_ftimwin    = l_tw;
    cfg.toi          = 0:0.01:1.2;
    cfg.polyremoval  = 1;
    TF_pre_rip{subI} = ft_freqanalysis(cfg,fus_pre_rip_demeaned{1,subI});
    TF_no_rip{subI}  = ft_freqanalysis(cfg,fus_no_rip_demeaned{1,subI});
    
    cfg.toi          = -1.4:0.01:-1;
    base_pre_TF{subI}= ft_freqanalysis(cfg,fus_pre_rip_baseline{1,subI});
    base_no_TF{subI} = ft_freqanalysis(cfg,fus_no_rip_baseline{1,subI});

    
    % fit glm with predictors
    cfg              = [];
    cfg.avgoverchan  = 'yes';
    tf_avg_chan_pre{subI} = ft_selectdata(cfg, TF_pre_rip{subI});
    cfg              = [];
    cfg.latency  = [0 1];
    tf_avg_chan_pre{subI} = ft_selectdata(cfg,tf_avg_chan_pre{subI});
    xs = pre_infos(:,1:3);
    for freqI = 1:size(tf_avg_chan_pre{subI}.freq,2)
        for timeI = 1:size(tf_avg_chan_pre{subI}.time,2)
            b(1:size(xs,2)+1,freqI,timeI) = glmfit(xs,log(squeeze(tf_avg_chan_pre{subI}.powspctrm(:,1,freqI,timeI))));
        end
    end
    
    entropy = b(2,:,:);
    surprise = b(3,:,:);
    TF_ent_pre{subI}.powspctrm = entropy;
    TF_ent_pre{subI}.time = tf_avg_chan_pre{subI}.time;
    TF_ent_pre{subI}.freq = tf_avg_chan_pre{subI}.freq;
    TF_ent_pre{subI}.label = {'fus'};
    TF_ent_pre{subI}.dimord = 'chan_freq_time';
    TF_surp_pre{subI}.powspctrm = surprise;
    TF_surp_pre{subI}.time = tf_avg_chan_pre{subI}.time;
    TF_surp_pre{subI}.freq = tf_avg_chan_pre{subI}.freq;
    TF_surp_pre{subI}.label = {'fus'};
    TF_surp_pre{subI}.dimord = 'chan_freq_time';
    
    cfg              = [];
    cfg.avgoverchan  = 'yes';
    tf_avg_chan_no{subI} = ft_selectdata(cfg, TF_no_rip{subI});
    cfg              = [];
    cfg.latency  = [0 1];
    tf_avg_chan_no{subI} = ft_selectdata(cfg,tf_avg_chan_no{subI});
    xs_no = no_infos(:,1:3);
    for freqI = 1:size(tf_avg_chan_no{subI}.freq,2)
        for timeI = 1:size(tf_avg_chan_no{subI}.time,2)
            b_no(1:size(xs,2)+1,freqI,timeI) = glmfit(xs_no,log(squeeze(tf_avg_chan_no{subI}.powspctrm(:,1,freqI,timeI))));
        end
    end
    
    entropy_no = b_no(2,:,:);
    surprise_no = b_no(3,:,:);
    TF_ent_no{subI}.powspctrm = entropy_no;
    TF_ent_no{subI}.time = tf_avg_chan_no{subI}.time;
    TF_ent_no{subI}.freq = tf_avg_chan_no{subI}.freq;
    TF_ent_no{subI}.label = {'fus'};
    TF_ent_no{subI}.dimord = 'chan_freq_time';
    TF_surp_no{subI}.powspctrm = surprise_no;
    TF_surp_no{subI}.time = tf_avg_chan_no{subI}.time;
    TF_surp_no{subI}.freq = tf_avg_chan_no{subI}.freq;
    TF_surp_no{subI}.label = {'fus'};
    TF_surp_no{subI}.dimord = 'chan_freq_time';
    
        cfg = [];
    cfg.avgoverrpt = 'yes';
    TF_pre_rip{subI} = ft_selectdata(cfg, TF_pre_rip{subI});
    TF_no_rip{subI} = ft_selectdata(cfg, TF_no_rip{subI});
    base_pre_TF{subI} = ft_selectdata(cfg, base_pre_TF{subI});
    base_no_TF{subI} = ft_selectdata(cfg, base_no_TF{subI});
    
    cfg = [];
    cfg.avgoverchan = 'yes';
    TF_pre_rip_avg{subI} = ft_selectdata(cfg, TF_pre_rip{subI});
    TF_no_rip_avg{subI} = ft_selectdata(cfg, TF_no_rip{subI});
    base_pre_TF_avg{subI} = ft_selectdata(cfg, base_pre_TF{subI});
    base_no_TF_avg{subI} = ft_selectdata(cfg, base_no_TF{subI});
    
    cfg=[];cfg.masknans = 'yes';
    mean_pre = nanmean(base_pre_TF_avg{subI}.powspctrm,3); % 4 %trlxfreq baseline
    data_pre = TF_pre_rip_avg{subI}.powspctrm;
    TF_pre_rip_avg_base{subI}.powspctrm = ((data_pre-mean_pre)./mean_pre);
    
    mean_no = nanmean(base_no_TF_avg{subI}.powspctrm,3); % 4 %trlxfreq baseline
    data_no = TF_no_rip_avg{subI}.powspctrm;
    TF_no_rip_avg_base{subI}.powspctrm = ((data_no-mean_no)./mean_no);
    
    TF_pre_rip_avg_base{subI}.label = {'fus'};
    TF_pre_rip_avg_base{subI}.freq = TF_pre_rip{subI}.freq;
    TF_pre_rip_avg_base{subI}.time = TF_pre_rip{subI}.time;
    TF_pre_rip_avg_base{subI}.dimord = 'chan_freq_time';
    
    TF_no_rip_avg_base{subI}.label = {'fus'};
    TF_no_rip_avg_base{subI}.freq = TF_no_rip{subI}.freq;
    TF_no_rip_avg_base{subI}.time = TF_no_rip{subI}.time;
    TF_no_rip_avg_base{subI}.dimord = 'chan_freq_time';
    

    % timeseries of filtered gamma power
    avg_beta_surp_pre(subI,:) = squeeze(mean(TF_surp_pre{subI}.powspctrm,2))';
    avg_beta_surp_no(subI,:) = squeeze(mean(TF_surp_no{subI}.powspctrm,2))';
    
    avg_gamma_pre(subI,:) = squeeze(nanmean(TF_pre_rip_avg_base{subI}.powspctrm))';
    avg_gamma_no(subI,:) = squeeze(nanmean(TF_no_rip_avg_base{subI}.powspctrm));
    
    clear mean_no data_no mean_pre data_pre mean_after data_after fus_idx_pre fus_idx pre_infos no_infos xs xs_no b
end

%%

figure('position',[10 10 500 400]);bounds_pre = nanstd(avg_beta_surp_pre)/sqrt(numel(hpc_fus_sub));
bounds_no =nanstd(avg_beta_surp_no)/sqrt(numel(hpc_fus_sub));
boundedline(1:101, nanmean(avg_beta_surp_pre), bounds_pre, 'alpha', 'transparency', 0.1,'r'); hold on;
boundedline(1:101, nanmean(avg_beta_surp_no), bounds_no, 'alpha', 'transparency', 0.1); hold on;
ylabel('Surprise beta value');xlabel('Time (S)');xticks(1:10:101);
xticklabels([0:0.1:1]);legend({'','Pre-stim ripple','','No ripple'});
set(gca,'FontSize',20);

figure('position',[10 10 500 400]);
bounds_pre_gamma = nanstd(avg_gamma_pre)/sqrt(numel(hpc_fus_sub));
bounds_no_gamma =nanstd(avg_gamma_no)/sqrt(numel(hpc_fus_sub));
boundedline(1:101, nanmean(avg_gamma_pre(:,1:101)), bounds_pre_gamma, 'alpha', 'transparency', 0.1,'r'); hold on;
boundedline(1:101, nanmean(avg_gamma_no(:,1:101)), bounds_no_gamma, 'alpha', 'transparency', 0.1); hold on;
ylabel('Mean gamma (42.5-90Hz) power');xlabel('Time (S)');xticks(1:10:121);
xticklabels([0:0.1:1.2]);legend({'','Pre-stim ripple','','No ripple'});
set(gca,'FontSize',20);
%% stats
alpha=0.05;nrand='all';latency=[0 1];freqi=[35 160];%[42.5 95];%
cfg =[];
cfg.keepindividual = 'yes';
TF_pre_grandavg = ft_freqgrandaverage(cfg, TF_pre_rip_avg_base{:});
TF_no_grandavg = ft_freqgrandaverage(cfg, TF_no_rip_avg_base{:});

figure(1);cfg=[]; ft_singleplotTFR(cfg,TF_pre_grandavg);
title({'Grand average pre-stim ripples'}); set(gca, 'FontSize',14);
ylabel('Frequency (Hz)'); xlabel('Time relative to stim onset (S)');

figure(2);cfg=[]; ft_singleplotTFR(cfg,TF_no_grandavg);
title({'Grand average no ripples'}); set(gca, 'FontSize',14);
ylabel('Frequency (Hz)'); xlabel('Time relative to stim onset(S)');

figure(4)
pre_no =TF_GroupStats(TF_pre_grandavg,TF_no_grandavg,alpha,alpha,nrand,latency,freqi,'Same trial vs no',3, size(hpc_fus_sub,2));

%%
cfg = [];
for subI = 1:size(fus,2)
    figure(6);subplot(3,3,subI);ft_singleplotTFR(cfg,TF_pre_rip_avg_base{subI})
    figure(7);subplot(3,3,subI);ft_singleplotTFR(cfg,TF_no_rip_avg_base{subI})
end
figure(6);sgtitle('Pre-stim ripples');
figure(7);sgtitle('No ripples');

%%

for subI=1:size(hpc_fus,2)
    null{subI}.powspctrm  = zeros(1,size(TF_ent_pre{1}.powspctrm,2),size(TF_ent_pre{1}.powspctrm,3));
    null{subI}.time = TF_ent_pre{subI}.time;
    null{subI}.freq = TF_ent_pre{subI}.freq;
    null{subI}.dimord = 'chan_freq_time';
    null{subI}.label = {'fus'};
end
cfg=[];
cfg.keepindividual = 'yes';
cfg.toilim = [0 1];
null_grandavg =  ft_freqgrandaverage(cfg, null{:});
ent_pre_avg = ft_freqgrandaverage(cfg,TF_ent_pre{:});
surp_pre_avg = ft_freqgrandaverage(cfg,TF_surp_pre{:});
ent_stats=TF_GroupStats(ent_pre_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Pre-stim ripples Entropy',1, size(hpc_fus_sub,2));
surp_stats=TF_GroupStats(surp_pre_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Pre-stim ripples Surprise',1, size(hpc_fus_sub,2));

ent_no_avg = ft_freqgrandaverage(cfg,TF_ent_no{:});
surp_no_avg = ft_freqgrandaverage(cfg,TF_surp_no{:});
ent_no_stats=TF_GroupStats(ent_no_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'No ripples Entropy',1, size(hpc_fus_sub,2));
surp_no_stats=TF_GroupStats(surp_no_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'No ripples Surprise',1, size(hpc_fus_sub,2));


pre_minus_no_surp = surp_pre_avg;
pre_minus_no_surp.powspctrm = surp_pre_avg.powspctrm - surp_no_avg.powspctrm;
figure();ft_singleplotTFR(cfg,pre_minus_no_surp);title('pre vs no beta surprise');set(gca,'FontSize',14);
surp_pre_no_stats=TF_GroupStats(pre_minus_no_surp,null_grandavg,alpha,alpha,nrand,latency,freqi,'Pre vs no ripples Surprise',1, size(hpc_fus_sub,2));

subplot(121);ft_singleplotTFR(cfg,surp_pre_avg);title('pre');subplot(122);ft_singleplotTFR(cfg,surp_no_avg);title('no ripples')

%%
% double dissociation of the clusters
rip_surp_mask_rip = surp_pre_avg.powspctrm(:,surp_stats.mask);
rip_surp_mask_no_rip = surp_pre_avg.powspctrm(:,surp_no_stats.mask);
no_rip_surp_mask_no_rip = surp_no_avg.powspctrm(:,surp_no_stats.mask);
no_rip_surp_mask_rip = surp_no_avg.powspctrm(:,surp_stats.mask);

avgs = [mean(rip_surp_mask_rip,2),mean(rip_surp_mask_no_rip,2),...
    mean(no_rip_surp_mask_no_rip,2),mean(no_rip_surp_mask_rip,2)];

% comparison of power in sig cluster
load('fus_stats.mat')
[r,c] = find(squeeze(stats.Surprise_high_post.mask));
min_time = min(r); max_time = max(r);
min_freq = min(c); max_freq = max(c);
rip_surp_mask_clust = surp_pre_avg.powspctrm(:,squeeze(stats.Surprise_high_post.mask));
no_rip_surp_mask_clust = surp_no_avg.powspctrm(:,squeeze(stats.Surprise_high_post.mask));

clust_avgs = [mean(rip_surp_mask_clust,2),mean(no_rip_surp_mask_clust,2)];


