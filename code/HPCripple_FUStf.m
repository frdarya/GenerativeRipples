clearvars

global ft_default
ft_default.spmversion = 'spm12';
ft_defaults
ft_hastoolbox('brewermap', 1);
fus_gamma_range = [40 80];

region = 'anterior'; baseline = 0; resp_lock = 0;
nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];
iszurich = logical([zeros(1,13),ones(1,3)]);
hpfilt = 200;
ripdur = 25;
fname = sprintf('HPCRipples/HPCAnterior_vaz_hpf%d_%dms_%dsubjs.mat',hpfilt, ripdur, numel(nsubs));
load(fname)

% create stp and patient_data structs
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);

for subI = 1:numel(nsubs)
    if nsubs(subI) == 22 % no ripples for sub 22 (and excluded anyways) - position 8
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
    sub(subI).clean_trials = clean.info;
    cd ../../../
    stp.poststim = 1;
    stp.prestim = 1;
    timest={-stp.prestim:0.002:stp.poststim};
    
    cleandata = preproc(subI,nsubs,stp,patient_data,clean.trl);
    if stp.zurich(subI) == 1 || nsubs(subI) == 31 || nsubs(subI) == 32 || ...
            nsubs(subI) == 36 || nsubs(subI) == 37
        cfg         = [];
        cfg.time    = repmat({timest{1}(1:end-1)},1,size(cleandata.trial,2));
        cleandata   = ft_resampledata(cfg,cleandata);
    end
    
    cfg             = [];
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 20;
    cfg.demean      = 'yes';
    hpc_demeaned_low{subI} = ft_preprocessing(cfg,cleandata);
    clear cleandata
end
hpc = sub;

%% load fus
clear sub

% get ts
region = 'fus'; baseline = 0; resp_lock = 0;
nsubs = [3,8,13,15,31,32,36];%37 too noisy
iszurich = logical([zeros(1,7)]);
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);
for subI = 1:numel(nsubs)
    if nsubs(subI) == 22
        continue
    end
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
    stp.poststim = 1;
    stp.prestim = 1.8;
    timest={-stp.prestim:0.002:stp.poststim};
    
    cleandata = preproc(subI,nsubs,stp,patient_data,clean.trl);
    cfg         = [];
    cfg.time    = repmat({timest{1}(1:end-1)},1,size(cleandata.trial,2));
    cleandata   = ft_resampledata(cfg,cleandata);
    
    cfg             = [];
    cfg.demean      = 'yes';
    fus_demeaned{subI} = ft_preprocessing(cfg,cleandata);
    
    clear cleandata
    
end
fus = sub;

%% subset hpc patients for those with fus contacts and get corresponding trls
hpc_fus=[1,3,5,6,10,11,12]; %37 too noisy
% indices to keep
hpc_fus_sub = hpc(hpc_fus);

hpc_demeaned_low_fus = hpc_demeaned_low(hpc_fus);
infos = {};fus_4tf = {}; fus_base_4tf = {}; TF_corrected = {};
hpc_rip_peak = {}; hpc_rip_trough = {}; hpc_timeseries_low = {};
% find hpc ripple trials and corresponding fus trials
for subI = 1:numel(nsubs)
    x = 1;paddings=cell(1,10);
    
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
    
    
    hpc_idx_trials = hpc_fus_sub(subI).unique_rips(:,2);%trial index based on clean_trials
    hpc_trials = hpc_fus_sub(subI).unique_rips(:,11); % trial number, not index!
    fus_trls = fus(subI).clean_trials(:,5);
    combined_trl = intersect(hpc_trials, fus_trls);
    comb(subI) = size(combined_trl,1);
    fprintf(['\n', 'combined trl ', num2str(size(combined_trl,1)), '\n'])
    [~,fus_idx] = ismember(combined_trl,fus_trls);
    
    % create two cells with the data for subsequent analysis
    % the fus_gamma has the time series of gamma-band activity
    % the hpc_rip has the ripple info ripple start & peak will be used as
    % indices to cut fus gamma time-series accordingly
    rip_trl= 1;po=1;pr=1;peaki = 1; troughi=1;pr8=1; peaki_pre=1;peaki_post=1;
    troughi_pre=1; troughi_post=1;prtf=1;potf=1;pr8tf=1;
    for trlI = 1:numel(combined_trl)
        fus_demeaned_trl{1,trlI} = fus_demeaned{1, subI}.trial{1, fus_idx(trlI)}(:,401:end);
        fus_baseline_trl{1,trlI} = fus_demeaned{1, subI}.trial{1, fus_idx(trlI)}(:,100:500); % for baseline -1.4 to -1, estimate -1,6 to -0.8
        hpc_low_trl{1,trlI} = hpc_demeaned_low_fus{1, subI}.trial{1, hpc_idx_trials(trlI)};
        hpc_idx(1,trlI) = {find(combined_trl(trlI)==hpc_trials)};
        hpc_idx{2,trlI} = combined_trl(trlI);
        hpc_rip_trl{1,trlI} = hpc_fus_sub(subI).unique_rips(hpc_idx{1,trlI},:);
        
        % sanity check
        if hpc_rip_trl{1,trlI}(1,11) ~= fus_trls(fus_idx(trlI))
            break
        end
        for ripI=1:size(hpc_rip_trl{1,trlI},1)
            % use hpc ripple time to cut fus gamma ts
            rip_info = hpc_rip_trl{1,trlI}(ripI,3:10);
            fus_start = rip_info(1,1)-100;
            fus_end = rip_info(1,2)+100;
            rip_peak = rip_info(1,3); %1=start, 2=end, 3 = peak
            rip_dur =rip_info(1,2)-rip_info(1,1)+2;
            rip_elec = hpc_rip_trl{1,trlI}(ripI,1);
            
            full_trl =linspace(-400,400,401)/1000;
            % get unified +- 400ms time vector for tf
            if rip_peak < 201 % ripple in the beginnig
                nan_pad = 201-rip_peak;
                time_idx = [1:rip_peak+200];
                fus_ts_4tf= [zeros(size(fus_demeaned_trl{1,trlI},1),nan_pad),fus_demeaned_trl{1,trlI}(:,time_idx)];
                paddings{rip_trl} = 1:nan_pad;
            elseif rip_peak > 800 % ripple in the end
                time_idx = [rip_peak-200:1000];
                nan_pad = 401-(length(time_idx)); %there are 201 timepoints - 400ms
                fus_ts_4tf = [fus_demeaned_trl{1,trlI}(:,time_idx),zeros(size(fus_demeaned_trl{1,trlI},1),nan_pad)];
                paddings{rip_trl} = 1000-nan_pad:1000;
            else
                time_idx = [rip_peak-200:rip_peak+200];
                fus_ts_4tf = fus_demeaned_trl{1,trlI}(:,time_idx);
                
            end
            
            fus_4tf{1,subI}.trial{rip_trl} =  fus_ts_4tf;
            fus_4tf{1,subI}.time{rip_trl} = full_trl;
            fus_4tf{1,subI}.pred(rip_trl,:) = rip_info(5:7);
            
            fus_base_4tf{1,subI}.trial{rip_trl} = fus_baseline_trl{1,trlI};
            fus_base_4tf{1,subI}.time{rip_trl} = linspace(-400,400,401)/1000;
            fus_base_4tf{1,subI}.pred(rip_trl,:) = rip_info(5:7);
            
            
            if rip_peak < 500
                fus_4tf_pre{1,subI}.trial{prtf} = fus_4tf{1,subI}.trial{rip_trl};
                fus_4tf_pre{1,subI}.time{prtf} = fus_4tf{1,subI}.time{rip_trl};
                fus_4tf_pre{1,subI}.pred{prtf} = fus_4tf{1,subI}.pred(rip_trl,:);
                fus_base_4tf_pre{1,subI}.trial{prtf} = fus_baseline_trl{1,trlI};
                fus_base_4tf_pre{1,subI}.time{prtf} = linspace(-400,400,401)/1000;
                fus_base_4tf_pre{1,subI}.pred(prtf,:) = rip_info(5:7);
                prtf=prtf+1;
            else
                fus_4tf_post{1,subI}.trial{potf} = fus_4tf{1,subI}.trial{rip_trl};
                fus_4tf_post{1,subI}.time{potf} = fus_4tf{1,subI}.time{rip_trl};
                fus_4tf_post{1,subI}.pred{potf} = fus_4tf{1,subI}.pred(rip_trl,:);
                fus_base_4tf_post{1,subI}.trial{potf} = fus_baseline_trl{1,trlI};
                fus_base_4tf_post{1,subI}.time{potf} = linspace(-400,400,401)/1000;
                fus_base_4tf_post{1,subI}.pred(potf,:) = rip_info(5:7);
                potf=potf+1;
            end
            
            infos{1,subI}.pred(rip_trl,:) = rip_info(5:7);
            
            x= x+1; rip_trl = rip_trl+1;
            clear rip_info fus_sart fus_end fus_ts fus_ts_4tf fus_ts_gamma fus_ts_gamma_short fus_ts_gamma_baseline hpc_rip_short
        end
    end
    
    %check no nans
    base_nan =cell2mat(cellfun(@find,(cellfun(@isnan,fus_4tf{1,subI}.trial,'UniformOutput',false)),'UniformOutput',false));
    base_inf =cell2mat(cellfun(@find,(cellfun(@isinf,fus_4tf{1,subI}.trial,'UniformOutput',false)),'UniformOutput',false));
    
    if ~isempty(base_nan) || ~isempty(base_inf)
        warning('nan/inf in data!')
    end
    
    fus_4tf{1,subI}.fsample = 500;
    fus_4tf{1,subI}.label = fus_demeaned{1,subI}.label;
    fus_4tf_pre{1,subI}.fsample = 500;
    fus_4tf_pre{1,subI}.label = fus_demeaned{1,subI}.label;
    fus_4tf_post{1,subI}.fsample = 500;
    fus_4tf_post{1,subI}.label = fus_demeaned{1,subI}.label;
    fus_base_4tf{1,subI}.fsample = 500;
    fus_base_4tf{1,subI}.label = fus_demeaned{1,subI}.label;
    fus_base_4tf_pre{1,subI}.fsample = 500;
    fus_base_4tf_pre{1,subI}.label = fus_demeaned{1,subI}.label;
    fus_base_4tf_post{1,subI}.fsample = 500;
    fus_base_4tf_post{1,subI}.label = fus_demeaned{1,subI}.label;
    
    
    l_foi =35:2.5:160;%2.5:2.5:32.5%;
    l_tw = 0.4*ones(length(l_foi),1)';
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.keeptrials   = 'yes';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'dpss';%
    cfg.foi          = l_foi;
    cfg.pad          = 'maxperlen';
    cfg.tapsmofrq    = 10; % specifies half the spectral concentration
    
    cfg.t_ftimwin    = l_tw;
    cfg.toi          = -0.2:0.01:0.2;
    cfg.polyremoval  = 1;
    TF{subI}         = ft_freqanalysis(cfg,fus_4tf{1,subI});
    preTF{subI}      = ft_freqanalysis(cfg,fus_4tf_pre{1,subI});
    postTF{subI}     = ft_freqanalysis(cfg,fus_4tf_post{1,subI});
    baseTF{subI}     = ft_freqanalysis(cfg,fus_base_4tf{1,subI});
    prebaseTF{subI}  = ft_freqanalysis(cfg,fus_base_4tf_pre{1,subI});
    postbaseTF{subI} = ft_freqanalysis(cfg,fus_base_4tf_post{1,subI});
        
    % fit glm with predictors
    cfg              = [];
    cfg.avgoverchan  = 'yes';
    pre_avg_chan{subI} = ft_selectdata(cfg, preTF{subI});
    post_avg_chan{subI} = ft_selectdata(cfg, postTF{subI});
    tf_avg_chan{subI} = ft_selectdata(cfg, TF{subI});
    
    xs = infos{1,subI}.pred(:,1:3);
    for i = 1:size(fus_4tf_pre{subI}.pred,2)
        pre_xs(i,:) = fus_4tf_pre{subI}.pred{i};
    end
    for i = 1:size(fus_4tf_post{subI}.pred,2)
        post_xs(i,:) = fus_4tf_post{subI}.pred{i};
    end
    
    for freqI = 1:size(tf_avg_chan{subI}.powspctrm,3)
        for timeI = 1:size(tf_avg_chan{subI}.powspctrm,4)
            b(1:size(xs,2)+1,freqI,timeI) = glmfit(xs,log(squeeze(tf_avg_chan{subI}.powspctrm(:,1,freqI,timeI))));
            b_pre(1:size(pre_xs,2)+1,freqI,timeI) = glmfit(pre_xs,log(squeeze(pre_avg_chan{subI}.powspctrm(:,1,freqI,timeI))));
            b_post(1:size(post_xs,2)+1,freqI,timeI) = glmfit(post_xs,log(squeeze(post_avg_chan{subI}.powspctrm(:,1,freqI,timeI))));
            
        end
    end
    entropy = b(2,:,:);
    surprise = b(3,:,:);
    pre_entropy = b_pre(2,:,:);
    pre_surprise = b_pre(3,:,:);
    post_entropy = b_post(2,:,:);
    post_surprise = b_post(3,:,:);
    
    TF_ent{subI}.powspctrm = entropy;
    TF_ent{subI}.time = tf_avg_chan{subI}.time;
    TF_ent{subI}.freq = tf_avg_chan{subI}.freq;
    TF_ent{subI}.label = {stp.region};
    TF_ent{subI}.dimord = 'chan_freq_time';
    preTF_ent{subI}.powspctrm = pre_entropy;
    preTF_ent{subI}.time = tf_avg_chan{subI}.time;
    preTF_ent{subI}.freq = tf_avg_chan{subI}.freq;
    preTF_ent{subI}.label = {stp.region};
    preTF_ent{subI}.dimord = 'chan_freq_time';
    postTF_ent{subI}.powspctrm = post_entropy;
    postTF_ent{subI}.time = tf_avg_chan{subI}.time;
    postTF_ent{subI}.freq = tf_avg_chan{subI}.freq;
    postTF_ent{subI}.label = {stp.region};
    postTF_ent{subI}.dimord = 'chan_freq_time';
    
    TF_surp{subI}.powspctrm = surprise;
    TF_surp{subI}.time = tf_avg_chan{subI}.time;
    TF_surp{subI}.freq = tf_avg_chan{subI}.freq;
    TF_surp{subI}.label = {stp.region};
    TF_surp{subI}.dimord = 'chan_freq_time';
    preTF_surp{subI}.powspctrm = pre_surprise;
    preTF_surp{subI}.time = tf_avg_chan{subI}.time;
    preTF_surp{subI}.freq = tf_avg_chan{subI}.freq;
    preTF_surp{subI}.label = {stp.region};
    preTF_surp{subI}.dimord = 'chan_freq_time';
    postTF_surp{subI}.powspctrm = post_surprise;
    postTF_surp{subI}.time = tf_avg_chan{subI}.time;
    postTF_surp{subI}.freq = tf_avg_chan{subI}.freq;
    postTF_surp{subI}.label = {stp.region};
    postTF_surp{subI}.dimord = 'chan_freq_time';
    
    clear b xs pre_xs post_xs
    
    cfg = [];
    cfg.avgoverrpt = 'yes';
    TF_corr{subI} = ft_selectdata(cfg, TF{subI});
    preTF_corr{subI} = ft_selectdata(cfg, preTF{subI});
    postTF_corr{subI} = ft_selectdata(cfg, postTF{subI});
    baseTF_corr{subI}     = ft_selectdata(cfg,baseTF{subI});
    prebaseTF_corr{subI}  = ft_selectdata(cfg,prebaseTF{subI});
    postbaseTF_corr{subI} = ft_selectdata(cfg,postbaseTF{subI});

    cfg              = [];
    cfg.avgoverchan  = 'yes';
    TF_avg{subI}     = ft_selectdata(cfg,TF_corr{subI});
    preTF_avg{subI}  = ft_selectdata(cfg,preTF_corr{subI});
    postTF_avg{subI} = ft_selectdata(cfg,postTF_corr{subI});
    baseTF_avg{subI}     = ft_selectdata(cfg,baseTF_corr{subI});
    prebaseTF_avg{subI}  = ft_selectdata(cfg,prebaseTF_corr{subI});
    postbaseTF_avg{subI} = ft_selectdata(cfg,postbaseTF_corr{subI});

    cfg=[];cfg.masknans = 'yes';
    mean = nanmean(baseTF_avg{subI}.powspctrm,3); % 4 %trlxfreq baseline
    data = TF_avg{subI}.powspctrm;
    TF_avg{subI}.powspctrm = ((data-mean)./mean);
    mean_pre = nanmean(prebaseTF_avg{subI}.powspctrm,3); %4 %trlxfreq baseline
    data_pre = preTF_avg{subI}.powspctrm;
    preTF_avg{subI}.powspctrm = ((data_pre-mean_pre)./mean_pre);%10*log10(data_pre ./ mean_pre);
    mean_post = nanmean(postbaseTF_avg{subI}.powspctrm,3); %4 %trlxfreq baseline
    data_post = postTF_avg{subI}.powspctrm;
    postTF_avg{subI}.powspctrm = ((data_post-mean_post)./mean_post);%10*log10(data_post ./ mean_post);

    TF_avg{subI}.label     = {'fus'};
    preTF_avg{subI}.label     = {'fus'};
    postTF_avg{subI}.label     = {'fus'};

    clear fus_gamma_trl hpc_low_trl hpc_idx hpc_rip_trl...
        combined_trl peak_trl trough_trl mean mean_pre mean_post...
        data data_pre data_post
    
end


%% TF agroupavg
for subI=1:numel(nsubs)
    null{subI}.powspctrm  = zeros(1,size(TF_avg{1}.freq,2),41);
    null{subI}.time = TF_avg{subI}.time;
    null{subI}.freq = TF_avg{subI}.freq;
    null{subI}.dimord = 'chan_freq_time';
    null{subI}.label = {'fus'};
end
alpha=0.05;nrand='all';latency=[-0.2 0.2];freqi=[35 160];
%% betas
cfg=[];
cfg.keepindividual = 'yes';
null_grandavg =  ft_freqgrandaverage(cfg, null{:});
ent_avg = ft_freqgrandaverage(cfg,TF_ent{:});
surp_avg = ft_freqgrandaverage(cfg,TF_surp{:});
pre_ent_avg = ft_freqgrandaverage(cfg,preTF_ent{:});
pre_surp_avg = ft_freqgrandaverage(cfg,preTF_surp{:});
post_ent_avg = ft_freqgrandaverage(cfg,postTF_ent{:});
post_surp_avg = ft_freqgrandaverage(cfg,postTF_surp{:});
ent_stats=TF_GroupStats(ent_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Entropy',1, size(fus,2));
surp_stats=TF_GroupStats(surp_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Surprise',1, size(fus,2));
pre_ent_stats=TF_GroupStats(pre_ent_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Entropy pre-stim',1, size(fus,2));
pre_surp_stats=TF_GroupStats(pre_surp_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Surprise pre-stim',1, size(fus,2));
post_ent_stats=TF_GroupStats(post_ent_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Entropy post-stim',1, size(fus,2));
post_surp_stats=TF_GroupStats(post_surp_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Surprisepost-stim ',1, size(fus,2));
pre_post_ent_stats=TF_GroupStats(pre_ent_avg,post_ent_avg,alpha,alpha,nrand,latency,freqi,'Entropy pre vs post',1, size(fus,2));
pre_post_surp_stats=TF_GroupStats(pre_surp_avg,post_surp_avg,alpha,alpha,nrand,latency,freqi,'Surprise pre vs post',1, size(fus,2));

%%
cfg =[];
cfg.keepindividual = 'yes';
% cfg.foilim = [50, 150]
TF_grandavg = ft_freqgrandaverage(cfg, TF_avg{:});
null_grandavg =  ft_freqgrandaverage(cfg, null{:});
figure(3);cfg=[]; ft_singleplotTFR(cfg,TF_grandavg);
title({'Grand average','using relative baseline -1.4 to -1'}); set(gca, 'FontSize',14);
ylabel('Frequency (Hz)'); xlabel('Time relative to HPC ripple (S)');

overall_stats=TF_GroupStats(TF_grandavg,null_grandavg,alpha,alpha,nrand,latency,freqi,'rip_fus_overall',1,size(fus,2));

cfg =[];
cfg.keepindividual = 'yes';
% cfg.foilim = [50, 150]
TF_grandavg_pre = ft_freqgrandaverage(cfg, preTF_avg{:});
figure(4);cfg=[]; ft_singleplotTFR(cfg,TF_grandavg_pre)
title({'Pre stim','using relative baseline -1.4 to -1'}); set(gca, 'FontSize',14);
ylabel('Frequency (Hz)'); xlabel('Time relative to HPC ripple (S)');
pre_stats=TF_GroupStats(TF_grandavg_pre,null_grandavg,alpha,alpha,nrand,latency,freqi,'rip_fus_pre',2, size(fus,2));


cfg =[];
cfg.keepindividual = 'yes';
% cfg.foilim = [50, 150]
TF_grandavg_post = ft_freqgrandaverage(cfg, postTF_avg{:});
figure(5);cfg=[]; ft_singleplotTFR(cfg,TF_grandavg_post)
title({'Post stim','using relative baseline -1.4 to -1'}); set(gca, 'FontSize',14);
ylabel('Frequency (Hz)'); xlabel('Time relative to HPC ripple (S)');
post_stats=TF_GroupStats(TF_grandavg_post,null_grandavg,alpha,alpha,nrand,latency,freqi,'rip_fus_post',3, size(fus,2));

pre_post =TF_GroupStats(TF_grandavg_pre,TF_grandavg_post,alpha,alpha,nrand,latency,freqi,'rip_fus_pre_post',3, size(fus,2));

pre_post_raw = TF_grandavg_post;
pre_post_raw.powspctrm = TF_grandavg_pre.powspctrm - TF_grandavg_post.powspctrm


%%
clear tvals pvals maxt_beta_sub
overall_mask = overall_stats.posclusterslabelmat==1;
tvals = overall_stats.stat(overall_mask);
maxt = max(abs(tvals));
idx_maxt = abs(tvals)==maxt;
pvals = overall_stats.prob(overall_mask);
maxp = pvals(idx_maxt);

for subI = 1:numel(hpc_fus_sub)
    ov = [squeeze(TF_grandavg.powspctrm(subI,:,:,:))];
    ov_mask = ov(overall_mask);
    maxt_beta_sub(subI) = ov_mask(idx_maxt);
end
cohens_d = mean(maxt_beta_sub)/std(maxt_beta_sub);

largest_t_p_d.overall = [maxt,maxp,cohens_d];
clear tvals pvals maxt_beta_sub

pre_mask = pre_stats.posclusterslabelmat==1;
tvals = pre_stats.stat(pre_mask);
maxt = max(abs(tvals));
idx_maxt = abs(tvals)==maxt;
pvals = pre_stats.prob(pre_mask);
maxp = pvals(idx_maxt);

for subI = 1:numel(hpc_fus_sub)
    ovp = [squeeze(TF_grandavg_pre.powspctrm(subI,:,:,:))];
    ovp_mask = ovp(pre_mask);
    maxt_beta_sub(subI) = ovp_mask(idx_maxt);
end
cohens_d = mean(maxt_beta_sub)/std(maxt_beta_sub);

largest_t_p_d.pre = [maxt,maxp,cohens_d];
