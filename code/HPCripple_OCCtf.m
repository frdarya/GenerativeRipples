clearvars

global ft_default
ft_default.spmversion = 'spm12';
ft_defaults
ft_hastoolbox('brewermap', 1);
occ_gamma_range = [60, 160];

load('HPCRipples/HPCAnterior_vaz_hpf200_25ms_16subjs_Jan22_reclean.mat')
region = 'anterior'; baseline = 0; resp_lock = 0;
nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];
hpc_nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];

millis = linspace(-1,1,1001)*1000;
iszurich = logical([zeros(1,13),ones(1,3)]);

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
        % downsample - for zurich & 31,36 ruber
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

%% load occ
clear sub
% get ripples - although prob not used?
load('OCCRipples/vaz_hpf200_25ms.mat')

% get ts
region = 'occ'; baseline = 0; resp_lock = 0;
nsubs = [3,8,9,16,25,31,36,37];
iszurich = logical([zeros(1,8)]);
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);
for subI = 1:numel(nsubs)
    if nsubs(subI) == 22 % no ripples for sub 22 (and excluded anyways) - position 8
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
    %timest = cleandata.time(1,1);
    %     if stp.zurich(subI) == 1 || nsubs(subI) == 31 || nsubs(subI) == 36
    % downsample - for zurich & 31,36 ruber
    cfg         = [];
    cfg.time    = repmat({timest{1}(1:end-1)},1,size(cleandata.trial,2));
    cleandata   = ft_resampledata(cfg,cleandata);
    %     end
    
    cfg             = [];
    cfg.demean      = 'yes';
    occ_demeaned{subI} = ft_preprocessing(cfg,cleandata);
    cfg             = [];
    cfg.avgoverchan      = 'yes';
    occ_demeaned{subI} = ft_selectdata(cfg,occ_demeaned{subI});
    cfg             = [];
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = occ_gamma_range;
    cfg.demean      = 'yes';
    cfg.hilbert     = 'abs';
    occ_demeaned_gamma{subI} = ft_preprocessing(cfg,cleandata);
    clear cleandata
end
occ = sub;

%% subset hpc patients for those with occ contacts and get corresponding trls
hpc_occ = [1,3,4,7,9,10,12,13]; % indices to keep
if sum(~ismember(hpc_nsubs(hpc_occ),nsubs)) ~= 0
    error('wrong patients subsetting!')
end
hpc_occ_sub = hpc(hpc_occ);

hpc_demeaned_low_occ = hpc_demeaned_low(hpc_occ);

infos = {};occ_4tf = {}; occ_base_4tf = {}; TF_corrected = {};
hpc_rip_peak = {}; hpc_rip_trough = {}; hpc_timeseries_low = {};
% find hpc ripple trials and corresponding occ trials
for subI = 1:size(hpc_occ_sub,2)
    x = 1;paddings=cell(1,10);
    k=0;invalid_rips(subI,2)=0;j=0;
    % find same ripple in different channels
    hpc_occ_sub(subI).unique_rips(:,13) = 1:size(hpc_occ_sub(subI).unique_rips,1);
    rips_sorted_by_chan_start{subI} = sortrows(hpc_occ_sub(subI).unique_rips,[2,3]);
    dif=diff(rips_sorted_by_chan_start{subI}(:,2:3));
    same_rip =find(abs(dif(:,2))<10 & dif(:,1)==0);
    rips_sorted_by_chan_start{subI}(same_rip,12) = 1;
    % return to original order to remove rip
    rips_sorted_by_chan_start{subI} = sortrows(rips_sorted_by_chan_start{subI},13);
    hpc_occ_sub(subI).same_rip=sum(rips_sorted_by_chan_start{subI}(:,12));
    % remove rip
    same_rip_rmv = find(hpc_occ_sub(subI).unique_rips(:,12));
    
    hpc_occ_sub(subI).unique_rips(same_rip_rmv,:) = [];
    
    
    hpc_idx_trials = hpc_occ_sub(subI).unique_rips(:,2);%trial index based on clean_trials
    if nsubs(subI) ~= 37
        hpc_trials = hpc_occ_sub(subI).unique_rips(:,11); % trial number, not index!
    else
        hpc_trials = hpc_occ_sub(subI).unique_rips(:,11)+1; % trial number, not index!
        
    end
    occ_trls = occ(subI).clean_trials(:,5);
    combined_trl = intersect(hpc_trials, occ_trls);
    comb(subI) = size(combined_trl,1);
    fprintf(['\n', 'combined trl ', num2str(size(combined_trl,1)), '\n'])
    [~,occ_idx] = ismember(combined_trl,occ_trls);
    
    % create two cells with the data for subsequent analysis
    % the occ_gamma has the time series of gamma-band activity
    % the hpc_rip has the ripple info ripple start & peak will be used as
    % indices to cut occ gamma time-series accordingly
    rip_trl= 1;po=1;pr=1;peaki = 1; troughi=1;pr8=1; peaki_pre=1;peaki_post=1;
    troughi_pre=1; troughi_post=1;prtf=1;potf=1;pr8tf=1;
    for trlI = 1:numel(combined_trl)
        occ_demeaned_trl{1,trlI} = occ_demeaned{1, subI}.trial{1, occ_idx(trlI)}(:,401:end);
        occ_baseline_trl{1,trlI} = occ_demeaned{1, subI}.trial{1, occ_idx(trlI)}(:,100:500); % for baseline -1.4 to -1, estimate -1,6 to -0.8
        occ_gamma_trl{1,trlI} = occ_demeaned_gamma{1, subI}.trial{1, occ_idx(trlI)};
        hpc_low_trl{1,trlI} = hpc_demeaned_low_occ{1, subI}.trial{1, hpc_idx_trials(trlI)};
        hpc_idx(1,trlI) = {find(combined_trl(trlI)==hpc_trials)};
        hpc_idx{2,trlI} = combined_trl(trlI);
        hpc_rip_trl{1,trlI} = hpc_occ_sub(subI).unique_rips(hpc_idx{1,trlI},:);
        invalid_rips(subI,2) =  invalid_rips(subI,2)+size(hpc_rip_trl{1,trlI},1);
        
        % sanity check
        if nsubs(subI) ~= 37
            if hpc_rip_trl{1,trlI}(1,11) ~= occ_trls(occ_idx(trlI))
                break
            end
        else
            if hpc_rip_trl{1,trlI}(1,11)+1 ~= occ_trls(occ_idx(trlI))
                break
            end
        end
        
        for ripI=1:size(hpc_rip_trl{1,trlI},1)
            
            %
            % use hpc ripple time to cut occ gamma ts
            rip_info = hpc_rip_trl{1,trlI}(ripI,3:10);
            rt = rip_info(1,end);
            
            if millis(rip_info(1,2)) > rt  || rip_info(1,2) < 500
                invalid_rips(subI,1) = k+1;
                k=k+1;
            else
                j=j+1;
            end
            occ_start = rip_info(1,1)-100;
            occ_end = rip_info(1,2)+100;
            rip_peak = rip_info(1,3); %1=start, 2=end, 3 = peak
            rip_dur =rip_info(1,2)-rip_info(1,1)+2;
            rip_elec = hpc_rip_trl{1,trlI}(ripI,1)
            
            full_trl =linspace(-400,400,401)/1000;
            % get unified +- 400ms time vector for tf
            if rip_peak < 201 % ripple in the beginnig
                zero_pad = 201-rip_peak;
                time_idx = [1:rip_peak+200];
                occ_ts_4tf= [zeros(size(occ_demeaned_trl{1,trlI},1),zero_pad),occ_demeaned_trl{1,trlI}(:,time_idx)];
                occ_ts_gamma = [zeros(size(occ_gamma_trl{1,trlI},1),zero_pad),occ_gamma_trl{1,trlI}(:,time_idx)];
                paddings{rip_trl} = 1:zero_pad;
            elseif rip_peak > 800 % ripple in the end
                time_idx = [rip_peak-200:1000];
                zero_pad = 401-(length(time_idx)); %there are 201 timepoints - 400ms
                occ_ts_4tf = [occ_demeaned_trl{1,trlI}(:,time_idx),zeros(size(occ_demeaned_trl{1,trlI},1),zero_pad)];
                occ_ts_gamma = [occ_gamma_trl{1,trlI}(:,time_idx),zeros(size(occ_gamma_trl{1,trlI},1),zero_pad)];
                
                paddings{rip_trl} = 1000-zero_pad:1000;
            else
                time_idx = [rip_peak-200:rip_peak+200];
                occ_ts_4tf = occ_demeaned_trl{1,trlI}(:,time_idx);
                occ_ts_gamma = occ_gamma_trl{1,trlI}(:,time_idx);
                
            end
            
            occ_4tf{1,subI}.trial{rip_trl} =  occ_ts_4tf;
            occ_4tf{1,subI}.time{rip_trl} = full_trl;
            occ_4tf{1,subI}.pred(rip_trl,:) = rip_info(5:7);
            
            occ_base_4tf{1,subI}.trial{rip_trl} = occ_baseline_trl{1,trlI};
            occ_base_4tf{1,subI}.time{rip_trl} = linspace(-400,400,401)/1000;
            occ_base_4tf{1,subI}.pred(rip_trl,:) = rip_info(5:7);
            
            if rip_peak < 500
                occ_4tf_pre{1,subI}.trial{prtf} = occ_4tf{1,subI}.trial{rip_trl};
                occ_4tf_pre{1,subI}.time{prtf} = occ_4tf{1,subI}.time{rip_trl};
                occ_4tf_pre{1,subI}.pred{prtf} = occ_4tf{1,subI}.pred(rip_trl,:);
                occ_base_4tf_pre{1,subI}.trial{prtf} = occ_baseline_trl{1,trlI};
                occ_base_4tf_pre{1,subI}.time{prtf} = linspace(-400,400,401)/1000;
                occ_base_4tf_pre{1,subI}.pred(prtf,:) = rip_info(5:7);
                prtf=prtf+1;
            else
                occ_4tf_post{1,subI}.trial{potf} = occ_4tf{1,subI}.trial{rip_trl};
                occ_4tf_post{1,subI}.time{potf} = occ_4tf{1,subI}.time{rip_trl};
                occ_4tf_post{1,subI}.pred{potf} = occ_4tf{1,subI}.pred(rip_trl,:);
                occ_base_4tf_post{1,subI}.trial{potf} = occ_baseline_trl{1,trlI};
                occ_base_4tf_post{1,subI}.time{potf} = linspace(-400,400,401)/1000;
                occ_base_4tf_post{1,subI}.pred(potf,:) = rip_info(5:7);
                potf=potf+1;
            end
            
            infos{1,subI}.pred(rip_trl,:) = rip_info(5:7);
            
            x= x+1; rip_trl = rip_trl+1;
            clear rip_info occ_sart occ_end occ_ts occ_ts_4tf occ_ts_gamma occ_ts_gamma_short occ_ts_gamma_baseline hpc_rip_short
        end
    end
    
    %check no nans
    base_nan =cell2mat(cellfun(@find,(cellfun(@isnan,occ_4tf{1,subI}.trial,'UniformOutput',false)),'UniformOutput',false));
    base_inf =cell2mat(cellfun(@find,(cellfun(@isinf,occ_4tf{1,subI}.trial,'UniformOutput',false)),'UniformOutput',false));
    
    if ~isempty(base_nan) || ~isempty(base_inf)
        warning('nan/inf in data!')
    end
    
    % check number of ripples with occ trls <= total ripples
    if size(hpc_occ_sub(subI).unique_rips,1) < invalid_rips(subI,2)
        error('wrong number of ripples in hpc with occ trls!')
    end
    
    occ_4tf{1,subI}.fsample = 500;
    occ_4tf{1,subI}.label = occ_demeaned{1,subI}.label;
    occ_4tf_pre{1,subI}.fsample = 500;
    occ_4tf_pre{1,subI}.label = occ_demeaned{1,subI}.label;
    occ_4tf_post{1,subI}.fsample = 500;
    occ_4tf_post{1,subI}.label = occ_demeaned{1,subI}.label;
    occ_base_4tf{1,subI}.fsample = 500;
    occ_base_4tf{1,subI}.label = occ_demeaned{1,subI}.label;
    occ_base_4tf_pre{1,subI}.fsample = 500;
    occ_base_4tf_pre{1,subI}.label = occ_demeaned{1,subI}.label;
    occ_base_4tf_post{1,subI}.fsample = 500;
    occ_base_4tf_post{1,subI}.label = occ_demeaned{1,subI}.label;
    
    
    l_foi = 35:5:160;%2.5:2.5:32.5;%
    l_tw = 0.2*ones(length(l_foi),1)';
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.keeptrials   = 'yes';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'dpss'; %hanning
    cfg.foi          = l_foi;
    cfg.pad          = 'maxperlen';
    cfg.tapsmofrq    = 10; % specifies half the spectral concentration
    
    cfg.t_ftimwin    = l_tw;
    cfg.toi          = -0.2:0.01:0.2;
    cfg.polyremoval  = 1;
    TF{subI}         = ft_freqanalysis(cfg,occ_4tf{1,subI});
    preTF{subI}      = ft_freqanalysis(cfg,occ_4tf_pre{1,subI});
    postTF{subI}     = ft_freqanalysis(cfg,occ_4tf_post{1,subI});
    baseTF{subI}     = ft_freqanalysis(cfg,occ_base_4tf{1,subI});
    prebaseTF{subI}  = ft_freqanalysis(cfg,occ_base_4tf_pre{1,subI});
    postbaseTF{subI} = ft_freqanalysis(cfg,occ_base_4tf_post{1,subI});
    
    % fit glm with predictors
    cfg              = [];
    cfg.avgoverchan  = 'yes';
    pre_avg_chan{subI} = ft_selectdata(cfg, preTF{subI});
    post_avg_chan{subI} = ft_selectdata(cfg, postTF{subI});
    tf_avg_chan{subI} = ft_selectdata(cfg, TF{subI});
    
    xs = infos{1,subI}.pred(:,1:3);
    for i = 1:size(occ_4tf_pre{subI}.pred,2)
        pre_xs(i,:) = occ_4tf_pre{subI}.pred{i};
    end
    for i = 1:size(occ_4tf_post{subI}.pred,2)
        post_xs(i,:) = occ_4tf_post{subI}.pred{i};
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
    
    TF_avg{subI}.label     = {'occ'};
    preTF_avg{subI}.label     = {'occ'};
    postTF_avg{subI}.label     = {'occ'};
    
    clear occ_gamma_trl hpc_low_trl hpc_idx hpc_rip_trl...
        combined_trl peak_trl trough_trl mean mean_pre mean_post...
        data data_pre data_post
    
end


%% betas groupavg

for subI=1:size(hpc_occ,2)
    null{subI}.powspctrm  = zeros(1,size(TF_avg{1}.freq,2),41);
    null{subI}.time = TF_avg{subI}.time;
    null{subI}.freq = TF_avg{subI}.freq;
    null{subI}.dimord = 'chan_freq_time';
    null{subI}.label = {'occ'};
end
alpha=0.05;nrand='all';latency=[-0.2 0.2];freqi=[35 160];

cfg=[];
cfg.keepindividual = 'yes';
null_grandavg =  ft_freqgrandaverage(cfg, null{:});
ent_avg = ft_freqgrandaverage(cfg,TF_ent{:});
surp_avg = ft_freqgrandaverage(cfg,TF_surp{:});
pre_ent_avg = ft_freqgrandaverage(cfg,preTF_ent{:});
pre_surp_avg = ft_freqgrandaverage(cfg,preTF_surp{:});
post_ent_avg = ft_freqgrandaverage(cfg,postTF_ent{:});
post_surp_avg = ft_freqgrandaverage(cfg,postTF_surp{:});
ent_stats=TF_GroupStats(ent_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Entropy',1, size(hpc_occ_sub,2));
surp_stats=TF_GroupStats(surp_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Surprise',1, size(hpc_occ_sub,2));
pre_ent_stats=TF_GroupStats(pre_ent_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Entropy pre-stim',1, size(hpc_occ_sub,2));
pre_surp_stats=TF_GroupStats(pre_surp_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Surprise pre-stim',1, size(hpc_occ_sub,2));
post_ent_stats=TF_GroupStats(post_ent_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Entropy post-stim',1, size(hpc_occ_sub,2));
post_surp_stats=TF_GroupStats(post_surp_avg,null_grandavg,alpha,alpha,nrand,latency,freqi,'Surprisepost-stim ',1, size(hpc_occ_sub,2));
pre_post_ent_stats=TF_GroupStats(pre_ent_avg,post_ent_avg,alpha,alpha,nrand,latency,freqi,'Entropy pre vs post',1, size(hpc_occ_sub,2));
pre_post_surp_stats=TF_GroupStats(pre_surp_avg,post_surp_avg,alpha,alpha,nrand,latency,freqi,'Surprise pre vs post',1, size(hpc_occ_sub,2));

%% baseline groupavg
cfg =[];
cfg.keepindividual = 'yes';
TF_grandavg = ft_freqgrandaverage(cfg, TF_avg{:});

figure(3);cfg=[]; ft_singleplotTFR(cfg,TF_grandavg);
title({'Grand average','using relative baseline -1.4 to -1'}); set(gca, 'FontSize',14);
ylabel('Frequency (Hz)'); xlabel('Time relative to HPC ripple (S)');
figure(31);
overall_stats=TF_GroupStats(TF_grandavg,null_grandavg,alpha,alpha,nrand,latency,freqi,'rip_occ_overall',1, size(hpc_occ_sub,2));

% TO GET MASKED RELCHANGE FIG NEED TO RMV KEEP INDIVIDUAL
% ADD MASK FIELD WITH POWSCTRM < 0 AND PLOT WITH PARAMETERS
% USED IN TF_GROUPSTATS
cfg =[];
cfg.keepindividual = 'yes';
TF_grandavg_pre = ft_freqgrandaverage(cfg, preTF_avg{:});
cfg=[];
cfg.zlim=[-0.1 0.3];
cfg.maskparameter = 'mask';
%  cfg.maskstyle = 'outline';
% cfg.maskalpha=0.7
%  TF_grandavg_pre.mask = TF_grandavg_pre.powspctrm < 0;
%   cfg.parameter = 'powspctrm';
%   figure('position',[10 10 500 400]);
% h=ft_singleplotTFR(cfg,TF_grandavg_pre);
% title({'Pre stim','using relative baseline -1.4 to -1'});
set(gca, 'FontSize',22);
ylabel('Frequency (Hz)'); xlabel('Time relative to HPC ripple (S)');
pre_stats=TF_GroupStats(TF_grandavg_pre,null_grandavg,alpha,alpha,nrand,latency,freqi,'rip_occ_pre',2, size(hpc_occ_sub,2));

cfg =[];
cfg.keepindividual = 'yes';
TF_grandavg_post = ft_freqgrandaverage(cfg, postTF_avg{:});
figure(5);cfg=[]; ft_singleplotTFR(cfg,TF_grandavg_post)
title({'Post stim','using relative baseline -1.4 to -1'}); set(gca, 'FontSize',14);
ylabel('Frequency (Hz)'); xlabel('Time relative to HPC ripple (S)');
post_stats=TF_GroupStats(TF_grandavg_post,null_grandavg,alpha,alpha,nrand,latency,freqi,'rip_occ_post',3, size(hpc_occ_sub,2));

pre_post =TF_GroupStats(TF_grandavg_pre,TF_grandavg_post,alpha,alpha,nrand,latency,freqi,'rip_occ_pre_post',3, size(hpc_occ_sub,2));


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

for subI = 1:numel(hpc_occ_sub)
    ov = [squeeze(TF_grandavg.powspctrm(subI,:,:,:))];
    ov_mask = ov(overall_mask);
    maxt_beta_sub(subI) = ov_mask(idx_maxt);
end
cohens_d = mean(maxt_beta_sub)/std(maxt_beta_sub);

largest_t_p_d.overall = [maxt,maxp,cohens_d]
clear tvals pvals maxt_beta_sub

post_mask = post_stats.posclusterslabelmat==1;
tvals = post_stats.stat(post_mask);
maxt = max(abs(tvals));
idx_maxt = abs(tvals)==maxt;
pvals = post_stats.prob(post_mask);
maxp = pvals(idx_maxt);

for subI = 1:numel(hpc_occ_sub)
    ovpr = [squeeze(TF_grandavg_pre.powspctrm(subI,:,:,:))];
    ovp_mask = ovpr(post_mask);
    maxt_beta_sub(subI) = ovp_mask(idx_maxt);
end
cohens_d = mean(maxt_beta_sub)/std(maxt_beta_sub);

largest_t_p_d.post = [maxt,maxp,cohens_d]

clear tvals pvals maxt_beta_sub

pre_post_mask = pre_post.negclusterslabelmat==1;
tvals = pre_post.stat(pre_post_mask);
maxt = max(abs(tvals));
idx_maxt = abs(tvals)==maxt;
pvals = pre_post.prob(pre_post_mask);
maxp = pvals(idx_maxt);

for subI = 1:numel(hpc_occ_sub)
    ovpr = [squeeze(TF_grandavg_pre.powspctrm(subI,:,:,:))];
    ovp_mask = ovpr(pre_post_mask);
    ovpo = [squeeze(TF_grandavg_post.powspctrm(subI,:,:,:))];
    ovpo_mask = ovpo(pre_post_mask);
    maxt_beta_sub(subI,1) = ovp_mask(idx_maxt);
    maxt_beta_sub(subI,2) = ovpo_mask(idx_maxt);
end
x1=maxt_beta_sub(:,1);
x2=maxt_beta_sub(:,2);
deltas   = x1 - x2;         % differences
sdDeltas = nanstd(deltas);
cohens_d = (mean(x1)-mean(x2))/sdDeltas;

largest_t_p_d.pre_post = [maxt,maxp,cohens_d];

