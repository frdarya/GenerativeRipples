%% STEP 3 of analysis - analyse iEEG data - HPC TF
% clean and preprocess data, construct TF and GLM, output group stats
% must run section 1) setup
% for cleaning  run sections 2-4
% for analysis and stats run sections 5-7
%
% Darya Frank 21-10-2020

%% 1) setup
clearvars
rng(111,'twister')
global ft_default
ft_default.spmversion = 'spm12';
ft_defaults
ft_hastoolbox('brewermap', 1);

region = 'anterior'; resp_lock = 0; baseline=0;
% for cleaning
if strcmpi(region,'head')
    nsubs = [3,9,13,15,16,22,25,36,37,6,10]; % 22 is bad
    iszurich = logical([zeros(1,9),ones(1,2)]);
elseif strcmpi(region,'body')
    nsubs = [3,6,8,13,31,32,6,8,10]; 
    iszurich = logical([zeros(1,6),ones(1,3)]);
elseif strcmpi(region,'tail')
    nsubs = [3,9,27,10];
    iszurich = logical([zeros(1,3),ones(1,1)]);
    % for analysis
elseif strcmpi(region,'anterior')
    nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10]; %22 not used poor data
    iszurich = logical([zeros(1,13),ones(1,3)]);
else
    error('wrong hippocampus axis spec')
end

% create stp and patient_data structs
[patient_data,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);

%% 2) choose behav trials and time-artf reject
% get valid behhavioural trials, load ieeg data and do artifact rejection
% for time domain.

% saves cleanBehavTime in each patient's folder (skips if folder exists)
[art_behav_time] = behav_time_clean(nsubs,patient_data,stp);

%% 3) Frequency domain artefact rejection
clear art_behav_time
% if need to reconstruct art_behav_time
%clear art_behav_time
% for subI = 1:numel(nsubs)
%     if stp.zurich(subI)
%         foldn = sprintf('Information/P%dz',nsubs(subI));
%     else
%         foldn = sprintf('Information/Patient%d+',nsubs(subI));
%     end
%     cd(foldn)
%     cd(sprintf('hpc_%s',stp.region))
%     load(sprintf('artf_behav_time_%s.mat',stp.montageD))
%     if isfield(art,'rmvtrials')
%         art = rmfield(art, 'rmvtrials');
%     end
%     art_behav_time(subI) = art;
%     cd ../../../
% end

% select only bipolar montage (cleaned artifacts already)
clear patient_data
patient_data = getMontage(nsubs,stp,0);

% calculate TF - high & low freqs
[TF_Low_all, TF_High_all] = TF_calculation(nsubs,stp,patient_data,art_behav_time,0);

% reject trials based on freq artifacts - saves artf_freq in each
% patient's folder - trial number represents trial number from 1:480.
[~, art_freq] = freq_clean(nsubs,stp,TF_High_all,TF_Low_all);

%% 4) create union of Behaviour, Time and TF clean trials.
% run if structs for all patients exist

for subI = 1:numel(nsubs)
    fprintf(['Starting Patient ',num2str(nsubs(subI)), '\n'])
    if stp.zurich(subI)
        foldn = sprintf('Information/P%dz',nsubs(subI));
    else
        foldn = sprintf('Information/Patient%d+',nsubs(subI));
    end
    cd(foldn)
    if strcmpi(stp.region,'anterior')
        cd(sprintf('hpc_%s',patient_data(subI).hpc_axis))
    else
        cd(sprintf('hpc_%s',stp.region))
    end
    % uncomment if need to reconstruct individually
    load(sprintf('artf_behav_time_%s.mat',stp.montageD))
    if isfield(art,'rmvtrials')
        art = rmfield(art, 'rmvtrials');
    end
    art_behav_time(subI) = art;
    
    load(sprintf('artf_freq_%s.mat',stp.montageD))
    art_freq(subI) = art;
    trl_num = art_behav_time(subI).info(:,5); % trial to keep after bahv & time cleaning
    if nsubs(subI) == 37 && trl_num(end) == 480 % missing first trigger
        trl_num = trl_num-1;
    end
    idx = ismember(trl_num,art_freq(subI).keeptrials);
    
    clean_trials(subI).nsub = nsubs(subI);
    clean_trials(subI).trl = art_freq(subI).keeptrials;
    clean_trials(subI).info = art_behav_time(subI).info(idx,:);
    
    clean = clean_trials(subI);
    save(sprintf('clean_trials_%s.mat',stp.montageD),'clean')
    
    cd ../../../
end

%% 5) Compute TF
% get bipolar montage
clear patient_data
% if strcmpi(stp.region,'anterior')
patient_data = getMontage(nsubs,stp,0);
% else
%     error('trying to run analysis but not with  HPC anterior!')
% end
% get clean trials
for subI = 1:numel(nsubs)
    fprintf(['Starting Patient ',num2str(nsubs(subI)), '\n'])
    if stp.zurich(subI) == 1
        foldn = sprintf('Information/P%dz',nsubs(subI));
    else
        foldn = sprintf('Information/Patient%d+',nsubs(subI));
    end
    cd(foldn)
    cd(sprintf('hpc_%s',patient_data(subI).hpc_axis))
    
    load('clean_trials_bipolar.mat','clean')
    clean_trials(subI) = clean;
    
    cd ../../../
end

% get TF
[TF_Low_all, TF_High_all] = TF_calculation(nsubs,stp,patient_data, clean_trials,0);

% Combine freqs
for subI = 1:numel(nsubs)
    TF_Low            = TF_Low_all(subI);
    TF_High           = TF_High_all(subI);
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.appendim      = 'freq';
    TF_Combined(subI) = ft_appendfreq(cfg,TF_Low,TF_High);
    
    if size(TF_Combined(subI).powspctrm,2) ~= 1 %channels
        TF_Combined(subI).powspctrm = nanmean(TF_Combined(subI).powspctrm,2);
    end
    
    TF_Combined(subI).powspctrm = squeeze(TF_Combined(subI).powspctrm);
    TF_Combined(subI).label = {stp.region};
end

%% 6) TF GLM
% subs to select from
subs = nsubs;
k=1;
stp.add_pred = 0;
% data transformation/GLM specification
stp.logcorr = 1;
stp.run_glm = 1; 

[TF_ent,TF_surp,TF_meanent,TF_const,TF_ent_grandavg,TF_surp_grandavg,TF_meanent_grandavg,TF_null_grandavg]...
    = TF_GLM(subs,stp,k,TF_Combined,clean_trials);

%% 7) Group stats
% rng(111)

clustalpha = 0.05; alpha = 0.05; numrand = 5000;
lat = 1.2;
if stp.resp_lock == 1
    prelat = [-0.1 0.4];
else
    prelat = [-1 0];
    postlat = [0 1];
end

% subset patients
nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];

hicks =[1:7,9:16];%patient 22 has less than 25% of trials
%
if stp.run_glm==0
    TF_comb.powspctrm = TF_comb.powspctrm(hicks,:,:,:);
    TF_null_grandavg.powspctrm = TF_null_grandavg.powspctrm(hicks,:,:,:);
    
    stats.all_pre_low = TF_GroupStats(TF_comb,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[2.5 32.5],sprintf('%s_pre_low_alltrls',stp.region),1, numel(hicks));
    
    stats.all_pre_high = TF_GroupStats(TF_comb,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[35 160],sprintf('%s_pre_high_alltrls',stp.region),1, numel(hicks));
    
    stats.all_post_low = TF_GroupStats(TF_comb,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[2.5 32.5],sprintf('%s_post_low_alltrls',stp.region),1, numel(hicks));
    
    stats.all_post_high = TF_GroupStats(TF_comb,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[35 160],sprintf('%s_post_high_alltrls',stp.region),1, numel(hicks));
    
    cfg =[];cfg.xlim=[-1 1]; cfg.ylim = [35 160];figure();ft_singleplotTFR(cfg,TF_comb)
else
    TF_ent_grandavg.powspctrm = TF_ent_grandavg.powspctrm(hicks,:,:,:);
    TF_surp_grandavg.powspctrm = TF_surp_grandavg.powspctrm(hicks,:,:,:);
    TF_meanent_grandavg.powspctrm = TF_meanent_grandavg.powspctrm(hicks,:,:,:);
    TF_null_grandavg.powspctrm = TF_null_grandavg.powspctrm(hicks,:,:,:);
end

if stp.resp_lock == 1
    % Entropy & Surprise
    stats.Entropy_high_resp = TF_GroupStats(TF_ent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[35 160],sprintf('%s_resp_entropy_high',stp.region),1, numel(hicks));
    stats.Entropy_low_resp = TF_GroupStats(TF_ent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[2.5 32.5],sprintf('%s_resp_entropy_low',stp.region),3, numel(hicks));
    stats.Surprise_high_resp = TF_GroupStats(TF_surp_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[35 160],sprintf('%s_resp_surprise_high',stp.region),2, numel(hicks));
    stats.Surprise_low_resp = TF_GroupStats(TF_surp_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[2.5 32.5],sprintf('%s_resp_surprise_low',stp.region),4, numel(hicks));
    
else
    % Entropy
    stats.Entropy_high_pre = TF_GroupStats(TF_ent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[35 160],sprintf('%s_pre_entropy_high',stp.region),1, numel(hicks));
    stats.Entropy_high_post = TF_GroupStats(TF_ent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[35 160],sprintf('%s_post_entropy_high',stp.region),2, numel(hicks));
    stats.Entropy_low_pre = TF_GroupStats(TF_ent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[2.5 32.5],sprintf('%s_pre_entropy_low',stp.region),3, numel(hicks));
    stats.Entropy_low_post = TF_GroupStats(TF_ent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[2.5 32.5],sprintf('%s_post_entropy_low',stp.region),4, numel(hicks));
    
    % Surprise
    stats.Surprise_high_pre = TF_GroupStats(TF_surp_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[35 160],sprintf('%s_pre_surprise_high',stp.region),1, numel(hicks));
    stats.Surprise_high_post = TF_GroupStats(TF_surp_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[35 160],sprintf('%s_post_surprise_high',stp.region),2, numel(hicks));
    stats.Surprise_low_pre = TF_GroupStats(TF_surp_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[2.5 32.5],sprintf('%s_pre_surprise_low',stp.region),3, numel(hicks));
    stats.Surprise_low_post = TF_GroupStats(TF_surp_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[2.5 32.5],sprintf('%s_post_surprise_low',stp.region),4, numel(hicks));
    
    % Mean Entropy
    stats.MeanEntropy_high_pre = TF_GroupStats(TF_meanent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[35 160],sprintf('%s_pre_meanent_high',stp.region),1, numel(hicks));
    stats.MeanEntropy_high_post = TF_GroupStats(TF_meanent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[35 160],sprintf('%s_post_meanent_high',stp.region),2, numel(hicks));
    stats.MeanEntropy_low_pre = TF_GroupStats(TF_meanent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,prelat,[2.5 32.5],sprintf('%s_pre_meanent_low',stp.region),3, numel(hicks));
    stats.MeanEntropy_low_post = TF_GroupStats(TF_meanent_grandavg,TF_null_grandavg,...
        alpha,clustalpha,numrand,postlat,[2.5 32.5],sprintf('%s_post_meanent_low',stp.region),4, numel(hicks));
end


%% largest t, pval, effect size
% ent
clear tvals pvals maxt_beta_sub

ent_high_pre_mask = stats.Entropy_high_pre.negclusterslabelmat==1;
tvals = stats.Entropy_high_pre.stat(ent_high_pre_mask);
maxt = max(abs(tvals));
idx_maxt = abs(tvals)==maxt;
pvals = stats.Entropy_high_pre.prob(ent_high_pre_mask);
maxp = pvals(idx_maxt);

for subI = 1:numel(hicks)
    ent_pre = [squeeze(TF_ent_grandavg.powspctrm(subI,:,14:end,1:101))];
    ent_pre_mask = ent_pre(ent_high_pre_mask);
    maxt_beta_sub(subI) = ent_pre_mask(idx_maxt);
end
cohens_d = mean(maxt_beta_sub)/std(maxt_beta_sub);

largest_t_p_d.ent_high_pre = [maxt,maxp,cohens_d];
clear tvals pvals maxt_beta_sub
% surp
surp_high_post_mask = stats.Surprise_high_post.posclusterslabelmat==1;
tvals = stats.Surprise_high_post.stat(surp_high_post_mask);
maxt = max(abs(tvals));
idx_maxt = abs(tvals)==maxt;
pvals = stats.Surprise_high_post.prob(surp_high_post_mask);
maxp = pvals(idx_maxt);

for subI = 1:numel(hicks)
    surp_post = [squeeze(TF_surp_grandavg.powspctrm(subI,:,14:end,101:end))];
    surp_post_mask = surp_post(surp_high_post_mask);
    maxt_beta_sub(subI) = surp_post_mask(idx_maxt);
end
cohens_d = mean(maxt_beta_sub)/std(maxt_beta_sub);

largest_t_p_d.surp_high_post = [maxt,maxp,cohens_d]
save(sprintf('TF_LFP/%s_largest_t_p_d.mat',stp.region),'largest_t_p_d')

%% plot betas
ent_high_pre_mask = stats.Entropy_high_pre.negclusterslabelmat==1;
for subI = 1:numel(nsubs)
    surp_high_pre_beta = TF_surp{1,subI}.powspctrm(:,14:end,1:101);
    surp_beta(subI,:) = mean2(surp_high_pre_beta(ent_high_pre_mask));
    ent_high_pre_beta = TF_ent{1,subI}.powspctrm(:,14:end,1:101);
    ent_beta(subI,:) = mean2(ent_high_pre_beta(ent_high_pre_mask));
    meanEnt_high_post_beta = TF_meanent{1,subI}.powspctrm(:,14:end,1:101);
    meanEnt_beta(subI,:) = mean2(meanEnt_high_post_beta(ent_high_pre_mask));
end

figure();
beeswarm([ones(15,1);ones(15,1).*2],...
    [ent_beta(hicks);surp_beta(hicks)],'sort_style','up',...
    'dot_size',2,'overlay_style','ci','use_current_axes',false)
xticklabels({'Entropy','Surprise'});xticks([1,2])

%%
surp_high_post_mask = stats.Surprise_high_post.posclusterslabelmat==1;
for subI = 1:numel(nsubs)
    surp_high_pre_beta = TF_surp{1,subI}.powspctrm(:,14:end,101:end);
    surp_beta(subI,:) = mean2(surp_high_pre_beta(surp_high_post_mask));
    ent_high_pre_beta = TF_ent{1,subI}.powspctrm(:,14:end,101:end);
    ent_beta(subI,:) = mean2(ent_high_pre_beta(surp_high_post_mask));
    meanEnt_high_post_beta = TF_meanent{1,subI}.powspctrm(:,14:end,101:end);
    meanEnt_beta(subI,:) = mean2(meanEnt_high_post_beta(surp_high_post_mask));
end
figure();
beeswarm([ones(15,1);ones(15,1).*2],...
    [ent_beta(hicks);surp_beta(hicks)],'sort_style','up',...
    'dot_size',2,'overlay_style','ci','use_current_axes',false)
xticklabels({'Entropy','Surprise'});xticks([1,2])

%% power as trial in block
ent_high_pre_mask = stats.Entropy_high_pre.negclusterslabelmat==1;
k=1;
for subI = 1:numel(nsubs)
    if nsubs(subI) == 22
        continue
    end
    % get betas to correlate later with occ
    surp_high_pre_beta = TF_ent{1,subI}.powspctrm(:,14:end,1:101);
    ent_beta(subI,:) = mean2(surp_high_pre_beta(ent_high_pre_mask));
    
    % get power to correlate now with entropy
    ent_high_pre_power = log(TF_High_all(subI).powspctrm(:,:,:,21:121));
    for trlI = 1:size(ent_high_pre_power,1)
        for chanI = 1:size(ent_high_pre_power,2)
            trl_ent_high_pre_power = ent_high_pre_power(trlI,chanI, :,:);
            trl_power(chanI,:) = trl_ent_high_pre_power(ent_high_pre_mask);
        end
        mean_trl_power = mean2(trl_power);
        preds = clean_trials(subI).info(trlI,:);
        trlNum = clean_trials(subI).trl(trlI);
        if trlNum < 40
            trial_block = trlNum;
        else
            new_t =  mod(trlNum,40);
            if new_t == 0
                new_t = 40;
            end
            trial_block = new_t;
        end
        to_correlate(trlI,:) = [mean_trl_power, preds];
        power_trlBlk(trlI,1:2) = [mean_trl_power, trial_block];
        power_trlBlk_sorted(trlI,trial_block) = mean_trl_power;
        all_subs_trl_power(k,1:3) = [mean_trl_power,trlNum,nsubs(subI)];
        k=k+1;
        clear preds mean_trl_power trl_power trl_ent_high_pre_power trial_block
    end
    display(['Patient ' num2str(nsubs(subI)), ' number of trials: ' num2str(trlI)])
    ent_high_coefficients1 = polyfit(to_correlate(:,2),to_correlate(:,1), 1);
    ent_high_coefficients2 = polyfit(to_correlate(:,2),to_correlate(:,1), 2);
    
    power_trlBlk_sorted(find(power_trlBlk_sorted==0)) = NaN;
    power_trlBlk_sorted_sub(subI,1:40) = nanmean(power_trlBlk_sorted);
    % Create a new x axis
    xFit = linspace(min(to_correlate(:,2)), max(to_correlate(:,2)), 10);
    % Get the estimated yFit value for each of new x locations.
    yFit1 = polyval(ent_high_coefficients1 , xFit);
    yFit2 = polyval(ent_high_coefficients2 , xFit);
    
    figure(123);subplot(4,4,subI);scatter(to_correlate(:,2),to_correlate(:,1));
    hold on;
    plot(xFit, yFit1, 'r--', 'LineWidth', 2);
    plot(xFit,yFit2,'k--','LineWidth', 2);
    xlabel('Entropy'); ylabel('log(power) of gamma cluster');
    hold off
    ent_correlation = corrcoef(to_correlate(:,1:2));
    power_ent_corr(subI) = ent_correlation(1,2);
    
    trl_high_coefficients1 = polyfit(power_trlBlk(:,2),power_trlBlk(:,1), 1);
    trl_high_coefficients2 = polyfit(power_trlBlk(:,2),power_trlBlk(:,1), 2);
    xFit = linspace(min(power_trlBlk(:,2)), max(power_trlBlk(:,2)), 10);
    % Get the estimated yFit value for each of new x locations.
    yFit1 = polyval(trl_high_coefficients1 , xFit);
    yFit2 = polyval(trl_high_coefficients2 , xFit);
    figure(999);subplot(4,4,subI);scatter(power_trlBlk(:,2),power_trlBlk(:,1));
    hold on
    plot(xFit, yFit1, 'r--', 'LineWidth', 2);
    plot(xFit,yFit2,'k--','LineWidth', 2);
    xlabel('Trial # in block'); ylabel('log(power) of gamma cluster');
    hold off
    trl_power_correlation = corrcoef(power_trlBlk(:,1:2));
    power_trl_corr(subI) = trl_power_correlation(1,2);
    yFits(subI,:) = yFit2;
    clear to_correlate power_trlBlk power_trlBlk_sorted
end

figure();
bounds = nanstd(power_trlBlk_sorted_sub)/sqrt(15);
boundedline(1:40, nanmean(power_trlBlk_sorted_sub), bounds, 'alpha', 'transparency', 0.1);
xlabel('Trial # in block'); ylabel('log(power) of gamma cluster');
title('Entropy cluster gamma power as a function of trial in block');
set(gca, 'FontSize', 14);

figure();beeswarm(ones(15,1),ent_beta(hicks),'sort_style','up',...
    'corral_style','gutter','dot_size',5,'overlay_style','ci',...
    'use_current_axes',false)

