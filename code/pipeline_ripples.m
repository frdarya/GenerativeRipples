clearvars -except patient_data nsubs sub stp
rng(111)
region = 'anterior'; baseline = 0; resp_lock = 0;
nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10];
iszurich = logical([zeros(1,13),ones(1,3)]);
% create stp and patient_data structs
[~,stp] = setup(nsubs,iszurich,baseline,resp_lock,region);
patient_data = getMontage(nsubs,stp,0);
hpfilt = 200;
ripdur = 25;
fname = sprintf('HPCRipples/HPCAnterior_vaz_hpf%d_%dms_%dsubjs.mat',hpfilt, ripdur, numel(nsubs));
res = 300;

load(fname)
ent = zeros(numel(nsubs),1);
start = zeros(numel(nsubs),1);
trial_block = zeros(numel(nsubs),1);
times = [-1:0.002:1]';
surp = zeros(1,1);
peak=zeros(numel(nsubs),1);
g=1; j=0.8:0.1:2;x=1; y=1; all_trls = []; all_rips = [];

ytick={'0.8-0.9', '0.9-1','1-1.1','1.1-1.2','1.2-1.3','1.3-1.4' '1.4-1.5',...
    '1.5-1.6','1.6-1.7','1.7-1.8','1.8-1.9','1.9-2'};
xtick={'-1 to -0.8','-0.8 to -0.6', '-0.6 to -0.4', '-0.4 to-0.2',...
    '-0.2 to 0', '0 to 0.2', '0.2 to 0.4','0.4 to 0.6', '0.6 to 0.8','0.8 to 1'};
surp_tick = {'0.5-1','1-1.5','1.5-2','2-2.5','2.5-3','3-3.5','3.5-4','4-4.5','4.5-5'};

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
    if strcmp(stp.region,'anterior')||strcmp(stp.region,'head')||strcmp(stp.region,'body')
        cd(sprintf('hpc_%s',patient_data(subI).hpc_axis))
    else
        cd(sprintf('%s',patient_data(subI).region))
    end
    load('clean_trials_bipolar.mat','clean')
    clean_trials(subI) = clean;
    cd ../../../
    ent_tot(g:g+size(clean_trials(subI).info(:,1),1)-1) = clean_trials(subI).info(:,1);
            % find same ripple in different channels
    sub(subI).unique_rips(:,13) = 1:size(sub(subI).unique_rips,1);
    rips_sorted_by_chan_start{subI} = sortrows(sub(subI).unique_rips,[2,3]);
    dif=diff(rips_sorted_by_chan_start{subI}(:,2:3));
    same_rip =find(abs(dif(:,2))<10 & dif(:,1)==0);
    rips_sorted_by_chan_start{subI}(same_rip,12) = 1;
    % return to original order to remove rip
    rips_sorted_by_chan_start{subI} = sortrows(rips_sorted_by_chan_start{subI},13);
    sub(subI).same_rip=sum(rips_sorted_by_chan_start{subI}(:,12));
    % remove rip
        same_rip_rmv = find(sub(subI).unique_rips(:,12));
    sub(subI).unique_rips(same_rip_rmv,:) = [];
         sub(subI).unique_rips(:,12:13)  = []; 
         
    for ripI = 1:size(sub(subI).unique_rips,1)
        t_idx = sub(subI).unique_rips(ripI,2);
        t = clean_trials(subI).trl(t_idx); % trial # insead of index
        ent(subI,ripI) = clean_trials(subI).info(t_idx,1);
        surp(subI,ripI) = clean_trials(subI).info(t_idx,2);
        start(subI,ripI)=times(sub(subI).unique_rips(ripI,3));
        peak(subI,ripI)=times(sub(subI).unique_rips(ripI,5));
        if t < 40
            trial_block(subI,ripI) = t;
        else
            new_t =  mod(t,40);
            if new_t == 0
                new_t = 40;
            end
            trial_block(subI,ripI) = new_t;
        end
    end
    g=g+size(clean_trials(subI).info(:,1));
    
    
    rip_rate(subI,1:2) = [size(clean_trials(subI).info(:,1),1), size(sub(subI).unique_rips,1)];
    rip_rate(subI,3) = rip_rate(subI,2)/rip_rate(subI,1); %rip per trial
    rip_rate(subI,4) = rip_rate(subI,2)/(rip_rate(subI,1)*2.2); %rip per sec (each trial is 2.2s)
    % ripple frequency as a function of time-bin
    c = histogram(times(sub(subI).unique_rips(:,3)),'BinEdges',[-1:0.2:1]);
    rip_times_prob(subI,1:10) = c.Values ./ sum(c.Values);
    
    % floats for hist
    ent_rounded = 10*round(sub(subI).unique_rips(:,7),1);
    surp_rounded = 10*round(sub(subI).unique_rips(:,8),1);
    e1d= histogram(ent_rounded, 'BinEdges',[8:1:20]);
    subj_count_ent_rips(subI,1:12) = e1d.Values;
    s1d=histogram(surp_rounded, 'BinEdges',[1:5:50]);
    
    subj_count_surp_rips(subI,1:9) = s1d.Values;
    if nsubs(subI) == 22
        subj_ent_bins(subI,:) = zeros(1,12);
        subj_surp_bins(subI,:) = zeros(1,9);
    end
    
    ent_trls = 10*round(clean_trials(subI).info(:,1),1);
    surp_trls = 10*round(clean_trials(subI).info(:,2),1);
    c1d = histogram(ent_trls, 'BinEdges', [8:1:20]);
    count_ent_trls = c1d.Values;
    cs1d  =histogram(surp_trls, 'BinEdges', [1:5:50]);
    count_surp_trls = cs1d.Values;
    subj_tot_ent(subI,:) = count_ent_trls;
    subj_tot_surp(subI,:) = count_surp_trls;
    % can be >1 because of multiple ripples per trial
    % y=1 for each entropy trial there was a correposnding ripple
    subj_ent_bins(subI,:) = subj_count_ent_rips(subI,:) ./ count_ent_trls;
    subj_surp_bins(subI,:) = subj_count_surp_rips(subI,:) ./ count_surp_trls;
    marks = find(isnan(subj_ent_bins(subI,:)));
    
    figure(2);subplot(4,4,subI);
    bar([0.9:0.1:2],subj_ent_bins(subI,:));
    xlabel('Entropy');
    title(sprintf('patient %d',nsubs(subI)));
    ylabel({'ripple conut/entropy trials','>1 because of mult. ripples per trl'});
    if ~isempty(marks); hold on; plot(j(marks),1,'r*'); hold off; end
    clear ent_trls surp_trls count_ent_trls count_surp_trls marks
    
    rip_time = 1000*round(times(sub(subI).unique_rips(:,5)),2);
    ent_rounded = 10*round(sub(subI).unique_rips(:,7),1);
    surp_rounded = 10*round(sub(subI).unique_rips(:,8),1);
    freq_norm_ent = subj_ent_bins(subI,:);
    freq_norm_surp = subj_surp_bins(subI,:);

    cr=histogram(rip_time,'BinEdges',[-1000:200:1000])
    time_dist(subI,:) = cr.Values;
    ce = histogram(ent_rounded,'BinEdges',[8:1:20])
        ent_dist(subI,:) = ce.Values;

    figure(3);sgtitle({sprintf('Patient %d',nsubs(subI)),'Ripple peak time/entropy distribution'});
    subplot(121);
    c=histogram2(rip_time,ent_rounded,'XBinEdges',[-1000:200:1000],'YBinEdges',[8:1:20]); % if want prob add 'Normalization','Probability', otherwise it's normalized count
    xlabel('Ripple peak time'); ylabel('Entropy'); zlabel('Count');title('Raw count')
    ent_time_hist = c.Values;
    cs=histogram2(rip_time,surp_rounded,'XBinEdges',[-1000:200:1000],'YBinEdges',[1:5:50]); % if want prob add 'Normalization','Probability', otherwise it's normalized count
    xlabel('Ripple peak time'); ylabel('Surprise'); zlabel('Count');title('Raw count')
    surp_time_hist = cs.Values;
    %     prop_ent_time_hist = ent_time_hist ./ sum(sum(ent_time_hist)); % convert to prob to avoid issues with different numbers of trials
    sub(subI).ent_prop_time = ent_time_hist .* freq_norm_ent;
    sub(subI).surp_prop_time = surp_time_hist .* freq_norm_surp;
    sub(subI).ent_time = ent_time_hist;
    sub(subI).surp_time = surp_time_hist;
    subplot(122); h=heatmap(sub(subI).ent_prop_time');
    ylabel('Entropy'); xlabel('Ripple peak time');
    h.XDisplayLabels = xtick; h.YDisplayLabels = ytick;
    h.CellLabelColor='none'; h.NodeChildren(3).YDir='normal';
    title('Normalised by proportion per entropy bin')

    ent_rip_prop_3d(:,:,subI) = sub(subI).ent_prop_time;
    ent_rip_3d(:,:,subI) = sub(subI).ent_time;
    
    surp_rip_prop_3d(:,:,subI) = sub(subI).surp_prop_time;
    surp_rip_3d(:,:,subI) = sub(subI).surp_time;
    clear rip_time rip_ent rip_surp freq_norm_ent freq_norm_surp ent_time_hist surp_time_hist
    
    % data for R - LME (get trials without ripples)
    t_indices = sub(subI).unique_rips(:,2);
    no_rips = setdiff(1:size(clean_trials(subI).info,1),t_indices);
    no_rips_mat = clean_trials(subI).info(no_rips,1:5);
    no_rips_mat = [nan(size(no_rips,2),6), no_rips_mat];
    no_rips_mat(:,2) = no_rips';
    all_rips(x:x+size(sub(subI).unique_rips,1)-1,:) = [sub(subI).unique_rips,...
        repmat(subI,[size(sub(subI).unique_rips,1),1])];
    
    all_len = size(sub(subI).unique_rips,1)+size(no_rips_mat,1);
    all_trls(y:y+all_len-1,:) = [[no_rips_mat;sub(subI).unique_rips],...
        repmat(subI,all_len,1)];
    
    x = x+size(sub(subI).unique_rips,1);
    y = y+all_len;
    clear no_rips_mat
end

%% time-ent heatmap

if numel(nsubs) > 7
    ent_rip_prop_3d(:,:,8) = [];
    surp_rip_prop_3d(:,:,8) = [];
end
%normalised
group_avg_ent_time = nanmean(ent_rip_prop_3d,3);
group_sum_ent_time = nansum(ent_rip_prop_3d,3);
figure( 'position',[10 10 500 400]); h1=heatmap(group_avg_ent_time'); %title({'Group average','Normalised by entropy bin'})
ylabel('Entropy'); xlabel('Ripple peak time');
h1.XDisplayLabels = xtick; h1.YDisplayLabels = ytick; h1.FontSize = 20;
h1.CellLabelColor='none'; h1.NodeChildren(3).YDir='normal';h1.Colormap = parula;
 print('-dsvg',fullfile('Manuscript/Figures/','normalised_rip_ent_time'),['-r' num2str(res)])

% raw
group_avg_ent_time_r = nanmean(ent_rip_3d,3);
figure( 'position',[10 10 900 700]); h1=heatmap(group_avg_ent_time_r'); %title({'Group average'})
ylabel('Entropy'); xlabel('Ripple peak time');
h1.XDisplayLabels = xtick; h1.YDisplayLabels = ytick; h1.FontSize = 20;
h1.CellLabelColor='none'; h1.NodeChildren(3).YDir='normal';h1.Colormap = parula;
% print('-dtiff',fullfile('Manuscript/Figures/','raw_rip_ent_time'),['-r' num2str(res)])

group_avg_surp_time = nanmean(surp_rip_prop_3d,3);
group_sum_surp_time = nansum(surp_rip_prop_3d,3);
figure( 'position',[10 10 900 700]); h1=heatmap(group_avg_surp_time'); %title({'Group average','Normalised by surprise bin'})
ylabel('Surprise'); xlabel('Ripple peak time');
h1.XDisplayLabels = xtick; h1.YDisplayLabels = surp_tick; h1.FontSize = 20;
h1.CellLabelColor='none'; h1.NodeChildren(3).YDir='normal';h1.Colormap = parula;
% print('-dtiff',fullfile('Manuscript/Figures/','normalised_rip_surp_time'),['-r' num2str(res)])

group_avg_surp_time_r = nanmean(surp_rip_3d,3);
figure( 'position',[10 10 900 700]); h1=heatmap(group_avg_surp_time_r'); %title({'Group average'})
ylabel('Surprise'); xlabel('Ripple peak time');
h1.XDisplayLabels = xtick; h1.YDisplayLabels = surp_tick; h1.FontSize = 20;
h1.CellLabelColor='none'; h1.NodeChildren(3).YDir='normal';h1.Colormap = parula;
% print('-dtiff',fullfile('Manuscript/Figures/','raw_rip_surp_time'),['-r' num2str(res)])

% perm test
if numel(nsubs) > 7
    nperm = numel(nsubs)-1;
else
    nperm = numel(nsubs);
end

for subI = 1:nperm
    sub_mat=ent_rip_prop_3d(:,:,subI);
    sub_mat = sub_mat(~isnan(sub_mat));
    
    nrand = 1000;
    for i = 1:nrand
        shuff_mat(:,i) = sub_mat(randperm(size(sub_mat,1)));
        corrs(subI,i) = corr(sub_mat,shuff_mat(:,i),'Type','Spearman');
        corrs_kend(subI,i) = corr(sub_mat,shuff_mat(:,i),'Type','Kendall');

    end
    clear sub_mat shuff_mat
end

% not sig meaning does not correlate with random noise.
avg_corr = mean(corrs,2);
[h_cor,p_cor,ci,stats_cor]=ttest(avg_corr);

% kendall - same result as spearman
avg_corr_kend = mean(corrs_kend,2);
[h_cor_k,p_cor_k,ci_k,stats_cor_k]=ttest(avg_corr_kend);

% % or x2 on the vectorised avg matrix
% [h_x2,p_x2,stats_x2]=chi2gof(group_avg_ent_time(:));

% ks on the avg matrix - sig --> not normally distributed
[h_ks,p_ks,ksstat,cv]=kstest(group_avg_ent_time);

%uniform distribution - sig --> not uniformly distributed
% https://math.stackexchange.com/questions/2435/is-there-a-simple-test-for-uniform-distributions
% dist=makedist('uniform',0,4);
% [h_ks_uni,p_ks_uni,ksstat_uni,cv_uni]=kstest(group_avg_ent_time,dist);

Xent = unifrnd(min(group_avg_ent_time(:)),max(group_avg_ent_time(:)),10,12);
[~,ent_p_ks_unif,ent_ksstat_uni]=kstest2(group_avg_ent_time(:),Xent(:));

Xsurp = unifrnd(min(group_avg_surp_time(:)),max(group_avg_surp_time(:)),10,12);
[~,surp_p_ks_unif,surp_ksstat_uni]=kstest2(group_avg_surp_time(:),Xsurp(:));

%% trial in block
idx=find(trial_block);
rips_trial_block = trial_block(find(trial_block));
rips_ent = ent(find(ent));
rips_surp = surp(idx);
rips_peak = peak(idx);

ytick={'< 0.8','0.8-1', '1-1.2','1.2-1.4','1.4-1.6','1.6-1.8','1.8-2'};
xtick={'-1 to -0.8','-0.8 to -0.6', '-0.6 to -0.4', '-0.4 to-0.2',...
    '-0.2 to 0', '0 to 0.2', '0.2 to 0.4','0.4 to 0.6', '0.6 to 0.8','0.8 to 1'};
% three way heatmap
x=discretize(rips_peak,10);
y=discretize(rips_ent,7);
tbl=array2table([x,y,rips_trial_block]);
figure(8);h=heatmap(tbl,'Var1','Var2','ColorVariable','Var3','ColorMethod','mean');
xlabel('Time'); ylabel('Entropy');title('Distribution of HPC ripples as a function of entropy, peri-stimulus time and trial # in block');
set(gca, 'FontSize', 20);
h.XDisplayLabels = xtick;
h.YDisplayLabels = ytick;
h.CellLabelColor='none'; h.NodeChildren(3).YDir='normal';

% distribution over trials in block
for subI = 1:numel(nsubs)
    if nsubs(subI)==22
        continue
    end
    trlIdx= find(trial_block(subI,:));
    blkI = trial_block(subI,trlIdx);
    c1d = histogram(blkI, 'BinEdges', [1:41],'Normalization', 'Probability');%,   
    count_ent_trls(subI,:) = c1d.Values;
        clean_trials(subI).info(:,6)=mod(clean_trials(subI).info(:,5),40);
    clean_trials(subI).info((clean_trials(subI).info(:,6)==0),6) = 40;
    trl_idx=clean_trials(subI).info(:,6)==1;%missing first trial
    if sum(trl_idx) ==0
        count_ent_trls(subI,1)=NaN;
    end

    
end
% fit exponential learning curve
% https://people.richland.edu/james/lecture/m116/logs/models.html
x = [1:40]';
g = fittype('b*(1-exp(-c*x))');
[f_exp, gof]=fit(x,nanmean(count_ent_trls)',g,'StartPoint',[1,0]);
figure( 'position',[10 10 900 700]);plot(f_exp,x,nanmean(count_ent_trls)');
xlabel('Trial # of block'); ylabel('p(ripple)');text(31,0.03,sprintf('adjusted R^2 = %.2f',gof.adjrsquare))
set(gca, 'FontSize', 20,'LineWidth',2);


[f_lin, gof_lin] = fit(x,nanmean(count_ent_trls)','poly2');
figure(10);plot(f_lin,x,nanmean(count_ent_trls)');
xlabel('Trial # of block'); ylabel('p(ripple)');
text(31,0.03,sprintf('adjusted R^2 = %.2f',gof_lin.adjrsquare))
set(gca, 'FontSize', 20);

err = nanstd(count_ent_trls)/sqrt(numel(nsubs)-1);
figure('position',[10 10 500 400]);bar(nanmean(count_ent_trls)); 
hold on;
er = errorbar([1:40],nanmean(count_ent_trls),err);
er.Color = [0 0 0];
er.LineStyle = 'none';
p=plot(f_exp,'-r');xlabel('Trial # in block');ylabel('Ripple probability'); ylim([0 0.045]);
p.LineWidth=2;hold off
set(gca, 'FontSize', 20);
print('-dtiff',fullfile('Manuscript/Figures/','rip_trl_in_blk'),['-r' num2str(res)])

% text(1,0.042,sprintf('adjusted R^2 = %.2f',gof.adjrsquare))hold off;

