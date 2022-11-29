function [statTF] = TF_GroupStats(TF_grndavg,TF_null,alpha,clustalpha,nrand,latency,freqi,nameFig,k, subs)
%run cluter permutation test and plot results
% TF_grandavg is the condition grandavg_freq
% TF_null is the corresponding matrices of 0
% alpha for cfg.alpha (0.025)
% clusteralpha for cfg.clusteralpha (0.05)
% nrand for numrandomization 
% latency is start and end points for stats
% name is string for plot
% k is location of fif in subplot

ns =subs;
if size(latency,2)~= 2
    error('need [begin end]')
end
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.design           = [ones(1,ns), ones(1,ns).*2;[1:ns], [1:ns]];
cfg.uvar             = 2;
cfg.ivar             = 1;
cfg.correctm         = 'cluster';%'no'%
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = alpha;
cfg.clusteralpha     = clustalpha; 
%cfg.clusterstatistic = 'wcm'%'maxsize'; %def maxsum
cfg.minnbchan        = 0;
cfg.correcttail      = 'alpha';
cfg.neighbours       = []; 
cfg.numrandomization = nrand;
%cfg.clusterthreshold = 'nonparametric_individual';
cfg.latency          = latency;
cfg.frequency        = freqi;
statTF = ft_freqstatistics(cfg, TF_grndavg, TF_null);

% %FT tutorial way
% if isfield(statTF, 'posclusters')
% if ~isempty(statTF.posclusters)
% pos_cluster_pvals = [statTF.posclusters(:).prob];
% pos_signif_clust = find(pos_cluster_pvals < alpha);
% pos = ismember(statTF.posclusterslabelmat, pos_signif_clust);
% else
% pos = 0;
% end
% else 
%     pos = 0;
% end
% 
% if isfield(statTF, 'negclusters')
% if ~isempty(statTF.negclusters)
% neg_cluster_pvals = [statTF.negclusters(:).prob];
% neg_signif_clust = find(neg_cluster_pvals < alpha);
% neg = ismember(statTF.negclusterslabelmat, neg_signif_clust);
% else
% neg = 0;
% end
% else
%     neg = 0;
% end
% statTF.mask2 = logical(pos+neg);
% 
% if sum(sum(statTF.mask2,3)) ~= 0
%     statTF.powspctrm = statTF.stat.*statTF.mask2;
% else
%     statTF.powspctrm = statTF.stat;
% end

cfg = [];
Plotclust = statTF;
Plotclust.mask = statTF.mask;%statTF.negclusterslabelmat==1;%;%%
Plotclust.powspctrm = statTF.stat ;
statTF.dimord = 'chan_freq_time';
cfg.parameter = 'powspctrm';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.maskalpha = 1;
cfg.zlim = [-4 4];
% cfg.colormap = colormap(brewermap(256, '*RdBu'));
cfg.figure = 'gcf';
cfg.title = ' ';
%  subplot(2,2,k);

figure('position',[10 10 500 400]);
ft_singleplotTFR(cfg, Plotclust); colorbar; %tit=title(nameFig); 
xlabel('Time relative to HPC ripple (S)'); ylabel('Frequency (Hz)');
set(gca,'FontSize',20);
res=300
% print('-dtiff',fullfile('Manuscript/Figures/',nameFig),['-r' num2str(res)])
end

