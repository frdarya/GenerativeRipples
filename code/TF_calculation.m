function [TF_Low_all, TF_High_all] = TF_calculation(nsubs,stp,patient_data,rejectedTrials,plot)
%Preprocess data and compute TF, outputs TF for low, high freqs separately.
% input arguments:
%   nsubs = array of patient numbers
%   stp = output from setup function
%   patient_data = output from setup function (hpcContacts function)
%   rejectedTrials = structure from artefact (time or tf) cleaning which
%   includes a keeptrials field with trl index to use.
%   plot = whether or not to plot TF per subj (1/0)

h_foi = 35:2.5:160; % Frequency of interest
l_foi = 2.5:2.5:32.5;
h_tw  = 0.4*ones(length(h_foi),1)'; % fRayleigh = 1/0.4s = 2.5 Hz
l_tw = 0.4*ones(length(l_foi),1)';
fw  = 10*ones(length(h_foi),1)';  %frequency smoothing - only for high frequencies


for subI = 1:numel(nsubs)
    fprintf(['Calculating TF for Patient ',num2str(nsubs(subI)), '\n'])
    keeptrials = rejectedTrials(subI).info; % to get trl num 1:480, and not index from cleanbehav
    cleandata = preproc(subI,nsubs,stp,patient_data,keeptrials);    
    % TF
    if stp.zurich(subI) == 1 || nsubs(subI) == 31 || nsubs(subI) == 32 ... 
            || nsubs(subI) == 36 || nsubs(subI) == 37
        % downsample - for zurich and p31 & 36 ruber
        cfg             = [];
        cfg.resamplefs  = 500;
        cfg.demean      = 'yes';
        cleandata       = ft_resampledata(cfg,cleandata);
    end
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.keeptrials   = 'yes';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'dpss';
    cfg.pad          = 'maxperlen';
    cfg.foi          = h_foi; % must contain integer multiples of the Rayleigh frequency and the spectral concentration
    cfg.tapsmofrq    = fw(1); % specifies half the spectral concentration
    cfg.toi          = -stp.lat:0.01:stp.lat;
    cfg.t_ftimwin    = h_tw; % overlap between time windows
    cfg.polyremoval  = 1;
    TF_High          = ft_freqanalysis(cfg, cleandata);
    
    % Low
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.keeptrials   = 'yes';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = l_foi;
    cfg.pad          = 'maxperlen';
    cfg.t_ftimwin    = l_tw;
    cfg.toi          = -stp.lat:0.01:stp.lat;
    cfg.polyremoval  = 1;
    TF_Low           = ft_freqanalysis(cfg,cleandata);
    
    if nsubs(subI) == 25 || nsubs(subI) == 36 || nsubs(subI) == 37
        TF_Low.trialinfo = zeros(size(TF_Low.powspctrm,1),1);
        TF_High.trialinfo = zeros(size(TF_High.powspctrm,1),1);
    end
    
    TF_Low_all(subI) = TF_Low;
    TF_High_all(subI) = TF_High;
    
    if plot
        h = figure();
        subplot(1,2,1);
        cfg = [];
        cfg.xlim =[-1 1];
        cfg.ylim = [2.5 32.5];
        cfg.colorbar = 'yes';
        ft_singleplotTFR(cfg,TF_Low)
        subplot(1,2,2);
        cfg = [];
        cfg.xlim =[-1 1];
        cfg.ylim = [35 150];
        cfg.colorbar = 'yes';
        ft_singleplotTFR(cfg,TF_High)
        pause(1)
        close(h)
    end
    clearvars -except nsubs h_foi l_foi fw lat l_tw h_tw plot...
        TF_Low_all TF_High_all patient_data iszurich stp rejectedTrials
end
end