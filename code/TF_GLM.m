function varargout = TF_GLM(nsubs,stp,k,TF_Combined,clean_T_TF,run_glm)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5
    run_glm = 1;
end

for subI = 1:numel(nsubs)
    
    fprintf(['Starting GLM Patient ',num2str(nsubs(subI)), '\n'])
    
    % remove NaNs from TF
    TF_Combined(subI).powspctrm = TF_Combined(subI).powspctrm(:,:,21:end-21);
    TF_Combined(subI).time = TF_Combined(subI).time(21:end-21);
    TF_Combined(subI).dimord = 'rpt_freq_time';
    
    % remove trials rejected during TF cleaning
    cleanbehav = clean_T_TF(subI).info;
    timepoints = size(TF_Combined(subI).powspctrm,3);
    freqs = size(TF_Combined(subI).powspctrm,2);
    trls = size(TF_Combined(subI).powspctrm,1);
    
    % check correct number of trials was removed from behaviour
    if size(cleanbehav,1) ~= trls
        error('mismatch between cleaned behavior trials and TF trials!')
    end
    
    
    % de-mean or logcorr
    if stp.demean == 1 % from July 19 demean per timepoint (following dlz afvice),
        % not per frequency as before because then we cannot compare between freqs
        for timeI = 1:size(TF_Combined(subI).powspctrm,3)
            mean_power(1:trls,timeI) = mean(TF_Combined(subI).powspctrm(:,:,timeI),2);
        end
        pow = permute(TF_Combined(subI).powspctrm,[1 3 2]);
        pow_de = pow - mean_power;
        TF_Combined(subI).powspctrm = permute(pow_de,[1 3 2]);
    elseif stp.logcorr == 1
        TF_Combined(subI).powspctrm = log(TF_Combined(subI).powspctrm);
    end
    
    
    
    fn_ent = cleanbehav(:,1); %entropy
    fn_surp = cleanbehav(:,2); %surprise
    fn_meanent = cleanbehav(:,3); %covar - mean entropy
    for trlI = 1:trls
        fn_meanpow(trlI,1) = mean2(TF_Combined(subI).powspctrm(trlI,:,:));
    end
    fn_trl = cleanbehav(:,5);
    
    % fit glm with predictors
    for freqI = 1:size(TF_Combined(subI).powspctrm,2)
        for timeI = 1:size(TF_Combined(subI).powspctrm,3)
            if stp.trlN == 0
                xs = [fn_ent,fn_surp,fn_meanent,fn_meanpow];
            else
                xs = [fn_ent,fn_surp,fn_meanent,fn_meanpow,fn_trl];
            end
            % run glm
            b(1:size(xs,2)+1,freqI,timeI) = glmfit(xs,TF_Combined(subI).powspctrm(:,freqI,timeI));
        end
    end
    
    const = b(1,:,:);
    entropy = b(2,:,:);
    surprise = b(3,:,:);
    meanent = b(4,:,:);
    
    TF_ent{subI}.powspctrm = entropy;
    TF_ent{subI}.time = TF_Combined(subI).time;
    TF_ent{subI}.freq = TF_Combined(subI).freq;
    TF_ent{subI}.label = {stp.region};
    TF_ent{subI}.dimord = 'chan_freq_time';
    
    TF_surp{subI}.powspctrm = surprise;
    TF_surp{subI}.time = TF_Combined(subI).time;
    TF_surp{subI}.freq = TF_Combined(subI).freq;
    TF_surp{subI}.label = {stp.region};
    TF_surp{subI}.dimord = 'chan_freq_time';
    
    TF_meanent{subI}.powspctrm = meanent;
    TF_meanent{subI}.time = TF_Combined(subI).time;
    TF_meanent{subI}.freq = TF_Combined(subI).freq;
    TF_meanent{subI}.label = {stp.region};
    TF_meanent{subI}.dimord = 'chan_freq_time';
    
    TF_const{subI}.powspctrm = const;
    TF_const{subI}.time = TF_Combined(subI).time;
    TF_const{subI}.freq = TF_Combined(subI).freq;
    TF_const{subI}.label = {stp.region};
    TF_const{subI}.dimord = 'chan_freq_time';
    
    null{subI}.powspctrm  = zeros(1,freqs,timepoints);
    null{subI}.time = TF_ent{subI}.time;
    null{subI}.freq = TF_ent{subI}.freq;
    null{subI}.dimord = 'chan_freq_time';
    null{subI}.label = {stp.region};
    
    k=k+1;
    clearvars -except analysisDir yh yl x clims nsubs aHPC Hicks k demean stp lat freqrange subI null run_glm ...
        logcorr freqrange Lnull Hnull TF_surp_Low TF_surp TF_ent TF_meanent TF_const TF_int clean_T_TF TF_Combined comb
    
end


%% grand average
fprintf('Getting group-averages ready \n')
cfg = [];
cfg.keepindividual = 'yes';
TF_ent_grandavg = ft_freqgrandaverage(cfg, TF_ent{:});
TF_surp_grandavg = ft_freqgrandaverage(cfg, TF_surp{:});
TF_meanent_grandavg = ft_freqgrandaverage(cfg, TF_meanent{:});
TF_null_grandavg = ft_freqgrandaverage(cfg, null{:});
varargout = {TF_ent,TF_surp,TF_meanent,TF_const,TF_ent_grandavg,TF_surp_grandavg,TF_meanent_grandavg,TF_null_grandavg};

end

