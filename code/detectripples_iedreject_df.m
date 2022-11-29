function [ripplelogic, iedlogic, zhilb, zied] = detectripples_iedreject_df(hilbamp_rip, hilbamp_ied, eeg, ripplewidth, ripthresh, ripmaxthresh,fsample)
%%
% ripple detection code from Vaz et al, Science, 2019
% written in MATLAB_R2019a
% please direct any questions to vazio92@gmail.com

%%
% INPUT:
% hilbamp_rip  = trials x duration (ms) 
%    - hilbert amplitude data in 80-120Hz frequency band

% hilbamp_ied  = trials x duration (ms)
%    - hilbert amplitude data in 250-500Hz frequency band

% eeg          = trials x duration (ms)
% ripplewidth  = 25ms %df added a line to asjust this for diff sampling
% rates
% ripthresh    = 2; (threshold detection)
% ripmaxthresh = 3; (ripple event must reach this maximum)
% fsample      = sampling rate of data to adjust for ms controls

% OUTPUT:
% ripplelogic  = trials x dur (same size mask for hilbamp_rip)
% iedlogic     = trials x dur (same size mask for hilbamp_ied)

%%
% number of trials
numtrials   = size(eeg,1);

% adjust ripple width for sampling rate DF
ripplewidth = round(fsample/1000*ripplewidth);

% z-score hilbert amplitudes for ripple band and ied band
hilbamp_rip = (hilbamp_rip - mean(hilbamp_rip(:))) / std(hilbamp_rip(:));
hilbamp_ied = (hilbamp_ied - mean(hilbamp_ied(:))) / std(hilbamp_ied(:));

zhilb = hilbamp_rip; %DF
zied = hilbamp_ied; %DF

% create logicals for ied and potential ripples
ripplelogic = hilbamp_rip>ripthresh;
broadlogic  = hilbamp_ied>5; % (Staresina, NatNeuro, 2015)

% measure eeg gradient for IED rejection (Staresina, NatNeuro, 2015)
eegdiff = diff(eeg,1,2); % first derivative over second dimension

% make logical arrays same size
eegdiff(:,end+1) = eegdiff(:,end);

% z-score gradient for artifact rejection
eegdiff = (eegdiff - mean(eegdiff(:))) / std(eegdiff(:));

% create logical for gradient based artifact rejection
difflogic = abs(eegdiff)>5;

% combine the logic from >250Hz (broadband) and gradient
iedlogic = broadlogic | difflogic;

% expand to +/- 100ms on either side of ied index
x = round(fsample/1000*200); % df to adjust for sampling rate
iedlogic = conv2(double(iedlogic), ones(1,x), 'same');
iedlogic = iedlogic > 0;

% remove overlapping indices between ied and ripple logicals
ripplelogic(iedlogic) = 0;

% loop through trials to remove false ripples
for trial = 1:numtrials
    
    ripplelogictrial = ripplelogic(trial,:);
    hilbamptrial     = hilbamp_rip(trial,:);
    
    % continue if no ripples
    if sum(ripplelogictrial) == 0
        continue
    end
    
    % extract ripples
    groups = bwconncomp(ripplelogictrial);
    inds   = groups.PixelIdxList;

    % loop through ripples
    for ripple = 1:length(inds)
        
        % remove ripple from trial if less than 25ms or if max amplitude is less than 3*std dev
        if      numel(inds{ripple}) < ripplewidth
            ripplelogictrial(inds{ripple}) = false;
        elseif  max(abs(hilbamptrial(inds{ripple}))) < ripmaxthresh
            ripplelogictrial(inds{ripple}) = false;
        end    

    end
    
    % reassign trial with removed ripples
    ripplelogic(trial,:) = ripplelogictrial;
    
end

% join ripples that are less than 15ms separated (Roux, NatNeuro, 2017)
for trial = 1:numtrials
    
    ripplelogictrial = ripplelogic(trial,:);

    % continue if no ripples
    if sum(ripplelogictrial) == 0
        continue
    end
    
    % extract ripples
    groups = bwconncomp(ripplelogictrial);
    inds   = groups.PixelIdxList;

    % continue if only 1 ripple
    if length(inds) == 1
        continue
    end
    
    % loop through ripples
    for ripple = 1:length(inds)-1
        if inds{ripple+1}(1) - inds{ripple}(end) < fsample/1000*15 %15ms between ripples
            ripplelogictrial(inds{ripple}(end):inds{ripple+1}(1)) = true;
        end
    end
    
    % reassign trial with joined ripples
    ripplelogic(trial,:) = ripplelogictrial;
    
end

end
