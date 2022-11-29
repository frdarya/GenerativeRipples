function [cleandata] = preproc(subI,nsubs,stp,patient_data,keeptrials)
%run ft_preprocessing for subI, output clean data (using keeptrials)
%   Detailed explanation goes here
if stp.zurich(subI)
    % load neuralynx files
    % can't deal with the OneDrive folder name for some reason
    fn                       = fullfile(sprintf('/Volumes/Promise_Pegasus/Darya/Zurich/P%dz/m',nsubs(subI))); %Users/daryafrank/Zurich/
    evR                      = ft_read_event(fn);
    cfg                      = [];
    cfg.dataset              = fn;
    cfg.trialdef.eventtype   = 'trigger';
    cfg.trialdef.eventvalue  = 10;
    cfg.trialdef.prestim     = stp.prestim; %prestim;
    cfg.trialdef.poststim    = stp.poststim; %poststim;
    cfg.trialfun             = 'ft_trialfun_general';
    cfg                      = ft_definetrial(cfg);
else
    % read .edf
    fn = sprintf('/Volumes/Promise_Pegasus/Darya/OneDrive - Universidad Politécnica de Madrid/WP4/Information/Patient%d+/s%d_Info.edf', nsubs(subI),nsubs(subI));%'/Users/daryafrank/OneDrive - Universidad Politécnica de Madrid/WP4/
    cfg             = [];
    cfg.dataset     = fn;
    cfg.prestim     = stp.prestim; %prestim;
    cfg.poststim    = stp.poststim; %poststim;
    cfg.triggers    = zeros(1,480);
    cfg.nsub        = nsubs(subI);
    if nsubs(subI) == 25 || nsubs(subI) == 36 || nsubs(subI) == 37
        cfg.trialfun = 'mytrialfun_IntraCranial';
    else
        cfg.trialfun = 'mytrialfun_IntraCranial2';
    end
    cfg              = ft_definetrial(cfg);
end

trl = cfg.trl;
if size(trl,1) ~= 480 && nsubs(subI) ~= 37
    if nsubs(subI) == 25
        cfg.trl = cfg.trl(34:513,:); % first 33 events are from previous failed launch (log ending in 09-08)
    elseif nsubs(subI) == 36
        cfg.trl = cfg.trl(40:519,:); % first 39 events are from previous failed launch
    else
        error('trial number isn''t equal to 480')
    end
end

%  Preprocessing
cfg.detrend          = 'no';  
cfg.demean           = 'yes'; 
cfg.continuous       = 'yes';
cfg.montage.tra      = patient_data(subI).bipolmat;
cfg.montage.labelnew = patient_data(subI).bipolabel;
cfg.montage.labelorg = patient_data(subI).labelorg;
cfg.channel          = patient_data(subI).chanofinterest;

trialdata = ft_preprocessing(cfg);

if stp.resp_lock == 1
    
    rtsec = keeptrials(:,end-1)/1000;
    offset = round(rtsec * trialdata.fsample);    
    cfg = [];
    cfg.offset = -offset;
    if nsubs(subI) == 37 && keeptrials(end,end) == 480
        cfg.trials = keeptrials(:,end)-1; %missing first trial
    else
        cfg.trials = keeptrials(:,end);
    end
    cleandata = ft_redefinetrial(cfg, trialdata);
else
    % Remove artifact trials
    cfg        = [];
    if nsubs(subI) == 37 && keeptrials(end,end) == 480
        cfg.trials = keeptrials(:,end)-1; %missing first trial

    else
        cfg.trials = keeptrials(:,end);
    end
    cleandata  = ft_selectdata(cfg, trialdata);
end
end

