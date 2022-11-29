function [behav_time] = behav_time_clean_ctx(nsubs,patient_data,stp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
tr = 40;
ep = 12;
for subI = 1:numel(nsubs)
    fprintf(['Starting behav & time cleaning for Patient ',num2str(nsubs(subI)), '\n'])
    % load behavioral data
    if stp.zurich(subI)
        foldn = sprintf('Information/P%dz',nsubs(subI));
    else
        foldn = sprintf('Information/Patient%d+',nsubs(subI));
    end
    cd(foldn)
    if ~exist(sprintf('%s',stp.region),'dir')
        load('results')
        load('Xfe')
        load('ProcessedBehav')
        
        rt = results(end-479:end,end-6:end);
                
        % these are fail-safes        
        neg_rt = find(rt(:,3)<=0 | rt(:,3)>=1000);
        results = results(end-479:end,end-6:end);
        incTrials = find(results(:,end)==1);
        neg_rt_unique = neg_rt(~ismember(neg_rt,incTrials));
        remove_idx = sortrows([neg_rt_unique; incTrials]);
        if (nsubs(subI) == 6 && stp.zurich(subI)) || nsubs(subI) ==27
            corrTrial = data(:,6);
            rtc = data(:,5);
        else
            corrTrial = setdiff(1:480,remove_idx)';
                    rtc = rt(corrTrial,3);

        end
        
        Xfe = Xfe(end-479:end,end-5:end);
        Xfec = Xfe(corrTrial,:);
        
        % calculate entropy and surprise
        karl = [1 1 -1 0 0 0;0 0 0 1 1 -1];
        X = (Xfec*karl');
        % get mean entropy (already corrected for trials);
        me = data(:,3);
        behav = [X, me,rtc, corrTrial];
        
        % load ieeg file
        if stp.zurich(subI)
            % load neuralynx files
            fn                       = sprintf('/Volumes/Promise Pegasus/Darya/WP4/Information/P%dz/m',nsubs(subI));
            evR                      = ft_read_event(fn);
            cfg                      = [];
            cfg.dataset              = fn;
            cfg.trialdef.eventtype   = 'trigger';
            cfg.trialdef.eventvalue  = 10;
            cfg.trialdef.prestim     = stp.clean_lat; %prestim;
            cfg.trialdef.poststim    = stp.clean_lat; %poststim;
            cfg.trialfun             = 'ft_trialfun_general';
            cfg                      = ft_definetrial(cfg);
        else
            % read .edf
            fn               = sprintf('s%d_Info.edf', nsubs(subI));
            cfg              = [];
            cfg.dataset      = fn;
            cfg.prestim      = stp.clean_lat; %prestim;
            cfg.poststim     = stp.clean_lat; %poststim;
            cfg.triggers     = zeros(1,480);
            cfg.nsub         = nsubs(subI);
            
            if nsubs(subI) == 25 || nsubs(subI) == 36 || nsubs(subI) == 37
                cfg.trialfun = 'mytrialfun_IntraCranial';
            else
                cfg.trialfun = 'mytrialfun_IntraCranial2';
            end
            cfg              = ft_definetrial(cfg);
        end
        
        trl = cfg.trl;
        if size(trl,1) ~= 480
            if nsubs(subI) == 25
                cfg.trl = cfg.trl(34:513,:); % first 33 events are from previous failed launch (log ending in 09-08)
            elseif nsubs(subI) == 36
                cfg.trl = cfg.trl(40:519,:); % first 39 events are from previous failed launch (log ending in 09-08)
                        elseif nsubs(subI) == 37
                cfg.trl = [zeros(1,3); cfg.trl]; % missing first trigger
                trl = cfg.trl;
            else
                error('trial number isn''t equal to 480')
            end
        end
        
        % keep only correct trials from behaviour
        cfg.trl = trl(corrTrial,:);
        
        %  Preprocessing
        cfg.detrend          = 'no';  % following dlz code for zurich
        cfg.demean           = 'yes'; % following dlz code for zurich
        cfg.continuous       = 'yes';
        cfg.montage.tra      = patient_data(subI).bipolmat;
        cfg.montage.labelnew = patient_data(subI).bipolabel;
        cfg.montage.labelorg = patient_data(subI).labelorg;
        cfg.channel          = patient_data(subI).chanofinterest;
        trialdata = ft_preprocessing(cfg);
        
        if stp.resp_lock == 1
            rtsec = behav(:,end-1)/1000;
            offset = round(rtsec * trialdata.fsample);
            %  shift the time axis
            cfg        = [];
            cfg.offset = -offset;
            trialdata  = ft_redefinetrial(cfg, trialdata);
        end
        
        % artefact rejection
        cfg                             = [];
        cfg.channel                     = 'all';
        cfg.viewmode                    = 'vertical';
        cfg.artifactalpha               = 0.8;
        cfg.artfctdef.contact1.artifact = [];
        cfg.artfctdef.contact2.artifact = [];
        cfg.artfctdef.contact3.artifact = [];
        artf                            = ft_databrowser(cfg,trialdata);
        
        % now reject the marked trials
        cfg                             = [];
        cfg.artfctdef.reject            = 'nan';
        cfg.artfctdef.miaccepttime      = 0;
        cfg.artfctdef.value             = 'nan';
        cfg.artfctdef.visual.artifact   = artf.artfctdef.visual.artifact;
        cfg.artfctdef.conract1.artifact   = artf.artfctdef.contact1.artifact;
        cfg.artfctdef.conract2.artifact   = artf.artfctdef.contact2.artifact;
        cfg.artfctdef.conract3.artifact   = artf.artfctdef.contact3.artifact;
        cleandata                       = ft_rejectartifact(cfg,trialdata);
        
        % remove noisy trials from behaviour too
        count = 1;
        countk = 1;
        for trialI = 1:size(cleandata.trial,2)
            if any(isnan(cleandata.trial{1,trialI}(1,:)))
                rmvtrials(count) = trialI;
                count = count +1;
            else
                keeptrials(countk) = trialI;
                countk=countk+1;
            end
        end
        
        % if no trials removed
        if count == 1
            rmvtrials = [];
        end
        
        cfg        = [];
        cfg.trials = keeptrials;
        cleandata  = ft_selectdata(cfg, cleandata);
        
        if stp.resp_lock == 1
            artf_bt = sprintf('artf_behav_time_resp_lock_%s.mat',stp.montageD);
        else
            artf_bt = sprintf('artf_behav_time_%s.mat',stp.montageD);
            
        end
        
        behav_time(subI).nsub = nsubs(subI);
        behav_time(subI).rmvtrials = rmvtrials;
        behav_time(subI).keeptrials = keeptrials;
        behav_time(subI).info = behav(keeptrials,:);
        
        art = behav_time(subI);
        mkdir(sprintf('%s',stp.region))
        cd(sprintf('%s',stp.region))
        save(artf_bt,'art');
        cd ../../../
        clearvars -except nsubs stp patient_data ep tr behav_time
        
    else %if folder already exists ~ already cleaned
        cd ../../
        clearvars -except nsubs stp patient_data ep tr behav_time
        
        continue
    end
    
    
end
end