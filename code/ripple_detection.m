% Vaz 2019,2020 Ripple detection method
% Assumes loaded nsubs stp fron pipeline
nsubs = [3,6,8,9,13,15,16,22,25,31,32,36,37,6,8,10]; %12 and 22 not used poor data
iszurich = logical([zeros(1,13),ones(1,3)]);
[patient_data,stp] = setup(nsubs,iszurich,0,0,'anterior');

clear patient_data
patient_data = getMontage(nsubs,stp,0);
hpfilt = 200;
ripdur = 25;
fname = sprintf('HPCRipples/HPCAnterior_vaz_hpf%d_%dms_%dsubjs.mat',hpfilt, ripdur, numel(nsubs));
%
for subI = 1:numel(nsubs)
    fprintf(['Getting clean trials for Patient ',num2str(nsubs(subI)), '\n'])
    if stp.zurich(subI) == 1
        foldn = sprintf('Information/P%dz',nsubs(subI));
    else
        foldn = sprintf('Information/Patient%d+',nsubs(subI));
    end
    cd(foldn)
    if strcmp(stp.region,'anterior') || strcmp(stp.region,'head')||strcmp(stp.region,'body')
        cd(sprintf('hpc_%s',patient_data(subI).hpc_axis))
    else
        cd(sprintf('%s',patient_data(subI).region))
    end
    load('clean_trials_bipolar.mat','clean')
    clean_trials(subI) = clean;
    cd ../../../
    
    stp.poststim = 1;
    stp.prestim = 1;
    timest={-stp.prestim:0.002:stp.poststim};
    
    cleandata = preproc(subI,nsubs,stp,patient_data,clean.trl);
    timest = cleandata.time(1,1);
    if stp.zurich(subI) == 1 || nsubs(subI) == 31|| nsubs(subI) == 32 ...
            || nsubs(subI) == 36 || nsubs(subI) == 37
%         downsample
        cfg = [];
        cfg.time = repmat({timest{1}(1:end-1)},1,size(cleandata.trial,2));
        cfg.demean = 'yes';
        cleandata{subI} = ft_resampledata(cfg,cleandata);
    end

    cfg=[];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [80 120];
    cfg.bpfilttype = 'but';
    cfg.bpfiltdir = 'twopass';
    cfg.bpfiltord = 2;
    cfg.demean = 'yes';
    signal = ft_preprocessing(cfg,cleandata);
    get Hilbetrt for ripple detection
    cfg.hilbert = 'abs';
    rip80_120 = ft_preprocessing(cfg,cleandata);
    
%     higher freqs for IED
    cfg=[];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = hpfilt;
    cfg.hpfilttype = 'but';
    cfg.hpfiltdir = 'twopass';
    cfg.hpfiltord = 2;
    cfg.demean = 'yes';
    ied_signal = ft_preprocessing(cfg,cleandata);
    
    cfg.hilbert = 'abs';
    ied = ft_preprocessing(cfg,cleandata);
    
    
    for elecI = 1:size(rip80_120.label,2)
        hilbmap_rip = [];
        hilbmap_ied = [];
        eeg = [];
%         get data in trials x durations structure
        for trialI = 1:size(rip80_120.trial,2)
            hilbmap_rip(trialI,:) = rip80_120.trial{1,trialI}(elecI,:);
            hilbmap_ied(trialI,:) = ied.trial{1,trialI}(elecI,:);
            eeg(trialI,:) = cleandata.trial{1,trialI}(elecI,:);
        end
        [sub(subI).ripplelogic(elecI,:,:),sub(subI).iedlogic(elecI,:,:),...
            sub(subI).hilb_sig(elecI,:,:),...
            sub(subI).hilb_ied(elecI,:,:)] = detectripples_iedreject_df(hilbmap_rip,...
            hilbmap_ied, eeg, ripdur, 2, 3, cleandata.fsample);
        
    end
%     sort everything such that unqiue_rips will be:
%     [channel, trial, start_ripple, end_ripple]
    [xtime,y,z]=ind2sub(size(sub(subI).ripplelogic),find(sub(subI).ripplelogic));
    [xi,yi,zi] = ind2sub(size(sub(subI).iedlogic),find(sub(subI).iedlogic));
    sub(subI).rips = [xtime,y,z];
    sub(subI).rips = sortrows(sub(subI).rips);
    sub(subI).ieds = [xi,yi,zi];
    sub(subI).ieds = sortrows(sub(subI).ieds);
    
%     find non-consecutive timepoints:
    unq = find(diff(sub(subI).rips(:,3)) ~= 1) + 1;
    if isempty(unq)
        continue
    else
        sub(subI).unique_rips = [sub(subI).rips(1,:);sub(subI).rips(unq,:)];
        sub(subI).unique_rips(:,4) = [sub(subI).rips((unq)-1,3);sub(subI).rips(end,3)];
    end
    unqi = find(diff(sub(subI).ieds(:,3)) ~= 1) + 1;
    if ~isempty(unqi)
        sub(subI).unique_ieds = [sub(subI).ieds(1,:);sub(subI).ieds(unqi,:)];
        sub(subI).unique_ieds(:,4) = [sub(subI).ieds((unqi)-1,3);sub(subI).ieds(end,3)];
    else
        sub(subI).unique_ieds = [];
    end
    xtime = timest{1,1}(1:end-1);
    
%     plot trials with ripples
    for ripI = 1:size(sub(subI).unique_rips,1)
%         get peak of ripples
        t_idx = sub(subI).unique_rips(ripI,2);
        t = clean_trials(subI).trl(t_idx); % trial # insead of index
        ent(subI,ripI) = clean_trials(subI).info(t_idx,1);
        surp(subI,ripI) = clean_trials(subI).info(t_idx,2);
        start_idx=sub(subI).unique_rips(ripI,3);
        end_idx=sub(subI).unique_rips(ripI,4);
        start(subI,ripI)=xtime(start_idx);
        peak = round((start_idx+end_idx)/2);
        el = sub(subI).unique_rips(ripI,1);
        
        ripple_peak = cleandata.trial{1,t_idx}(el,[start_idx:end_idx]);
        
        [a,b]=findpeaks(ripple_peak);
        if ~isempty(b)
            if (mod(length(a),2)==0) %If the number of the peaks is even
                [~,d]=max(a(length(a)/2:length(a)/2+1));
                peak=start_idx+b(length(a)/2+d-1)-1;
            else %If the number of the peaks is odd
                peak=start_idx+b(ceil(length(a)/2))-1;
            end
        end
        sub(subI).unique_rips(ripI,5) = peak;
        
%         for multiple ripple per trial
        if ripI < size(sub(subI).unique_rips,1)-1
            t_idx = sub(subI).unique_rips(ripI,2);
        else
            continue
        end
        
                t = clean_trials(subI).trl(t_idx); % trial # insead of index
                if t_idx == sub(subI).unique_rips(ripI+1,2)
                    flag = 'mult';
                else
                    flag = '';
                end
                el = sub(subI).unique_rips(ripI,1);
                xbars = [xtime(sub(subI).unique_rips(ripI,3)) xtime(sub(subI).unique_rips(ripI,4))]; %for shading
                figure('units','normalized','outerposition',[0 0 1 1]);
                sgtitle(sprintf('Patient %d, trial %d channel %s',nsubs(subI),t,cleandata.label{el}))
                subplot(311);
                [coef,f]=cwt(cleandata.trial{1,t_idx}(el,:),'amor',cleandata.fsample);
                imagesc(xtime,f,abs(coef));
                set(gca,'YDir','normal')
                hold on
                rm = ones(1,size(xtime,2))*130;
                plot(peak,rm,'vr','MarkerFaceColor','r');
                ylabel('Frequency (Hz)');
                hold off
        
                subplot(312);
                y = signal.trial{sub(subI).unique_rips(ripI,2)}(sub(subI).unique_rips(ripI,1),:);
                hp = plot(xtime, y);
                hold on
                patch([xbars(1) xbars(1), xbars(2) xbars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.9 0.9 0.9])
                plot(xtime, y)
                ylabel({"Ripple Band";80+"-"+120+" Hz";"(uV)"});
                hold off
        
                subplot(313)
                hilb_sig = squeeze(sub(subI).hilb_sig(el,t_idx,:));
                hp2 = plot(xtime, hilb_sig);
                hold on
                patch([xbars(1) xbars(1), xbars(2) xbars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.9 0.9 0.9])
                plot(xtime, hilb_sig)
                ripmin = ones(1,length(xtime))*2;
                ripmaxthresh = ones(1,length(xtime))*3;
                plot(xtime,ripmin,'--r');
                plot(xtime,ripmaxthresh,'--b')
                ylabel('Standardised Hilbert Amplitude')
                xlabel('Time (s)')
                set(findall(gcf,'-property','FontSize'),'FontSize',18)
        pause(2)
        
                cd(sprintf('HPCRipples/Vaz_%dms_hpfilt%d_perContact',ripdur,hpfilt))
                if stp.zurich(subI)
                    fn = sprintf('Rs%dzurich_%s_%s_t%d%s',nsubs(subI),stp.region,cleandata.label{el},t,flag);
                else
                    fn = sprintf('Rs%d_%s_%s_t%d%s',nsubs(subI),stp.region,cleandata.label{el},t,flag);
                end
                saveas(gcf,fn,'png');
                close all
                cd ../../
    end
    
%        plot IEDs
        for iedI = 1:size(sub(subI).unique_ieds,1)
            % for multiple ieds per trial
            if iedI < size(sub(subI).unique_ieds,1)-1
                t_idx = sub(subI).unique_ieds(iedI,2);
            else
                continue
            end
    
            t = clean_trials(subI).trl(t_idx); % trial # insead of index
            if t_idx == sub(subI).unique_ieds(iedI+1,2)
                flag = 'mult';
            else
                flag = '';
            end
            el = sub(subI).unique_ieds(iedI,1);
            xbars = [xtime(sub(subI).unique_ieds(iedI,3)) xtime(sub(subI).unique_ieds(iedI,4))]; %for shading
            midpoint = (xbars(2)+xbars(1))/2;
    
            figure('units','normalized','outerposition',[0 0 1 1]);
            sgtitle(sprintf('Patient %d, trial %d channel %s',nsubs(subI),t,cleandata.label{el}))
            subplot(3,2,[1,2]);
            [coef,f]=cwt(cleandata.trial{1,t_idx}(el,:),'amor',cleandata.fsample);
            imagesc(xtime,f,abs(coef));
            set(gca,'YDir','normal')
            hold on
            rm = ones(1,size(xtime,2))*200;
            plot(midpoint,rm,'vr','MarkerFaceColor','r');
            ylabel('Frequency (Hz)');
            hold off
    
            subplot(323);
            y = ied_signal.trial{sub(subI).unique_ieds(iedI,2)}(sub(subI).unique_ieds(iedI,1),:);
            hp = plot(xtime, y);
            hold on
            patch([xbars(1) xbars(1), xbars(2) xbars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.9 0.9 0.9])
            plot(xtime, y)
            ylabel({"IED Band";hpfilt+"-"+250+" Hz";"(uV)"});
            hold off
    
            subplot(325)
            hilb_sig = squeeze(sub(subI).hilb_ied(el,t_idx,:));
            hp2 = plot(xtime, hilb_sig);
            hold on
            patch([xbars(1) xbars(1), xbars(2) xbars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.9 0.9 0.9])
            plot(xtime, hilb_sig)
            ylabel('Standardised Hilbert Amp (IED)')
            xlabel('Time (s)')
    
            subplot(324)
            y = signal.trial{sub(subI).unique_ieds(iedI,2)}(sub(subI).unique_ieds(iedI,1),:);
            hp = plot(xtime, y);
            hold on
            patch([xbars(1) xbars(1), xbars(2) xbars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.9 0.9 0.9])
            plot(xtime, y)
            ylabel({"Ripple Band";80+"-"+120+" Hz";"(uV)"});
            hold off
    
            subplot(326)
            hilb_sig = squeeze(sub(subI).hilb_sig(el,t_idx,:));
            hp2 = plot(xtime, hilb_sig);
            hold on
            patch([xbars(1) xbars(1), xbars(2) xbars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.9 0.9 0.9])
            plot(xtime, hilb_sig)
            ripmin = ones(1,length(xtime))*2;
            ripmaxthresh = ones(1,length(xtime))*3;
            plot(xtime,ripmin,'--r');
            plot(xtime,ripmaxthresh,'--b')
            ylabel('Standardised Hilbert Amp (Ripple)')
            xlabel('Time (s)')
    
            set(findall(gcf,'-property','FontSize'),'FontSize',18)
    pause(2)
    
            cd(sprintf('HPCRipples/Vaz_%dms_hpfilt%d_perContact',ripdur,hpfilt))
            if stp.zurich(subI)
                fn = sprintf('Is%dzurich_%s_%s_t%d%s',nsubs(subI),stp.region,cleandata.label{el},t,flag);
            else
                fn = sprintf('Is%d_%s_%s_t%d%s',nsubs(subI),stp.region,cleandata.label{el},t,flag);
            end
            saveas(gcf,fn,'png');
            close all
            cd ../../
    
        end
    
    clearvars -except fname nsubs stp patient_data sub hpfilt ripdur clean_trials
    
end
save(fname, 'sub')
% average ripples
fname = sprintf('HPCRipples/HPCAnterior_vaz_hpf%d_%dms_medial_%dsubjs.mat',hpfilt, ripdur, numel(nsubs));
load(fname)

patient_data = getMontage(nsubs,stp,0);
ent = zeros(1,1);
start = zeros(1,1);
times = [-1:0.002:1]';
surp = zeros(1,1);
g=1;
j=1;
m=1;
figure(1);figure(2)
for subI = 1:numel(nsubs)
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
    
    stp.poststim = 1;
    stp.prestim = 1;
    timest={-stp.prestim:0.002:stp.poststim};
    cleandata = preproc(subI,nsubs,stp,patient_data,clean.trl);
    if stp.zurich(subI) == 1 || nsubs(subI) == 31 || nsubs(subI) == 32 || ...
            nsubs(subI) == 36 || nsubs(subI) == 37
%         downsample - for zurich & ruber
        cfg = [];
        cfg.time = repmat({timest{1}(1:end-1)},1,size(cleandata.trial,2));
        cleandata = ft_resampledata(cfg,cleandata);
    end
    if nsubs(subI) == 22 % no ripples for sub 22 (and excluded anyways)
        continue
    end
    ent_tot(g:g+size(clean_trials(subI).info(:,1),1)-1) = clean_trials(subI).info(:,1);
    
    tfr = [];
    for ripI = 1:size(sub(subI).unique_rips,1)
        t_idx = sub(subI).unique_rips(ripI,2);
        t = clean_trials(subI).trl(t_idx); % trial # insead of index
        
        ent(subI,ripI) = clean_trials(subI).info(t_idx,1);
        surp(subI,ripI) = clean_trials(subI).info(t_idx,2);
        rt(subI,ripI) = clean_trials(subI).info(t_idx,4);
        meanent(subI,ripI) = clean_trials(subI).info(t_idx,3);
        start_idx=sub(subI).unique_rips(ripI,3);
        end_idx=sub(subI).unique_rips(ripI,4);
        start(subI,ripI)=times(start_idx);
        midpoint = round((start_idx+end_idx)/2);
        ripdur = times(end_idx+1)-times(start_idx-1);
        
        elecI = sub(subI).unique_rips(ripI,1);
        
        ripple_peak = cleandata.trial{1,t_idx}(elecI,[start_idx:end_idx]);
        
        [~,peak_idx] = findpeaks(ripple_peak);
        if ~isempty(peak_idx)
            midpoint =start_idx+ peak_idx(end);
        end
        
        from marcos v6 staresina
        [a,b]=findpeaks(ripple_peak);
        if ~isempty(b)
            if (mod(length(a),2)==0) %If the number of the peaks is even
                [~,d]=max(a(length(a)/2:length(a)/2+1));
                midpoint=start_idx+b(length(a)/2+d-1)-1;
            else %If the number of the peaks is odd
                midpoint=start_idx+b(ceil(length(a)/2))-1;
            end
        end
        
        sub(subI).unique_rips(ripI,5) = midpoint;
        sub(subI).unique_rips(ripI,6) = ripdur;
        sub(subI).unique_rips(ripI,7) = ent(subI, ripI);
        sub(subI).unique_rips(ripI,8) = surp(subI, ripI);
        sub(subI).unique_rips(ripI,9) = meanent(subI, ripI);
        sub(subI).unique_rips(ripI,10) = rt(subI, ripI);
        sub(subI).unique_rips(ripI,11) = t; %trialnum not idx
        
        time_idx = [];
        if midpoint < 201 % ripple in the very beginnig
            nan_pad = 201-midpoint;
            time_idx = [1:midpoint+200];
            rip_signal(j,:) = [nan(1,nan_pad),cleandata.trial{1,t_idx}(elecI,time_idx)];
            subj_ripl(k,:) = [nan(1,nan_pad),cleandata.trial{1,t_idx}(elecI,time_idx)];
            tfr.trial{k} = [zeros(1,nan_pad),cleandata.trial{1,t_idx}(elecI,time_idx)];
            tfr.time{k} = linspace(-400,400,401)/1000;
            
            j=j+1; k =k+1;
        elseif midpoint > 800 % ripple in the very end
            time_idx = [midpoint-200:1000];
            nan_pad = 401-(length(time_idx)); %there are 201 timepoints - 400ms
            rip_signal(j,:) = [cleandata.trial{1,t_idx}(elecI,time_idx),nan(1,nan_pad)];
            subj_ripl(k,:) = [cleandata.trial{1,t_idx}(elecI,time_idx),nan(1,nan_pad)];
            tfr.trial{k} = [cleandata.trial{1,t_idx}(elecI,time_idx),zeros(1,nan_pad)];
            tfr.time{k} = linspace(-400,400,401)/1000;
            
            j=j+1; k=k+1;
        else
            time_idx = [midpoint-200:midpoint+200];
            rip_signal(j,:) = cleandata.trial{1,t_idx}(elecI,time_idx);
            subj_ripl(k,:) = cleandata.trial{1,t_idx}(elecI,time_idx);
            tfr.trial{k} = cleandata.trial{1,t_idx}(elecI,time_idx);
            tfr.time{k} = linspace(-400,400,401)/1000;
            j=j+1; k=k+1;
            
        end
        sanity check get power per trial avg later
        [coef_trl{k-1},~]=cwt(tfr.trial{k-1},'amor',500);
        
    end
    g=g+size(clean_trials(subI).info(:,1));
    
    figure(1)
    subplot(4,4,m);
    plot(linspace(-400,400,401),nanmean(subj_ripl,1))
    title({sprintf('mean ripples patient %d',nsubs(subI)), sprintf('N = %d',size(subj_ripl,1))})
    xtick = [-200 -150 -100 -50 0 50 100 150 200];
    set(gca,'xtick',xtick,'TickDir','out');
    ylabel('Voltage (uV)');
    xlabel('Time from ripple midpoint (ms)');
    
    figure(3)
    subplot(4,4,m);
    [coefs(:,:,subI),f]=cwt(nanmean(subj_ripl,1),'amor',500);
    surface (linspace(-400,400,401),f,abs(coefs(:,:,subI)),'EdgeColor','none');
    axis xy
    title({sprintf('TF patient %d',nsubs(subI)), sprintf('N = %d',size(subj_ripl,1))})
    ylabel('Frequency (Hz)');
    xlabel('Time from ripple peak (ms)');
    colorbar();
    
    sanity check of trial-by-trial power then avg
    for i = 1:size(coef_trl,2); coef_trl_mtrx = abs(coef_trl{i}); coef_3d(:,:,i) = coef_trl_mtrx; end
    figure(123);subplot(4,4,m);
    surface(linspace(-400,400,401),f,mean(coef_3d,3),'EdgeColor','none');axis xy; colorbar();
    subj_ripl_avg(subI,:) = nanmean(subj_ripl);
    m=m+1;
    clear subj_ripl  coef_3d
    
end


subj_ripl_avg(8,:) = [];
figure('position',[10 10 1100 700]);
subplot(121)
plot(linspace(-400,400,401),nanmean(subj_ripl_avg,1),'LineWidth',2)
title({'mean ripples across patients', sprintf('N = %d',size(rip_signal,1))})
ylabel('Voltage (uV)');
xlabel('Time relative to ripple peak (ms)');
set(gca,'FontSize',20);
subplot(122)
[coef,f]=cwt(nanmean(subj_ripl_avg,1),'amor',500); %exact same as abs(mean(coefs,3)) from all patients, just faster
surface (linspace(-400,400,401),f,abs(coef),'EdgeColor','none');
set(gca,'YDir','normal','FontSize',20);
ylabel('Frequency (Hz)');
xlabel('Time relative to ripple peak (ms)');
ylim([3 185])
colorbar();
print('-dtiff',fullfile('Manuscript/Figures/','avg_rip_pow_vol_row'),['-r' '300'])


save(fname, 'sub')
