%% STEP 1 of analysis - Recereate behavioural data from log file
%
% This creates the same values as y.mat (negative and correct), + positive
% values over the time limit the are represented as 0s in y.mat
% Darya Frank 07-09-2020

clearvars
nsubs = [3,6,8,9,12,13,15,25,27,16,31,32,36,37,6,8,10];
iszurich = logical([zeros(1,14),ones(1,3)]);
analysisDir = pwd;

for subI =1:numel(nsubs)
    % load behavioral data
    if iszurich(subI) == 1
        foldn = sprintf('Information/P%dz',nsubs(subI));
    else
        foldn = sprintf('Information/Patient%d+',nsubs(subI));
    end
cd(foldn)
load('results.mat')
d = dir('*.log');
log_file = d(end).name;
s=importdata(log_file);
vld = find(~cellfun(@isempty, s.textdata(:,2)));
cog_start = find(contains(s.textdata(:,1),'COGENT START'));
task_start = cog_start+1;
txt=s.textdata(task_start:vld(end),[2,4,5,6]);
t1=str2double(s.textdata(task_start:vld(end),1));
first_press = find(strcmp(txt(:,4),'DOWN'));
beg_trail = vld(1)-cog_start-1;
first_trls = s.textdata(cog_start+1:cog_start+beg_trail,1);

nums = regexp(first_trls,'[\d*\.]*\d*','match');
trls = size(nums,1);
add = [];
for trlI = 1: trls
 add = [add; str2num(nums{trlI,1}{1,1})];
end

t1(1:beg_trail,1) = add;
txt(1:beg_trail,2) = repmat({'colX'},beg_trail,1);

end_trail = length(txt)-length(s.data)-beg_trail;
time = num2cell([zeros(first_press(1)-1,1);s.data;zeros(end_trail,1)]);

log_tmp = [num2cell(t1),txt,time];
if strcmp(log_tmp(end,3),'COGENT STOP')
log_tmp(end,:) = [];
end

down_ind=find(strcmp(log_tmp(:,5),'DOWN')); % index for key press
up_ind=find(strcmp(log_tmp(:,5),'UP')); % index for key press

pic_ind=find(~contains(log_tmp(1:end-1,3),'Key')); %index for the corresponding stimulus
% pic number as defined by coegent
pic_num = results(:,end-2);

% get the corresponding keys (2,13,14,22) for each pic num
if subI == 15 && nsubs(subI) == 6 %P6z
    key1_shape = [74,124,48,50,14,12,66,10,20,100];
    key2_shape = [128,96,54,62,78,122,32,108,126,60,102,106];
    key3_shape = [80,88,22,46,98,94,58];
    key4_shape = [72,76,40,18,6,4,2,42,36];
else
key1_corr = find(results(:,end-3)==2 & results(:,end) == 0);
key1_shape = unique(results(key1_corr,end-2));

key2_corr = find(results(:,end-3)==13 & results(:,end) == 0);
key2_shape = unique(results(key2_corr,end-2));

key3_corr = find(results(:,end-3)==14 & results(:,end) == 0);
key3_shape = unique(results(key3_corr,end-2));

key4_corr = find(results(:,end-3)==22 & results(:,end) == 0);
key4_shape = unique(results(key4_corr,end-2));
end
% fill in the new log
log_rec = zeros(480,6);
trl = 1;
for stimI = 1:numel(pic_ind)
    % find closest down press
    stim_idx = pic_ind(stimI);
    
    if ismember(stim_idx+1,down_ind) %pressed down immediately after pic (_ pic - down - up_)
        idx_down = find(stim_idx+1==down_ind);
        down_press = down_ind(idx_down);
        RT = cell2mat(log_tmp(down_press,6))-cell2mat(log_tmp(stim_idx,1));
        keyp = str2double(cell2mat(log_tmp(down_press,4)));
        
    elseif ismember(stim_idx+2,down_ind) %had an upkey pressed before down (down - pic - up - down - up)
        idx_down = find(stim_idx+2==down_ind);
        down_press = down_ind(idx_down);
        RT = cell2mat(log_tmp(down_press,6))-cell2mat(log_tmp(stim_idx,1));
        keyp = str2double(cell2mat(log_tmp(down_press,4)));

    else
        RT = 0;
        keyp = 0;
    end

    log_rec(trl,1:4) = [trl,RT,pic_num(stimI),keyp];
    trl=trl+1;
    
end

% convert key presses based on log and results.mat
key1_pres = ismember(log_rec(:,3),key1_shape);
key2_pres = ismember(log_rec(:,3),key2_shape);
key3_pres = ismember(log_rec(:,3),key3_shape);
key4_pres = ismember(log_rec(:,3),key4_shape);

log_rec(key1_pres,5) = 2;
log_rec(key2_pres,5) = 13;
log_rec(key3_pres,5) = 14;
log_rec(key4_pres,5) = 22;

corr_key = find(log_rec(:,4) == log_rec(:,5));
log_rec(corr_key,6) = 1;

% only keep correct trials
log_corr = log_rec(corr_key,[1,2,4,6]);

% remover RT above 2200ms (ISI)
log_corr(log_corr(:,2)>2200,:) = [];

save('recreated_log.mat','log_rec','log_corr')
clearvars -except nsubs analysisDir iszurich
cd(analysisDir)
end
