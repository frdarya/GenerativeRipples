%% STEP 2 of analysis - Get behavioural data
%
% Use entropy, surprise and mean entropy per block as predictors of RT.
% Adapted from code by A. Duggins (Strange et al., 2005)
% Darya Frank 21/02/2020

clearvars
subs = [3,6,8,9,13,15,25,16,31,32,36,37,6,8,10];
iszurich = logical([zeros(1,12),ones(1,3)]);

d = dir();
d = d(7:end);
analysisDir = pwd;
exclude  = 0;

%subjs with large Xfe = [1 (2), 3 (2), 13 (3),24 (2) ,25 (2)]
for subI = 1:numel(subs)
    if iszurich(subI) ==1
        foldn = sprintf('Information/P%dz',subs(subI));
    else
        foldn = sprintf('Information/Patient%d+',subs(subI));
    end
    cd(foldn)
    subname = subs(subI);
    load Xfe;
    load y;
    load results;
    load recreated_log;
    %load new_y
    tr = 40;
    ep = 12;
    
    % Xfe is [entropycolor entropyshape mutualinfo selfinfocolor selfinfoshape
    % mutualselfinfo];
    % karl is a weights vector for a single subject
    % that combines the colour/shape entropy and mutual info columns
    % to generate an entropy column. similarly with self info columns.
    % This relies on the formula:
    % H(S intersection C) = H(S) + H(C) - I(S,C)
    karl = [1 1 -1 0 0 0;0 0 0 1 1 -1];
    
    if size(Xfe,2) ~= 6
        Xfe = Xfe(end-479:end,end-5:end);
    end
    X = (Xfe*karl');
    
    % get mean entropy across block
    blocks = [];
    for k = 1:ep
        blocks = [blocks eye(tr)];
    end
    avinfo = (blocks*reshape(X,(tr*ep),(2)))/ep;
    % The columns in the left half of reshape(X,(tr*ep),(subjects*2)) are
    % entropy main effects for each subject. The columns in the right half
    % are self information main effects. Multiplying on the left by blocks
    % adds the respective elements from different epochs within subject.
    % Dividing by the number of epochs, we get avinfo, which is
    % tr x (subjects*2) with the means of each of the information theoretic quantities
    % in column
    
    avinfo = repmat(avinfo,ep,1);
    
    if size(results,2) ~= 7
        results = results(:,end-6:end);
    end
    if size(y,1) == 480
        rt = y;
    else
        rt = y(end-479:end);
    end
    
    % columns are entropy, surprise, mean entropy (block), respCorr - where 0 = correct,
    % 1 = incorrect, RT
    data = [X, avinfo(:,1), results(:,7),rt, [1:480]' ];
    data(:,5) = log_rec(:,2);
    keeptrls = find((log_corr(:,2)>0) & log_corr(:,2)<1000);
    keeptrls = log_corr(keeptrls,1);
    data = data(keeptrls,:);

    %fit glm with predictors - E, S, Mean E (covariate), to RT;
    meanrt(subI) = mean(data(:,5));
    b = glmfit(data(:,1:3),data(:,5));
    bitsl(subI,1:3) = [subname,num2cell([b(2),b(3)])];
    blm = fitlm(data(:,1:3),data(:,5)); %same as GLM
    bitsl(subI,1:3) = [subs(subI) num2cell(blm.Coefficients.Estimate([2,3])')];
%     save('ProcessedBehav.mat', 'data');

cd(analysisDir)
    
    clearvars -except subs d corr_surp_ent analysisDir exclude subjects bits_ms bitsl iszurich meanrt
end

%% plots
% group average
idx = find(~cellfun(@isempty,bitsl(:,1)));
bitsl = cell2mat(bitsl(idx,2:3));
err_ent = std(bitsl(:,1))/sqrt(numel(bitsl));
err_surp = std(bitsl(:,2))/sqrt(numel(bitsl));
err = [err_ent,err_surp];
figure();bar(mean(bitsl)); xticklabels({'Entropy','Surprise'}); title('Group average'); ylabel('ms/bit');
hold on; er = errorbar(1:2,mean(bitsl),err); er.LineStyle = 'none'; er.Color = [0 0 0];
hold off


