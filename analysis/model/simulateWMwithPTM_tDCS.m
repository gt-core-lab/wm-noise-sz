clear all;
close all;

%% paths and data

% paths to PTM
paths.local = fullfile('/Users', 'wpark78', 'Documents', 'code', 'ptm-wm-sz');
paths.csv   = fullfile(paths.local, 'csv');
ptm_csv_path = fullfile(paths.csv,'tdcs.csv');

% paths to subject data files  
paths.remote = fullfile('/Users', 'wpark78', 'Dropbox', 'WoonJuPark', 'Projects', 'active');
paths.project = fullfile(paths.remote, 'vandy-wm', 'tdcs-wm');
paths.data = fullfile(paths.project, 'Sham');
files = dir(fullfile(paths.data, '*_shortNEW.csv'));

% get subject list
subList = {};
for i = 1:numel(files)
    fname = files(i).name;
    tok = regexp(fname, '(\d{5})-33_VWMp_shortNEW\.csv', 'tokens', 'once');
    if ~isempty(tok)
        subList{end+1} = tok{1}; %#ok<SAGROW>
    end
end

% get WM trial data (errors + set sizes) via user function getData(paths, subid)
Trials = table();
for i = 1:numel(subList)
    subid = subList{i};
    [errs, Ns] = getData_tDCS(paths, subid);
    nTrials = length(errs);
    Trials = [Trials; table(repmat(string(subid), nTrials, 1), errs, Ns, ...
        'VariableNames', {'SubID', 'ErrorRad', 'N'})];
end

trialSubIDs = unique(Trials.SubID, 'stable');
if ~exist('trialSubIDs','var'), error('trialSubIDs not defined.'); end

% group labels
groupLabels = strings(size(trialSubIDs));
for i = 1:numel(trialSubIDs)
    groupLabels(i) = "SZ";
end

%% model settings 

% spatial layout
numLocs         = 8;          % total possible stimulus locations on a circle
neuronsPerPool  = 30;         % neurons per location/pool 
prefAngles      = linspace(0, 2*pi, neuronsPerPool+1); prefAngles(end) = [];

% WM
setSizes        = [2 4];
nTrialsPerSet   = 10000;
delaySec        = 1;
probeRandomItem = true;
rng(2025);

% orientation population
kappa_tune_base = 10;         % base von Mises tuning sharpness
baseline_rate   = 1;

% global set-size effect on encoding (resource division):
gamma_setsize   = 1;          % 1=classic 1/K division; 0=no division

% noise 
poisson_like_mult = 0.15;     % multiplicative noise ~ variance proportional to rate

% maintenance (diffusion)
diffusion_extGain = 0.60;     % inflation by external filtering (phi)
diffusion_setGain = 0.35;     % inflation per extra item

% tuning
tune_filter_gain  = 0.7;      % tuning weight (external noise filtering)

% PTM handling
useZ = true;                  % true if z-scoring PTM across matched subs

% calibrated mapping parameters (from simulateWMwithPTM)
map_sigma_int_base = 0.2015;
map_sigma_int_slope = 0.0635;
map_ext_filt_base = 1.5105;
map_ext_filt_slope = 1.0007;
gain_base = 51.6253;
diffusion_base = 0.0944; 

% raw scaling (if useZ = false);
raw_sigma_scale     = 0.05;
raw_phi_scale       = 1.0;

%% read PTM

PTMtdcs = readtable(ptm_csv_path);

% ---------- choose stimulation site here ----------
stimSite = 'FCZ';     % 'PZ' or 'FCZ'
% --------------------------------------------------

% build column names for the chosen site
stimNaCol   = sprintf('%s_Na',   stimSite);   % e.g., 'PZ_Na' or 'FCZ_Na'
stimAfCol   = sprintf('%s_Af',   stimSite);   % e.g., 'PZ_Af' or 'FCZ_Af'
stimVar2Col = sprintf('%s_Var2', stimSite);   % e.g., 'PZ_Var2' or 'FCZ_Var2'
stimVar4Col = sprintf('%s_Var4', stimSite);   % e.g., 'PZ_Var4' or 'FCZ_Var4'
shPhiCol = 'SH_Af';
shNaCol   = 'SH_Na';
shVar2Col = 'Sh_Var2';
shVar4Col = 'Sh_Var4';

% check that all required columns exist
neededCols = {stimNaCol, stimAfCol, stimVar2Col, stimVar4Col, ...
              shNaCol, shPhiCol, shVar2Col, shVar4Col};

if ~all(ismember(neededCols, PTMtdcs.Properties.VariableNames))
    missing = neededCols(~ismember(neededCols, PTMtdcs.Properties.VariableNames));
    error('PTMtdcs table is missing required columns: %s', strjoin(missing, ', '));
end

subIDs = string(PTMtdcs.SID);
nSub   = numel(subIDs);

% long format table
conds    = ["Sh", stimSite];   % second condition label = 'PZ' or 'FCZ'
nCond    = numel(conds);
nModel   = nSub * nCond;
nK       = numel(setSizes);

PTM_long = table( ...
    strings(nModel,1), ...  % SubID
    strings(nModel,1), ...  % Cond
    nan(nModel,1), ...      % sigma_add
    nan(nModel,1), ...      % phi
    nan(nModel,1), ...      % invK2
    nan(nModel,1), ...      % invK4
    'VariableNames', {'SubID','Cond','sigma_add','phi','invK2','invK4'});

idx = 0;
for s = 1:nSub
    sid = subIDs(s);

    % -------- Sham row --------
    idx = idx + 1;
    PTM_long.SubID(idx)     = sid;
    PTM_long.Cond(idx)      = "Sh";
    PTM_long.sigma_add(idx) = PTMtdcs.(shNaCol)(s);
    PTM_long.phi(idx)       = PTMtdcs.(shPhiCol)(s);
    PTM_long.invK2(idx)     = PTMtdcs.(shVar2Col)(s);
    PTM_long.invK4(idx)     = PTMtdcs.(shVar4Col)(s);

    % -------- stimSite row (PZ or FCZ) --------
    idx = idx + 1;
    PTM_long.SubID(idx)     = sid;
    PTM_long.Cond(idx)      = string(stimSite);   % "PZ" or "FCZ"
    PTM_long.sigma_add(idx) = PTMtdcs.(stimNaCol)(s);
    PTM_long.phi(idx)       = PTMtdcs.(stimAfCol)(s);
    PTM_long.invK2(idx)     = PTMtdcs.(stimVar2Col)(s);
    PTM_long.invK4(idx)     = PTMtdcs.(stimVar4Col)(s);
end

% Empirical 1/kappa matrix: [nModel x nK]
empirical_invK = [PTM_long.invK2, PTM_long.invK4];


%% z scoring

valid = ~isnan(PTM_long.sigma_add) & ~isnan(PTM_long.phi);

z_sigma = nan(nModel,1);
z_phi   = nan(nModel,1);
z_dom   = nan(nModel,1);

if useZ
    mu_s = mean(PTM_long.sigma_add(valid),'omitnan');
    sd_s = std( PTM_long.sigma_add(valid),'omitnan');
    mu_f = mean(PTM_long.phi(valid),'omitnan');
    sd_f = std( PTM_long.phi(valid),'omitnan');

    z_sigma(valid) = (PTM_long.sigma_add(valid) - mu_s) ./ max(sd_s,eps);
    z_phi(valid)   = (PTM_long.phi(valid)        - mu_f) ./ max(sd_f,eps);
else
    z_sigma = raw_sigma_scale * PTM_long.sigma_add;
    z_phi   = raw_phi_scale   * PTM_long.phi;
end

% Dominant noise per subject×cond (which dimension deviates more)
for i = 1:nModel
    if ~valid(i), continue; end
    if isnan(z_sigma(i)) && isnan(z_phi(i))
        z_dom(i) = NaN;
    else
        z_dom(i) = z_sigma(i);
    end
end

%% map to model parameters

[sigma_map_all, tune_map_all, decay_map_all] = map_subject_params( ...
    valid, useZ, ...
    map_sigma_int_base, map_sigma_int_slope, ...
    map_ext_filt_base,  map_ext_filt_slope, ...
    raw_sigma_scale, raw_phi_scale, ...
    z_sigma, z_phi, z_dom, PTM_long);

%% simulate WM

fprintf('Simulating WM for %d SZ subjects × 2 conds (Sh, PZ).\n', nSub);

pred_invK = simulate_pred_invK_core( ...
    setSizes, nTrialsPerSet, neuronsPerPool, prefAngles, ...
    kappa_tune_base, baseline_rate, gain_base, gamma_setsize, ...
    diffusion_base, diffusion_extGain, diffusion_setGain, delaySec, ...
    poisson_like_mult, probeRandomItem, ...
    tune_filter_gain, ...
    valid, sigma_map_all, tune_map_all, decay_map_all, ...
    false, false, true);  % useNoise, useTune, useDecay

%% condition differences per subject

% Indices for K=2 and K=4 in setSizes
idxK2 = find(setSizes==2, 1);
idxK4 = find(setSizes==4, 1);

noise_Sh      = nan(nSub,1);
noise_stim    = nan(nSub,1);    % PZ or FCZ
emp_Sh_K2     = nan(nSub,1); emp_stim_K2  = nan(nSub,1);
emp_Sh_K4     = nan(nSub,1); emp_stim_K4  = nan(nSub,1);
pred_Sh_K2    = nan(nSub,1); pred_stim_K2 = nan(nSub,1);
pred_Sh_K4    = nan(nSub,1); pred_stim_K4 = nan(nSub,1);

for s = 1:nSub
    sid  = subIDs(s);
    idxSh   = find(PTM_long.SubID==sid & PTM_long.Cond=="Sh");
    idxStim = find(PTM_long.SubID==sid & PTM_long.Cond==string(stimSite));

    if isempty(idxSh) || isempty(idxStim)
        warning('Missing Sh or %s row for subject %s', stimSite, sid);
        continue;
    end

    noise_Sh(s)   = PTM_long.sigma_add(idxSh); % z_dom? z_sigma?
    noise_stim(s) = PTM_long.sigma_add(idxStim);

    emp_Sh_K2(s)    = empirical_invK(idxSh,   idxK2);
    emp_stim_K2(s)  = empirical_invK(idxStim, idxK2);
    emp_Sh_K4(s)    = empirical_invK(idxSh,   idxK4);
    emp_stim_K4(s)  = empirical_invK(idxStim, idxK4);

    pred_Sh_K2(s)   = pred_invK(idxSh,   idxK2);
    pred_stim_K2(s) = pred_invK(idxStim, idxK2);
    pred_Sh_K4(s)   = pred_invK(idxSh,   idxK4);
    pred_stim_K4(s) = pred_invK(idxStim, idxK4);
end

% Δ noise (stimSite - Sh); Δ WM (stimSite - Sh) (positive = more noise, worse precision)
dNoise_dom  = noise_stim   - noise_Sh;

dEmp_K2     = emp_stim_K2  - emp_Sh_K2;
dEmp_K4     = emp_stim_K4  - emp_Sh_K4;

dPred_K2    = pred_stim_K2 - pred_Sh_K2;
dPred_K4    = pred_stim_K4 - pred_Sh_K4;

%% plot predicted performance vs set size (individual + mean)

conds = ["Sh", string(stimSite)];
cmap  = lines(numel(conds));

figure('Color','w','Position',[100 100 640 420]); hold on; box on; grid on;

for ci = 1:numel(conds)
    cond = conds(ci);
    clr  = cmap(ci,:);

    % rows for this condition
    condMask = (PTM_long.Cond == cond);
    subIDs_cond = unique(PTM_long.SubID(condMask),'stable');

    % ---- individual subject lines ----
    for s = 1:numel(subIDs_cond)
        sid  = subIDs_cond(s);
        ridx = find(PTM_long.SubID == sid & PTM_long.Cond == cond);
        if isempty(ridx), continue; end
        % each ridx should correspond to a single row in pred_invK
        plot(setSizes, pred_invK(ridx,:), '-', ...
            'Color', clr, 'LineWidth', 0.5, 'HandleVisibility','off');
    end

    % ---- condition mean ± SEM ----
    preds_cond = pred_invK(condMask, :);   % [nSub x nK] for this cond
    meanPred   = mean(preds_cond, 1, 'omitnan');
    semPred    = std(preds_cond, 0, 1, 'omitnan') ./ sqrt(sum(condMask));

    errorbar(setSizes, meanPred, semPred, '-o', ...
        'Color', clr, 'MarkerFaceColor', clr, ...
        'LineWidth', 2, 'DisplayName', sprintf('%s mean', cond));
end

set(gca,'XTick', setSizes);
xlabel('Set size (K)');
ylabel('Predicted 1/\kappa');
title(sprintf('Predicted WM precision vs set size (Sh vs %s)', stimSite));
legend('Location','northwest');

%% Average empirical & predicted WM across set sizes (2 & 4)

% dEmp_avg  = nanmean([dEmp_K2,  dEmp_K4], 2);   % Δ empirical = stimSite - Sh
% dPred_avg = nanmean([dPred_K2, dPred_K4], 2);  % Δ predicted = stimSite - Sh
dEmp_avg = 1./((1./emp_stim_K2 + 1./emp_stim_K4)./2) - 1./((1./emp_Sh_K2+1./emp_Sh_K4)./2);
dPred_avg = 1./((1./pred_stim_K2 + 1./pred_stim_K4)./2) - 1./((1./pred_Sh_K2+1./pred_Sh_K4)./2);

okNoise = isfinite(dNoise_dom);
okEmp   = isfinite(dEmp_avg);
okPred  = isfinite(dPred_avg);

%% 1) Correlation: Δ noise vs Δ empirical WM
mask = okNoise & okEmp;

[r1, p1] = corr(dNoise_dom(mask), dEmp_avg(mask), 'Type','Pearson');

fprintf('\n1) ΔNoise → ΔEmpirical (%s - Sh)\n', stimSite);
fprintf('r = %.3f, p = %.3g, N = %d\n', r1, p1, sum(mask));

%% 2) Correlation: Δ noise vs Δ predicted WM
mask = okNoise & okPred;

[r2, p2] = corr(dNoise_dom(mask), dPred_avg(mask), 'Type','Pearson');

fprintf('\n2) ΔNoise → ΔPredicted (%s - Sh)\n', stimSite);
fprintf('r = %.3f, p = %.3g, N = %d\n', r2, p2, sum(mask));

%% 3) Correlation: Δ predicted WM vs Δ empirical WM
mask = okPred & okEmp;

[r3, p3] = corr(dPred_avg(mask), dEmp_avg(mask), 'Type','Pearson');

fprintf('\n3) ΔPredicted → ΔEmpirical (%s - Sh)\n', stimSite);
fprintf('r = %.3f, p = %.3g, N = %d\n', r3, p3, sum(mask));

figure; scatter(dNoise_dom(mask), dEmp_avg(mask), 80, 'filled'); lsline;
title(sprintf('ΔNoise vs ΔEmpirical (%s - Sh)', stimSite));
xlabel('ΔNoise'); ylabel('ΔEmp');

figure; scatter(dNoise_dom(mask), dPred_avg(mask), 80, 'filled'); lsline;
title(sprintf('ΔNoise vs ΔPredicted (%s - Sh)', stimSite));
xlabel('ΔNoise'); ylabel('ΔPred');

figure; scatter(dEmp_avg(mask), dPred_avg(mask), 80, 'filled'); lsline;
title(sprintf('ΔEmpirical vs ΔPredicted (%s - Sh)', stimSite));
xlabel('ΔEmp'); ylabel('ΔPred');
