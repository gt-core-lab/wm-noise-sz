clear all;
close all;

%% paths and data

% paths to PTM
paths.local = fullfile('/Users', 'wpark78', 'Documents', 'code', 'ptm-wm-sz');
paths.csv   = fullfile(paths.local, 'csv');
ptm_csv_path = fullfile(paths.csv,'ptm-wm.csv');

% paths to subject data files  
paths.remote = fullfile('/Users', 'wpark78', 'Dropbox', 'WoonJuPark', 'Projects', 'active');
paths.project = fullfile(paths.remote, 'vandy-wm', 'wm');
paths.data = fullfile(paths.project, '@Data');
files = dir(fullfile(paths.data, '*_VWMp.csv'));

% output folder for CSV exports
paths.out = fullfile(paths.local, 'csv');   
if ~exist(paths.out, 'dir')
    mkdir(paths.out);
end

% get subject list
subList = {};
for i = 1:numel(files)
    fname = files(i).name;
    tok = regexp(fname, '([23]\d{3}[A-Za-z]{2})_VWMp', 'tokens', 'once');
    if ~isempty(tok)
        subList{end+1} = tok{1}; %#ok<SAGROW>
    end
end

% get WM trial data (errors + set sizes) 
Trials = table();
for i = 1:numel(subList)
    subid = subList{i};
    [errs, Ns] = getData(paths, subid);
    nTrials = length(errs);
    Trials = [Trials; table(repmat(string(subid), nTrials, 1), errs, Ns, ...
        'VariableNames', {'SubID', 'ErrorRad', 'N'})];
end

trialSubIDs = unique(Trials.SubID, 'stable');
if ~exist('trialSubIDs','var'), error('trialSubIDs not defined.'); end

% group labels
groupLabels = strings(size(trialSubIDs));
for i = 1:numel(trialSubIDs)
    idnum = str2double(extractBefore(trialSubIDs(i),5));
    if idnum < 3000
        groupLabels(i) = "SZ";
    else
        groupLabels(i) = "NT";
    end
end

%% model settings 

% spatial layout
numLocs         = 8;          % total possible stimulus locations on a circle
neuronsPerPool  = 30;         % neurons per location/pool 
prefAngles      = linspace(0, 2*pi, neuronsPerPool+1); prefAngles(end) = [];

% WM
setSizes        = [1 2 4];
nTrialsPerSet   = 1000;
delaySec        = 1;
probeRandomItem = true;
rng(2025);

% orientation population
kappa_tune_base = 10;         % base von Mises tuning sharpness
gain_base       = 40;         % base response gain (set size = 1)
baseline_rate   = 1;

% global set-size effect on encoding (resource division):
gamma_setsize   = 1;          % 1=classic 1/K division; 0=no division

% noise 
poisson_like_mult = 0.15;     % multiplicative noise ~ variance proportional to rate

% maintenance (diffusion)
diffusion_base    = 0.015;    % rad^2/s at K=1, ext_filter=0
diffusion_extGain = 0.60;     % inflation by external filtering (phi)
diffusion_setGain = 0.35;     % inflation per extra item

% tuning
tune_filter_gain  = 0.7;      % tuning weight (external noise filtering)

% PTM handling
useZ = true;                  % true if z-scoring PTM across matched subs

% mapping parameters to MODEL UNITS
% NOTE: sigma_int now means ANGULAR SD (radians) added AFTER decoding.
map_sigma_int_base  = 0.25;   % ~8.6°; re-calibrated by optimizer anyway
map_sigma_int_slope = 0.06;   % z-score slope
map_ext_filt_base   = 1.5;
map_ext_filt_slope  = 1.0;

% if not z-scoring (useZ=false), optionally scale raw PTM to model units:
raw_sigma_scale     = 0.05;   % sigma_add * scale -> model angular SD (radians)
raw_phi_scale       = 1.0;    % phi * scale -> model external filtering

%% read PTM

PTMtbl = getDataPTM(ptm_csv_path, trialSubIDs);  % SubID, sigma_add, phi, k1var,k2var,k4var

% align to trialSubIDs order
[tf, loc] = ismember(string(trialSubIDs(:)), PTMtbl.SubID);
if ~all(tf)
    warning('Missing PTM rows for %d SubID(s). Those rows will be NaN and skipped.', sum(~tf));
end
PTMtbl = PTMtbl(loc,:);  % align rows

%% z-scoring

valid = tf & ~isnan(PTMtbl.sigma_add) & ~isnan(PTMtbl.phi);
gl     = string(groupLabels(:));
groups_u = unique(gl,'stable');

% GLOBAL z across everyone
z_sigma_global = nan(height(PTMtbl),1);
z_phi_global   = nan(height(PTMtbl),1);
if useZ
    mu_sig_g = mean(PTMtbl.sigma_add(valid),'omitnan');
    sd_sig_g = std( PTMtbl.sigma_add(valid),'omitnan');
    mu_phi_g = mean(PTMtbl.phi(valid),'omitnan');
    sd_phi_g = std( PTMtbl.phi(valid),'omitnan');

    z_sigma_global(valid) = (PTMtbl.sigma_add(valid) - mu_sig_g) ./ max(sd_sig_g, eps);
    z_phi_global(valid)   = (PTMtbl.phi(valid)       - mu_phi_g) ./ max(sd_phi_g, eps);
else
    % if not z-scoring, just put raw scaled values into the "global" slots
    z_sigma_global = raw_sigma_scale * PTMtbl.sigma_add;
    z_phi_global   = raw_phi_scale   * PTMtbl.phi;
end

% WITHIN-GROUP z (for deciding dominance only) 
z_sigma_wg = nan(height(PTMtbl),1);
z_phi_wg   = nan(height(PTMtbl),1);
if useZ
    for gi = 1:numel(groups_u)
        gmask = valid & (gl == groups_u(gi));
        if ~any(gmask), continue; end
        mu_s = mean(PTMtbl.sigma_add(gmask),'omitnan');
        sd_s = std( PTMtbl.sigma_add(gmask),'omitnan');
        mu_f = mean(PTMtbl.phi(gmask),'omitnan');
        sd_f = std( PTMtbl.phi(gmask),'omitnan');
        z_sigma_wg(gmask) = (PTMtbl.sigma_add(gmask) - mu_s) ./ max(sd_s, eps);
        z_phi_wg(gmask)   = (PTMtbl.phi(gmask)       - mu_f) ./ max(sd_f, eps);
    end
else
    % when not using z, dominance by raw scaled values within-group
    for gi = 1:numel(groups_u)
        gmask = valid & (gl == groups_u(gi));
        if ~any(gmask), continue; end
        s_raw = raw_sigma_scale * PTMtbl.sigma_add(gmask);
        f_raw = raw_phi_scale   * PTMtbl.phi(gmask);
        % store them as "wg" just to reuse the same selection logic below
        z_sigma_wg(gmask) = s_raw;
        z_phi_wg(gmask)   = f_raw;
    end
end

% decide the winner using WITHIN-GROUP z, but take magnitude from GLOBAL z
z_sigma = z_sigma_global;  % these are the values to feed to mapping for sigma
z_phi   = z_phi_global;    % and for phi (tuning)
z_dom   = nan(height(PTMtbl),1);

for i = 1:height(PTMtbl)
    if ~valid(i), continue; end
    a_w = z_sigma_wg(i); b_w = z_phi_wg(i);

    if isnan(a_w) && isnan(b_w)
        z_dom(i) = NaN;                      % no info
    elseif isnan(b_w) || (a_w >= b_w)
        z_dom(i) = z_sigma_global(i);        % sigma wins -> take GLOBAL sigma value
    else
        z_dom(i) = z_phi_global(i);          % phi wins -> take GLOBAL phi value
    end
end

% empirical 1/kappa matrix (nSub x 3)
empirical_invK = [PTMtbl.k1var, PTMtbl.k2var, PTMtbl.k4var];

%% parameter calibration

target_invK   = [0.20, 0.32, 0.55]; % mean empirical 1/kappa
calib_nTrials = 1000;            
orig_nTrials  = nTrialsPerSet;
nTrialsPerSet = calib_nTrials;

% pack params (optimize logs to keep positivity)
p0 = log([ ...
    map_sigma_int_base,  ... % p1
    map_sigma_int_slope, ... % p2
    map_ext_filt_base,   ... % p3
    map_ext_filt_slope,  ... % p4
    gain_base,           ... % p5
    diffusion_base       ... % p6
]);

w_mean  = 1.0;   % match mean 1/kappa per K
w_slope = 1.0;   % slope toward unity (groupwise per K)

obj = @(p) calib_loss_map_clean(p, ...
    trialSubIDs, PTMtbl, valid, useZ, z_sigma, z_phi, z_dom, ...
    setSizes, neuronsPerPool, kappa_tune_base, baseline_rate, ...
    poisson_like_mult, gamma_setsize, diffusion_extGain, diffusion_setGain, ...
    delaySec, probeRandomItem, raw_sigma_scale, raw_phi_scale, ...
    target_invK, groupLabels, w_mean, w_slope, ...
    tune_filter_gain);

opts = optimset('Display','iter','TolX',1e-2,'TolFun',1e-3,'MaxIter',30);
phat = fminsearch(obj, p0, opts);

% unpack
map_sigma_int_base  = exp(phat(1));
map_sigma_int_slope = exp(phat(2));
map_ext_filt_base   = exp(phat(3));
map_ext_filt_slope  = exp(phat(4));
gain_base           = exp(phat(5));
diffusion_base      = exp(phat(6));

% restore full trials
nTrialsPerSet = orig_nTrials;

fprintf('Calibrated map:\n  sigma_base=%.3f, sigma_slope=%.3f, phi_base=%.3f, phi_slope=%.3f,\n  gain_base=%.3f, diffusion_base=%.4f\n', ...
    map_sigma_int_base, map_sigma_int_slope, map_ext_filt_base, map_ext_filt_slope, gain_base, diffusion_base);

%% main simulation 

nSub = numel(trialSubIDs);
nSet = numel(setSizes);

fprintf('Simulating %d subs; %d locations; %d neurons/pool; %d trials/cond.\n', ...
        nSub, numLocs, neuronsPerPool, nTrialsPerSet);

% subject-mapped parameters
[sigma_map_all, tune_map_all, decay_map_all] = map_subject_params( ...
    valid, useZ, ...
    map_sigma_int_base, map_sigma_int_slope, map_ext_filt_base, map_ext_filt_slope, ...
    raw_sigma_scale, raw_phi_scale, ...
    z_sigma, z_phi, z_dom, PTMtbl);

% run unified core simulator (full model ON)
pred_invK = simulate_pred_invK_core( ...
    setSizes, nTrialsPerSet, neuronsPerPool, prefAngles, ...
    kappa_tune_base, baseline_rate, gain_base, gamma_setsize, ...
    diffusion_base, diffusion_extGain, diffusion_setGain, delaySec, ...
    poisson_like_mult, probeRandomItem, ...
    tune_filter_gain, ...
    valid, sigma_map_all, tune_map_all, decay_map_all, ...
    true, true, true);

% extract by-products for reporting (R and kappa)
kappa_est = 1 ./ pred_invK;
Rvals = nan(nSub,nSet);
for s=1:nSub
    for si=1:nSet
        kappa = kappa_est(s,si);
        Rvals(s,si) = NaN;
    end
end

%% export CSVs

% one row per subject × set size)
nSubX = numel(trialSubIDs);
nSetX = numel(setSizes);

SubID_long   = repelem(string(trialSubIDs(:)), nSetX, 1);
Group_long   = repelem(string(groupLabels(:)), nSetX, 1);
Valid_long   = repelem(valid(:), nSetX, 1);
SetSize_long = repmat(setSizes(:), nSubX, 1); SetSize_long = SetSize_long(:);

Emp_long  = empirical_invK.';  Emp_long  = Emp_long(:);
Pred_long = pred_invK.';       Pred_long = Pred_long(:);

Tlong = table( ...
    SubID_long, Group_long, Valid_long, SetSize_long, ...
    Emp_long, Pred_long, ...
    'VariableNames', { ...
        'SubID','Group','Valid','SetSize', ...
        'Empirical_invK','Predicted_invK', ...
    });

out_long = fullfile(paths.out, 'wm_model_main.csv');
writetable(Tlong, out_long);
fprintf('Saved: %s\n', out_long);

%% plot group means

figure('Color','w','Position',[100 100 800 500]); hold on; 
cols = lines(numel(unique(groupLabels)));

uniqueGroups = unique(groupLabels, 'stable');
for gi = 1:numel(uniqueGroups)
    g = uniqueGroups(gi);
    idx = strcmp(groupLabels, g);
    plot(setSizes, pred_invK(idx,:)', '-', 'Color', cols(gi,:), 'LineWidth',0.5, 'HandleVisibility','off');
    meanK = mean(pred_invK(idx,:), 'omitnan');
    semK  = std(pred_invK(idx,:),0,'omitnan') / sqrt(sum(idx));
    errorbar(setSizes, meanK, semK, '-o', 'Color', cols(gi,:), ...
        'LineWidth',2, 'MarkerFaceColor',cols(gi,:), 'DisplayName',g);
end
set(gca,'XScale','log'); xticks(setSizes);
xlabel('Set size (K)'); ylabel('Predicted 1/\kappa');
title('Predicted WM Precision (Group Means ± SEM, log-scaled x)');
legend('Location','northeast'); grid on; box on; hold off;

%% correlation by group

if exist('empirical_invK','var') && ~isempty(empirical_invK)
    Klist = setSizes(:)'; nSet = numel(Klist);
    cmap  = lines(nSet);
    glower = lower(string(groupLabels));
    wantGroups = ["nt","sz"];
    for gi = 1:numel(wantGroups)
        gname = wantGroups(gi);
        gmask = (glower == gname);
        if ~any(gmask)
            warning('No subjects for %s.', upper(gname)); continue;
        end
        figure('Color','w','Position',[200+520*(gi-1) 200 520 440]); hold on; grid on;
        for si = 1:nSet
            xi = empirical_invK(gmask, si);
            yi = pred_invK(gmask, si);
            ok = isfinite(xi) & isfinite(yi);
            scatter(xi(ok), yi(ok), 80, 'filled', ...
                'MarkerFaceColor', cmap(si,:), 'MarkerFaceAlpha', 0.8, ...
                'DisplayName', sprintf('K=%d', Klist(si)));
            if sum(ok) >= 3
                P  = polyfit(xi(ok), yi(ok), 1);
                xx = linspace(min(xi(ok)), max(xi(ok)), 50);
                plot(xx, polyval(P, xx), '--', 'Color', cmap(si,:), 'LineWidth', 1.5);
                [r,p] = corr(xi(ok), yi(ok), 'Type','Pearson');
                fprintf('%s | K=%d: r=%.3f, p=%.3g, N=%d\n', upper(gname), Klist(si), r, p, sum(ok));
            end
        end
        xlabel('Empirical 1/\kappa'); ylabel('Predicted 1/\kappa');
        title(sprintf('Empirical vs Predicted 1/\\kappa — %s', upper(gname)));
        legend('Location','northwestoutside'); axis square; box on;
    end
else
    warning('empirical_invK not available — skipping correlation plots.');
end

%% mixed effects model 

[nSubX, nSetX] = size(empirical_invK);
Klist = setSizes(:)';
SubID = repelem((1:nSubX)', nSetX, 1);
Group = repmat(string(groupLabels(:)), nSetX, 1);
SetSize = repmat(Klist, nSubX, 1); SetSize = SetSize(:);
Empirical = empirical_invK(:);
Predicted = pred_invK(:);

ok = isfinite(Empirical) & isfinite(Predicted);
T = table(SubID(ok), Group(ok), SetSize(ok), Empirical(ok), Predicted(ok), ...
    'VariableNames', {'SubID','Group','SetSize','Empirical','Predicted'});

groups = unique(T.Group, 'stable');
models = cell(numel(groups),1);
for gi = 1:numel(groups)
    gname = groups(gi);
    Tsub  = T(strcmp(T.Group, gname), :);
    fprintf('\nGroup: %s  (n=%d rows)\n', gname, height(Tsub));
    mdl = fitlme(Tsub, 'Predicted ~ Empirical + (1|SetSize) + (1|SubID)');
    models{gi} = mdl;
    disp(anova(mdl)); disp(mdl.Coefficients);
    ci = coefCI(mdl); b  = mdl.Coefficients;
    slopeRow = strcmp(b.Name,'Empirical');
    slope = b.Estimate(slopeRow); slopeCI = ci(slopeRow,:);
    fprintf('Slope = %.3f  (95%% CI [%.3f, %.3f])\n', slope, slopeCI(1), slopeCI(2));
end

%% export to csv

gridN = 120;                 % points per line
alphaCI = 0.32;              % 68% CI -> alpha = 1-0.68 = 0.32
useMarginal = false;         % marginal prediction (fixed effects only)

rows = table();

groups = unique(T.Group, 'stable');
Klist = unique(T.SetSize, 'stable');    % should be [1 2 4] (as numeric)

for gi = 1:numel(groups)
    gname = groups(gi);
    Tsub  = T(strcmp(T.Group, gname), :);

    % Fit the same model 
    mdl = fitlme(Tsub, 'Predicted ~ Empirical + (1|SetSize) + (1|SubID)');

    for ki = 1:numel(Klist)
        k = Klist(ki);

        Tk = Tsub(Tsub.SetSize == k, :);
        if height(Tk) < 3, continue; end

        xmin = min(Tk.Empirical);
        xmax = max(Tk.Empirical);
        xgrid = linspace(xmin, xmax, gridN)';

        dummySub = Tk.SubID(1);

        newT = table( ...
            repmat(dummySub, gridN, 1), ...   % SubID
            repmat(k, gridN, 1), ...          % SetSize
            xgrid, ...                        % Empirical
            'VariableNames', {'SubID','SetSize','Empirical'});

        % Predict (marginal = fixed effects only) + CI if available
        try
            if useMarginal
                [yhat, yCI] = predict(mdl, newT, 'Conditional', false, 'Alpha', alphaCI);
            else
                [yhat, yCI] = predict(mdl, newT, 'Alpha', alphaCI);
            end
            ci_lo = yCI(:,1);
            ci_hi = yCI(:,2);
        end

        Ti = table( ...
            repmat(string(gname), gridN, 1), ...
            repmat(k, gridN, 1), ...
            xgrid, yhat, ci_lo, ci_hi, ...
            'VariableNames', {'Group','SetSize','Empirical_grid','Predicted_hat','CI_lo','CI_hi'});

        rows = [rows; Ti]; %#ok<AGROW>
    end
end

out_lines = fullfile(paths.out, 'wm_model_full_LME_lines.csv');
writetable(rows, out_lines);
fprintf('Saved LME lines CSV: %s\n', out_lines);

%% permutation test

Nperm = 1000; 
perm_nTrialsPerSet = min(400, nTrialsPerSet); 
rng(2025);
groups_u = unique(T.Group,'stable');
score_real = 0;
for gi = 1:numel(groups_u)
    mdl = models{gi};
    try
        A = anova(mdl); row = strcmp(string(A.Term),"Empirical");
        F = A.FStat(row); df1=A.DF1(row); df2=A.DF2(row);
        s = (F*df1)/(F*df1+df2);
    catch
        b = mdl.Coefficients; j = strcmp(b.Name,'Empirical');
        s = max(0,b.Estimate(j))*max(0,b.tStat(j));
    end
    score_real = score_real + s;
end
fprintf('\nREAL mixed-model score (sum partial R^2) = %.4f\n', score_real);

scores_perm = nan(Nperm,1);
lbl_lower = lower(string(groupLabels));
grp_tags = unique(lbl_lower,'stable');

for pidx = 1:Nperm

    if useZ
        % Shuffle RAW PTM within group
        PTM_s = PTMtbl.sigma_add;
        PTM_f = PTMtbl.phi;
        for gi = 1:numel(grp_tags)
            idx  = find(lbl_lower == grp_tags(gi));
            perm = idx(randperm(numel(idx)));
            PTM_s(idx) = PTMtbl.sigma_add(perm);
            PTM_f(idx) = PTMtbl.phi(perm);
        end
        PTMtbl_sh = PTMtbl;
        PTMtbl_sh.sigma_add = PTM_s;
        PTMtbl_sh.phi       = PTM_f;
    
        % Recompute GLOBAL z on shuffled raws
        vmask = ~isnan(PTM_s) & ~isnan(PTM_f);
        mu_sg = mean(PTM_s(vmask),'omitnan'); sd_sg = std(PTM_s(vmask),'omitnan');
        mu_fg = mean(PTM_f(vmask),'omitnan'); sd_fg = std(PTM_f(vmask),'omitnan');
        z_sigma_global_sh = (PTM_s - mu_sg) ./ max(sd_sg, eps);
        z_phi_global_sh   = (PTM_f - mu_fg) ./ max(sd_fg, eps);
    
        % Recompute WITHIN-GROUP z (for winner)
        z_sigma_wg_sh = nan(size(PTM_s));
        z_phi_wg_sh   = nan(size(PTM_f));
        for gi = 1:numel(grp_tags)
            idx = find(lbl_lower == grp_tags(gi));
            mu_s = mean(PTM_s(idx),'omitnan'); sd_s = std(PTM_s(idx),'omitnan');
            mu_f = mean(PTM_f(idx),'omitnan'); sd_f = std(PTM_f(idx),'omitnan');
            z_sigma_wg_sh(idx) = (PTM_s(idx) - mu_s) ./ max(sd_s, eps);
            z_phi_wg_sh(idx)   = (PTM_f(idx) - mu_f) ./ max(sd_f, eps);
        end
    
        % Winner by within-group; magnitude from GLOBAL
        z_dom_sh = nan(size(PTM_s));
        for i = 1:numel(PTM_s)
            aw = z_sigma_wg_sh(i); bw = z_phi_wg_sh(i);
            if isnan(aw) && isnan(bw)
                z_dom_sh(i) = NaN;
            elseif isnan(bw) || (aw >= bw)
                z_dom_sh(i) = z_sigma_global_sh(i);
            else
                z_dom_sh(i) = z_phi_global_sh(i);
            end
        end
    
        % Map to model space (three outputs)
        [sigma_sh, tune_sh, decay_sh] = map_subject_params( ...
            valid, true, ...
            map_sigma_int_base, map_sigma_int_slope, map_ext_filt_base, map_ext_filt_slope, ...
            raw_sigma_scale, raw_phi_scale, ...
            z_sigma_global_sh, z_phi_global_sh, z_dom_sh, PTMtbl_sh);
    else
        % shuffle RAW PTM within group
        PTM_s = PTMtbl.sigma_add;
        PTM_f = PTMtbl.phi;
        for gi = 1:numel(grp_tags)
            idx  = find(lbl_lower == grp_tags(gi));
            perm = idx(randperm(numel(idx)));
            PTM_s(idx) = PTMtbl.sigma_add(perm);
            PTM_f(idx) = PTMtbl.phi(perm);
        end
        PTMtbl_sh = PTMtbl;
        PTMtbl_sh.sigma_add = PTM_s;
        PTMtbl_sh.phi       = PTM_f;

        % recompute within-group z for the SHUFFLED raws so dominance is shuffled too
        z_sigma_sh = nan(size(PTM_s));
        z_phi_sh   = nan(size(PTM_f));
        for gi = 1:numel(grp_tags)
            idx = find(lbl_lower == grp_tags(gi));
            mu_s = mean(PTM_s(idx),'omitnan'); sd_s = std(PTM_s(idx),'omitnan');
            mu_f = mean(PTM_f(idx),'omitnan'); sd_f = std(PTM_f(idx),'omitnan');
            z_sigma_sh(idx) = (PTM_s(idx) - mu_s) ./ max(sd_s, eps);
            z_phi_sh(idx)   = (PTM_f(idx) - mu_f) ./ max(sd_f, eps);
        end
        z_dom_sh = max(z_sigma_sh, z_phi_sh);

        % map to model-space arrays (useZ=false branch)
        [sigma_sh, tune_sh, decay_sh] = map_subject_params( ...
            valid, false, ...
            map_sigma_int_base, map_sigma_int_slope, ...
            map_ext_filt_base,  map_ext_filt_slope, ...
            raw_sigma_scale, raw_phi_scale, ...
            z_sigma_sh, z_phi_sh, z_dom_sh, PTMtbl_sh);
    end

    % simulate with SHUFFLED, MAPPED arrays (not z-scores)
    pred_invK_perm = simulate_pred_invK_core( ...
        setSizes, perm_nTrialsPerSet, neuronsPerPool, prefAngles, ...
        kappa_tune_base, baseline_rate, gain_base, gamma_setsize, ...
        diffusion_base, diffusion_extGain, diffusion_setGain, delaySec, ...
        poisson_like_mult, probeRandomItem, ...
        tune_filter_gain, ...
        valid, sigma_sh, tune_sh, decay_sh, ...
        true, true, true);

    % mixed models per group
    Emp_perm = empirical_invK(:); Pred_perm = pred_invK_perm(:);
    okp = isfinite(Emp_perm) & isfinite(Pred_perm);
    Tperm = table(SubID(okp), Group(okp), SetSize(okp), Emp_perm(okp), Pred_perm(okp), ...
        'VariableNames', {'SubID','Group','SetSize','Empirical','Predicted'});

    score_p = 0; gu = unique(Tperm.Group,'stable');
    for gi = 1:numel(gu)
        Tsub = Tperm(strcmp(Tperm.Group, gu(gi)), :);
        if height(Tsub) < 6, continue; end
        try
            mdlp = fitlme(Tsub, 'Predicted ~ Empirical + (1|SetSize) + (1|SubID)');
            A = anova(mdlp); row = strcmp(string(A.Term),"Empirical");
            F=A.FStat(row); df1=A.DF1(row); df2=A.DF2(row);
            s=(F*df1)/(F*df1+df2);
        catch
            b = mdlp.Coefficients; j=strcmp(b.Name,'Empirical');
            s = max(0,b.Estimate(j))*max(0,b.tStat(j));
        end
        score_p = score_p + s;
    end
    scores_perm(pidx) = score_p;
end

num_ge = sum(scores_perm >= score_real);
p_perm = (1 + num_ge) / (Nperm + 1);
fprintf('\nPermutation test:\n  Real score = %.4f\n  Mean shuffle = %.4f (SD=%.4f)\n  p = %.4f\n', ...
    score_real, mean(scores_perm,'omitnan'), std(scores_perm,'omitnan'), p_perm);

figure('Color','w','Position',[200 200 560 420]); hold on; 
histogram(scores_perm, 'Normalization','pdf');
yl = ylim; plot([score_real score_real], yl, 'r-', 'LineWidth',2);
xlabel('Shuffle mixed-model score'); ylabel('PDF');
title(sprintf('Permutation test (N=%d): p=%.4f', Nperm, p_perm));
legend({'Shuffles','Real'},'Location','best');

%% ablation

rng(250);
ablate_nTrialsPerSet = min(1000, nTrialsPerSet);

variants = { ...
    struct('name','FULL',       'useNoise',true,  'useTune',true,  'useDecay',true), ...
    struct('name','NOISE_ONLY', 'useNoise',true,  'useTune',false, 'useDecay',false), ...
    struct('name','TUNE_ONLY',  'useNoise',false, 'useTune',true,  'useDecay',false), ...
    struct('name','DRIFT_ONLY', 'useNoise',false, 'useTune',false, 'useDecay',true) ...
    };

gu = unique(T.Group,'stable'); 
ngr = numel(gu);

% Storage for predictions (for LME plots) and GOF tables
variant_preds = cell(numel(variants),1);
fit_tables    = cell(numel(variants),1);

% For GOF bar charts (RMSE)
rmse_overall = nan(numel(variants),1);
rmse_byGroup = nan(numel(variants), ngr);

for vi = 1:numel(variants)
    v = variants{vi};

    % simulate this variant
    pred_invK_var = simulate_pred_invK_core( ...
        setSizes, ablate_nTrialsPerSet, neuronsPerPool, prefAngles, ...
        kappa_tune_base, baseline_rate, gain_base, gamma_setsize, ...
        diffusion_base, diffusion_extGain, diffusion_setGain, delaySec, ...
        poisson_like_mult, probeRandomItem, ...
        tune_filter_gain, ...
        valid, sigma_map_all, tune_map_all, decay_map_all, ...
        v.useNoise, v.useTune, v.useDecay);

    variant_preds{vi} = pred_invK_var;

    % goodness-of-fit (no mixed model)
    % overall_method: 'sqmean' (default), 'mean', or 'median'
    fit_tables{vi} = score_model_fit(empirical_invK, pred_invK_var, groupLabels, setSizes, v.name, 'sqmean');

    % choose what to inspect/print
    Ti_overall = fit_tables{vi}.overall;   % aggregated across K
    Ti_perK    = fit_tables{vi}.perK;      % per-K rows
    
    % print 
    disp(Ti_perK);
    disp(Ti_overall);

    % mixed-effects model only for slope interpretation and plotting fits
    % build long table 
    [nSub2, nSet2] = size(empirical_invK);
    SubID_long   = repelem((1:nSub2)', nSet2, 1);
    Group_long   = repmat(string(groupLabels(:)), nSet2, 1);
    SetSize_long = repmat(setSizes(:)', nSub2, 1); SetSize_long = SetSize_long(:);
    Emp  = empirical_invK(:);
    Pred = pred_invK_var(:);
    okL  = isfinite(Emp) & isfinite(Pred);

    Tvar = table(SubID_long(okL), Group_long(okL), SetSize_long(okL), Emp(okL), Pred(okL), ...
        'VariableNames', {'SubID','Group','SetSize','Empirical','Predicted'});

    % fit and display LME per group
    for gi = 1:numel(gu)
        Tsub = Tvar(strcmp(Tvar.Group, gu(gi)), :);
        if height(Tsub) < 6, continue; end
        mdlv = fitlme(Tsub, 'Predicted ~ Empirical + (1|SetSize) + (1|SubID)');

        % Display results cleanly
        fprintf('\n──────────────────────────────────────────────\n');
        fprintf('Variant: %-12s | Group: %s\n', v.name, string(gu(gi)));
        disp(anova(mdlv));
        disp(mdlv.Coefficients);
        try
            b  = mdlv.Coefficients; 
            ci = coefCI(mdlv);
            [rowEmp, ~] = find_term_rows(toTable_compat(b), "Empirical");
            if any(rowEmp)
                slope   = b.Estimate(rowEmp);
                slopeCI = ci(rowEmp,:);
                fprintf('→ LME slope (Empirical) = %.3f  [%.3f, %.3f]\n', slope, slopeCI(1), slopeCI(2));
            end
        catch
        end
        fprintf('──────────────────────────────────────────────\n\n');
    end
end

% GOF bar charts (bootstrap)
variantNames = string(cellfun(@(x)x.name, variants,'uni',false)); 
boot = bootstrap_model_fit( ...
    empirical_invK, variant_preds, groupLabels, setSizes, variantNames, ...
    'NBoot', 1000, 'RMSEAgg', 'sqmean', 'RAgg', 'fisher', 'Seed', 2025, 'UseParallel', false);

% Pearson r panels (0..1)
bootplot_bars(boot.r.mean, boot.r.ci_lo, boot.r.ci_hi, boot.variants, boot.groups, ...
    'YLabel','Overall Pearson r (agg across K)', ...
    'Title','Model comparison with bootstrap CI — Pearson r', ...
    'YLim',[-0.1 1], 'BarColor',[0.65 0.85 0.65], 'Fmt','%.2f');

%% Paired bootstrap (Δr vs FULL) + export CSV

ci68 = [16 84];
ci95 = [2.5 97.5];

variants = string(boot.variants(:));
groups   = string(boot.groups(:));

r_dist = boot.r.dist;  % nVar x nGrp x NBoot

% index of FULL
iFULL = find(variants == "FULL", 1);
if isempty(iFULL), error('FULL not found in boot.variants'); end

% choose which comparisons to export 
compare = ["NOISE_ONLY","TUNE_ONLY","DRIFT_ONLY"];
iComp = nan(size(compare));
for k = 1:numel(compare)
    iComp(k) = find(variants == compare(k), 1);
    if isnan(iComp(k)), error('Variant %s not found in boot.variants', compare(k)); end
end

% build tidy long table: one row per (Group × Comparison)
nVar = numel(variants);
nGrp = numel(groups);
nC = numel(compare);
rows = nGrp * nC;

Group = strings(rows,1);
Ablation = strings(rows,1);  % interpreted as FULL - Ablation
delta_mean = nan(rows,1);
delta_ci68_lo = nan(rows,1);
delta_ci68_hi = nan(rows,1);
delta_ci95_lo = nan(rows,1);
delta_ci95_hi = nan(rows,1);
p_two = nan(rows,1);

row = 0;
for gi = 1:nGrp
    for ci = 1:nC
        row = row + 1;

        d = squeeze(r_dist(iFULL, gi, :) - r_dist(iComp(ci), gi, :)); % NBoot x 1
        d = d(isfinite(d));

        Group(row) = groups(gi);
        Ablation(row) = compare(ci);

        delta_mean(row) = mean(d, 'omitnan');
        delta_ci68_lo(row) = prctile(d, ci68(1));
        delta_ci68_hi(row) = prctile(d, ci68(2));
        delta_ci95_lo(row) = prctile(d, ci95(1));
        delta_ci95_hi(row) = prctile(d, ci95(2));

        % simple two-sided sign-based bootstrap p-value
        % (fraction of deltas <= 0, doubled)
        p_two(row) = 2 * min(mean(d <= 0), mean(d >= 0));
    end
end

Tdelta = table(Group, Ablation, delta_mean, delta_ci68_lo, delta_ci68_hi, delta_ci95_lo, delta_ci95_hi, p_two);

out_csv_delta = fullfile(paths.out, 'wm_model_ablation_paired.csv');
writetable(Tdelta, out_csv_delta);
fprintf('Saved paired-bootstrap Δr CSV: %s\n', out_csv_delta);

%% export raw r summary (68% + 95%)

r_mean = mean(r_dist, 3, 'omitnan');
r_ci68_lo = prctile(r_dist, ci68(1), 3);
r_ci68_hi = prctile(r_dist, ci68(2), 3);
r_ci95_lo = prctile(r_dist, ci95(1), 3);
r_ci95_hi = prctile(r_dist, ci95(2), 3);

Variant = strings(nVar*nGrp,1);
Group2  = strings(nVar*nGrp,1);
rMean   = nan(nVar*nGrp,1);
r68lo   = nan(nVar*nGrp,1);
r68hi   = nan(nVar*nGrp,1);
r95lo   = nan(nVar*nGrp,1);
r95hi   = nan(nVar*nGrp,1);

row = 0;
for vi = 1:nVar
    for gi = 1:nGrp
        row = row + 1;
        Variant(row) = variants(vi);
        Group2(row)  = groups(gi);
        rMean(row)   = r_mean(vi,gi);
        r68lo(row)   = r_ci68_lo(vi,gi);
        r68hi(row)   = r_ci68_hi(vi,gi);
        r95lo(row)   = r_ci95_lo(vi,gi);
        r95hi(row)   = r_ci95_hi(vi,gi);
    end
end

Tr = table(Variant, Group2, rMean, r68lo, r68hi, r95lo, r95hi, ...
    'VariableNames', {'Variant','Group','r_mean','r_ci68_lo','r_ci68_hi','r_ci95_lo','r_ci95_hi'});

out_csv_r = fullfile(paths.out, 'wm_model_ablation_bootstrap_r_summary.csv');
writetable(Tr, out_csv_r);
fprintf('Saved raw r summary CSV: %s\n', out_csv_r);


%% ablation plots (emprical vs predicted)

if exist('variant_preds','var') && ~isempty(variant_preds)
    cmapK = lines(numel(setSizes));         % consistent color per K
    gu    = unique(string(groupLabels), 'stable');

    for vi = 1:numel(variants)
        vname = variants{vi};
        PredMat = variant_preds{vi};        % [nSub x nSet]

        % Build long table for filtering/labels
        [nSub2, nSet2] = size(empirical_invK);
        SubID_long   = repelem((1:nSub2)', nSet2, 1);
        Group_long   = repmat(string(groupLabels(:)), nSet2, 1);
        SetSize_long = repmat(setSizes(:)', nSub2, 1); SetSize_long = SetSize_long(:);
        Emp  = empirical_invK(:);
        Pred = PredMat(:);
        okL  = isfinite(Emp) & isfinite(Pred);

        Tvar = table(SubID_long(okL), Group_long(okL), SetSize_long(okL), Emp(okL), Pred(okL), ...
            'VariableNames', {'SubID','Group','SetSize','Empirical','Predicted'});

        % One figure per group for this variant (with LME fit lines)
        for gi = 1:numel(gu)
            gname = gu(gi);
            Tsub  = Tvar(strcmp(Tvar.Group, gname), :);
            if height(Tsub) < 6, continue; end

            % Fit mixed-effects model (random intercept by SetSize)
            mdl = fitlme(Tsub, 'Predicted ~ Empirical + (1|SetSize) + (1|SubID)');

            % Extract slope and CI 
            b   = mdl.Coefficients;
            ci  = coefCI(mdl);
            ir  = strcmp(b.Name,'Empirical');
            slope = b.Estimate(ir);
            slopeCI = ci(ir,:);  % [low, high]

            % Plot scatter + LME prediction lines per K
            figure('Color','w','Position',[200+600*(vi-1) 140+480*(gi-1) 560 440]); 
            hold on; grid on; box on;

            for si = 1:numel(setSizes)
                k = setSizes(si);
                m = (Tsub.SetSize == k);
                xi = Tsub.Empirical(m);
                yi = Tsub.Predicted(m);

                scatter(xi, yi, 70, 'filled', ...
                        'MarkerFaceColor', cmapK(si,:), 'MarkerFaceAlpha', 0.85, ...
                        'DisplayName', sprintf('K=%d', k));

                % LME-predicted line for this K
                if numel(xi) >= 2
                    xgrid = linspace(min(xi), max(xi), 100)';
                else
                    % fallback range if too few points
                    xgrid = linspace(min(Tsub.Empirical), max(Tsub.Empirical), 100)';
                end
                % pick any existing SubID from the dataset
                dummySub = Tsub.SubID(1);

                newT = table( ...
                    repmat(k, numel(xgrid), 1), ...     % SetSize
                    xgrid, ...                          % Empirical
                    repmat(dummySub, numel(xgrid), 1), ... % SubID (dummy)
                    'VariableNames', {'SetSize','Empirical','SubID'});
                
                yhat = predict(mdl, newT);
                plot(xgrid, yhat, '--', 'Color', cmapK(si,:), 'LineWidth', 1.8);

            end

            xlabel('Empirical 1/\kappa');
            ylabel('Predicted 1/\kappa');
            title(sprintf('Ablation: %s — %s (LME slope = %.2f [%.2f, %.2f])', ...
                  vname, upper(string(gname)), slope, slopeCI(1), slopeCI(2)));
            legend('Location','northwestoutside');
            axis square;
        end
    end
else
    warning('variant_preds not found — ablation correlation plots (LME) skipped.');
end
