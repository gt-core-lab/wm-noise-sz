function boot = bootstrap_model_fit(empirical_invK, variant_preds, groupLabels, setSizes, variantNames, varargin)
% BOOTSTRAP_MODEL_FIT
% Bootstrap overall GOF metrics (RMSE and Fisher-aggregated Pearson r) per
% (variant × group), resampling subjects WITHIN each group.
%
% Inputs
%   empirical_invK : [nSub x nSet] empirical matrix
%   variant_preds  : {nVar x 1} cell, each [nSub x nSet] predictions for a variant
%   groupLabels    : [nSub x 1] string/cellstr group label per subject (e.g., "NT","SZ")
%   setSizes       : [1 x nSet] list of K (e.g., [1 2 4])
%   variantNames   : [nVar x 1] string array of variant names
%
% Name-value (optional)
%   'NBoot'          : default 1000
%   'RMSEAgg'        : 'sqmean'|'mean'|'median'      (default 'sqmean')
%   'RAgg'           : 'fisher'|'mean'|'median'      (default 'fisher')
%   'Seed'           : numeric RNG seed (default 42)
%   'UseParallel'    : true/false (default false; requires PCT)
%
% Output struct 'boot' fields:
%   .groups              : string ["ALL", unique(groups)...]
%   .variants            : string variant names
%   .r.dist              : [nVar x nGrp x nBoot] Fisher-aggregated Pearson r
%   .rmse.dist           : [nVar x nGrp x nBoot] overall RMSE
%   .r.mean / .rmse.mean : [nVar x nGrp]
%   .r.ci_lo/ci_hi       : [nVar x nGrp] 2.5/97.5 percentiles
%   .rmse.ci_lo/ci_hi    : same
%
% Requires helper:
%   score_model_fit(Emp, Pred, groupLabels, setSizes, variantName, RMSEAgg, RAgg)

    p = inputParser;
    addParameter(p, 'NBoot', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p, 'RMSEAgg', 'sqmean', @(s)ischar(s)||isstring(s));
    addParameter(p, 'RAgg', 'fisher', @(s)ischar(s)||isstring(s));
    addParameter(p, 'Seed', 42, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'UseParallel', false, @(x)islogical(x)&&isscalar(x));
    parse(p, varargin{:});
    NBOOT       = p.Results.NBoot;
    RMSE_AGG    = string(p.Results.RMSEAgg);
    R_AGG       = string(p.Results.RAgg);
    SEED        = p.Results.Seed;
    USE_PAR     = p.Results.UseParallel;

    rng(SEED);
    groupLabels = string(groupLabels(:));
    gu = unique(groupLabels,'stable');
    nSub = size(empirical_invK,1);
    nVar = numel(variant_preds);

    % Precompute per-group subject indices
    grpIdx = cell(numel(gu),1);
    for gi = 1:numel(gu)
        grpIdx{gi} = find(groupLabels == gu(gi));
    end

    % Report an "ALL" group composed by pooling resampled subjects across groups
    groupNames = ["ALL"; gu];
    nGrp = numel(groupNames);

    % Allocate distributions
    r_dist    = nan(nVar, nGrp, NBOOT);
    rmse_dist = nan(nVar, nGrp, NBOOT);

    % Optional parallelization
    if USE_PAR && license('test','Distrib_Computing_Toolbox')
        parforArg = Inf; %#ok<NASGU> (just to silence lint)
    end

    % Bootstrap loop
    if USE_PAR && license('test','Distrib_Computing_Toolbox')
        parfor b = 1:NBOOT
            [r_b, rmse_b] = one_boot_round(empirical_invK, variant_preds, groupLabels, setSizes, variantNames, gu, grpIdx, RMSE_AGG, R_AGG);
            r_dist(:, :, b)    = r_b;
            rmse_dist(:, :, b) = rmse_b;
        end
    else
        for b = 1:NBOOT
            [r_b, rmse_b] = one_boot_round(empirical_invK, variant_preds, groupLabels, setSizes, variantNames, gu, grpIdx, RMSE_AGG, R_AGG);
            r_dist(:, :, b)    = r_b;
            rmse_dist(:, :, b) = rmse_b;
        end
    end

    % Summaries (68% CI instead of 95%)
    ciBounds = [16 84];   % percentiles for central 68%
    prc = @(A,pct) prctile(A, pct, 3);  % along bootstrap dim
    boot.groups   = groupNames;
    boot.variants = variantNames(:);

    boot.r.dist   = r_dist;
    boot.r.mean   = mean(r_dist, 3, 'omitnan');
    boot.r.ci_lo  = prc(r_dist, ciBounds(1));
    boot.r.ci_hi  = prc(r_dist, ciBounds(2));

    boot.rmse.dist  = rmse_dist;
    boot.rmse.mean  = mean(rmse_dist, 3, 'omitnan');
    boot.rmse.ci_lo = prc(rmse_dist, ciBounds(1));
    boot.rmse.ci_hi = prc(rmse_dist, ciBounds(2));

    % Also compute 95% CI (just in case)
    ciBounds95 = [2.5 97.5];

    boot.r.ci95_lo = prc(r_dist, ciBounds95(1));
    boot.r.ci95_hi = prc(r_dist, ciBounds95(2));

    boot.rmse.ci95_lo = prc(rmse_dist, ciBounds95(1));
    boot.rmse.ci95_hi = prc(rmse_dist, ciBounds95(2));


end