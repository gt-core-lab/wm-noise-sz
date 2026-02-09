function ceiling = noise_ceiling_split_half(Trials, setSizes, nBoot, seed)
% Returns a struct with:
%   ceiling.rel_mean, rel_ci_lo, rel_ci_hi (reliability of full data)
%   ceiling.rmax_mean, rmax_ci_lo, rmax_ci_hi (noise ceiling on correlation)
% Reliability computed across subjects, aggregated across set sizes using Fisher z.

if nargin < 3, nBoot = 1000; end
if nargin < 4, seed = 2025; end
rng(seed);

subIDs = unique(string(Trials.SubID), 'stable');
nSub   = numel(subIDs);

% Pre-allocate bootstrap distribution of reliability (overall agg across K)
rel_dist = nan(nBoot,1);

for b = 1:nBoot
    invK_A = nan(nSub, numel(setSizes));
    invK_B = nan(nSub, numel(setSizes));

    for s = 1:nSub
        sid = subIDs(s);

        for ki = 1:numel(setSizes)
            K = setSizes(ki);

            err = Trials.ErrorRad(string(Trials.SubID)==sid & Trials.N==K);
            err = err(isfinite(err));
            if numel(err) < 10, continue; end

            % random split
            idx = randperm(numel(err));
            half = floor(numel(err)/2);
            eA = err(idx(1:half));
            eB = err(idx(half+1:2*half));

            invK_A(s,ki) = invK_from_errors(eA);
            invK_B(s,ki) = invK_from_errors(eB);
        end
    end

    % correlate across subjects, per K, then Fisher-aggregate across K
    rK = nan(1,numel(setSizes));
    for ki = 1:numel(setSizes)
        x = invK_A(:,ki); y = invK_B(:,ki);
        ok = isfinite(x) & isfinite(y);
        if sum(ok) >= 5
            rK(ki) = corr(x(ok), y(ok), 'Type', 'Pearson');
        end
    end

    % Fisher agg across K (like your RAgg='fisher')
    ok = isfinite(rK) & abs(rK) < 0.9999;
    if any(ok)
        z = atanh(rK(ok));
        rel_dist(b) = tanh(mean(z));
    end
end

% Split-half -> full reliability via Spearman–Brown
rel_full_dist = (2*rel_dist) ./ (1 + rel_dist);

% Noise ceiling on correlation with a noiseless predictor
rmax_dist = sqrt(max(rel_full_dist, 0));

% Summaries
ceiling.rel_mean = mean(rel_full_dist, 'omitnan');
ceiling.rel_ci_lo = prctile(rel_full_dist, 16);
ceiling.rel_ci_hi = prctile(rel_full_dist, 84);

ceiling.rmax_mean = mean(rmax_dist, 'omitnan');
ceiling.rmax_ci_lo = prctile(rmax_dist, 16);
ceiling.rmax_ci_hi = prctile(rmax_dist, 84);

ceiling.rel_dist = rel_full_dist;
ceiling.rmax_dist = rmax_dist;
end

function invK = invK_from_errors(err)
% err in radians, assumed in [-pi, pi] (or wrapped)
R = abs(mean(exp(1i*err)));
kappa = vm_kappa_from_R(R);   % you already have this
invK = 1 ./ kappa;
end
