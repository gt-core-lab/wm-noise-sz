function J = calib_loss_map_clean(p, ...
    trialSubIDs, PTMtbl, valid, useZ, z_sigma, z_phi, z_dom, ...
    setSizes, neuronsPerPool, kappa_tune_base, baseline_rate, ...
    poisson_like_mult, gamma_setsize, diffusion_extGain, diffusion_setGain, ...
    delaySec, probeRandomItem, raw_sigma_scale, raw_phi_scale, ...
    target_invK, groupLabels, w_mean, w_slope, ...
    tune_filter_gain)

    % unpack
    map_sigma_int_base  = exp(p(1));
    map_sigma_int_slope = exp(p(2));
    map_ext_filt_base   = exp(p(3));
    map_ext_filt_slope  = exp(p(4));
    gain_base           = exp(p(5));
    diffusion_base      = exp(p(6));

    % subject params
    [sigma_map_all, tune_map_all, decay_map_all] = map_subject_params( ...
        valid, useZ, ...
        map_sigma_int_base, map_sigma_int_slope, map_ext_filt_base, map_ext_filt_slope, ...
        raw_sigma_scale, raw_phi_scale, ...
        z_sigma, z_phi, z_dom, PTMtbl);


    % quick sim
    prefAngles = linspace(0, 2*pi, neuronsPerPool+1); prefAngles(end) = [];
    pred_invK = simulate_pred_invK_core( ...
        setSizes, 250, neuronsPerPool, prefAngles, ...
        kappa_tune_base, baseline_rate, gain_base, gamma_setsize, ...
        diffusion_base, diffusion_extGain, diffusion_setGain, delaySec, ...
        poisson_like_mult, probeRandomItem, ...
        tune_filter_gain, ...
        valid, sigma_map_all, tune_map_all, decay_map_all, ...
        true, true, true);

    % (a) mean matching per K
    sim_invK = mean(pred_invK, 1, 'omitnan');
    E_mean = sum((sim_invK - target_invK).^2);

    % (b) slope toward 1 within each group & K
    pred_mat = pred_invK;          % [nSub x nSet]
    nSet = numel(setSizes);
    % build empirical matrix from PTMtbl
    kcols = arrayfun(@(K) sprintf('k%dvar', K), setSizes(:)', 'uni', false);
    empirical_mat = nan(numel(trialSubIDs), nSet);
    for si = 1:nSet
        empirical_mat(:,si) = PTMtbl.(kcols{si});
    end

    groups = unique(string(groupLabels), 'stable');
    E_slope = 0; W_slope = 0;
    for gi = 1:numel(groups)
        gmask = strcmp(string(groupLabels), groups(gi));
        for si = 1:nSet
            x = empirical_mat(gmask, si); y = pred_mat(gmask, si);
            ok = isfinite(x) & isfinite(y);
            if sum(ok) < 3, continue; end
            xc = x(ok) - mean(x(ok));
            yc = y(ok) - mean(y(ok));
            denom = max(xc' * xc, eps);
            slope = (xc' * yc) / denom;
            E_slope = E_slope + (slope - 1)^2;
            W_slope = W_slope + 1;
        end
    end
    if W_slope > 0, E_slope = E_slope / W_slope; end

    % total objective + mild regularization
    J = w_mean * E_mean + w_slope * E_slope;
    J = J + 1e-4 * sum((p - log([0.15,0.10,1.5,1.0,40,0.015])).^2);
end