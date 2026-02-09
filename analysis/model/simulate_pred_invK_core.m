function pred_invK = simulate_pred_invK_core( ...
    setSizes, nTrialsPerSet, neuronsPerPool, prefAngles, ...
    kappa_tune_base, baseline_rate, gain_base, gamma_setsize, ...
    diffusion_base, diffusion_extGain, diffusion_setGain, delaySec, ...
    poisson_like_mult, probeRandomItem, ...
    tune_filter_gain, ...
    valid, sigma_map_all, tune_map_all, decay_map_all, ...
    useNoiseEffect, useTuneEffect, useDecayEffect)

    nSub = numel(valid);
    nSet = numel(setSizes);
    pred_invK = nan(nSub, nSet);

    % constants for OFF toggles (median across valid subjects)
    sigma_const = median(sigma_map_all(valid),'omitnan');
    tune_const  = median(tune_map_all(valid),'omitnan');
    decay_const = median(decay_map_all(valid),'omitnan');

    for s = 1:nSub
        if ~valid(s), continue; end

        % subject-specific mapped parameters
        sigma_int_sub   = sigma_map_all(s);   % post-decoding additive SD
        ext_tune_sub    = tune_map_all(s);    % φ-based mapping (tuning)
        ext_decay_sub   = decay_map_all(s);   % dominant-noise-based mapping (decay)

        % toggle ON/OFF per stream
        sigma_int  = useNoiseEffect * sigma_int_sub  + (~useNoiseEffect) * sigma_const;
        ext_tune   = useTuneEffect  * ext_tune_sub   + (~useTuneEffect)  * tune_const;
        ext_decay  = useDecayEffect * ext_decay_sub  + (~useDecayEffect) * decay_const;

        % tuning sharpness
        kappa_tune = kappa_tune_base;
        if useTuneEffect
            kappa_tune = kappa_tune_base / (1 + tune_filter_gain * ext_tune);
        end

        for si = 1:nSet
            K = setSizes(si);
            gain_item = gain_base / (K^gamma_setsize);

            % diffusion (with/without ext influence)
            if useDecayEffect
                D = diffusion_base * (1 + diffusion_extGain * ext_decay) ...
                                  * (1 + diffusion_setGain * (K - 1));
            else
                D = diffusion_base * (1 + diffusion_setGain * (K - 1));
            end
            sigmaTheta = sqrt(2 * D * delaySec);

            err = nan(nTrialsPerSet,1);
            for t = 1:nTrialsPerSet
                % ground truth and encoding
                theta_true = 2*pi*rand(K,1);
                theta_enc  = nan(K,1);

                for kk = 1:K
                    phi = prefAngles(:);

                    % mean rate
                    r_mu = baseline_rate + gain_item * exp(kappa_tune * cos(theta_true(kk) - phi));

                    % Poisson-like ONLY at neuron level
                    var_i = poisson_like_mult * max(r_mu, 0);
                    r = r_mu + sqrt(var_i) .* randn(neuronsPerPool,1);

                    % population vector decoding
                    th = circ_popvec(phi, r);
                    theta_enc(kk) = mod(th, 2*pi);
                end

                % post-decoding additive noise (angular)
                if useNoiseEffect && sigma_int > 0
                    theta_enc = add_circ_gauss(theta_enc, sigma_int);
                end

                % diffusion (always produce a memory angle vector)
                theta_mem = add_circ_gauss(theta_enc, sigmaTheta);

                % probe
                if probeRandomItem
                    j = randi(K);
                else
                    j = 1;
                end
                % use theta_mem consistently here
                err(t) = angle(exp(1i*(theta_mem(j) - theta_true(j)))); % [-pi, pi]
            end

            % R -> kappa -> invK
            R = abs(mean(exp(1i*err)));
            kappa = vm_kappa_from_R(R);
            pred_invK(s,si) = 1 / kappa;
        end
    end
end