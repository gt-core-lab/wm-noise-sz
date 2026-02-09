function [sigma_map_all, tune_map_all, decay_map_all] = map_subject_params( ...
    valid, useZ, ...
    map_sigma_int_base, map_sigma_int_slope, map_ext_filt_base, map_ext_filt_slope, ...
    raw_sigma_scale, raw_phi_scale, ...
    z_sigma, z_phi, z_dom, PTMtbl)

    nSub = height(PTMtbl);

    sigma_map_all = nan(nSub,1);   % post-decoding additive angular SD
    tune_map_all  = nan(nSub,1);   % external filtering for tuning (from φ / z_phi)
    decay_map_all = nan(nSub,1);   % external filtering for decay (from z_dom)

    if useZ
        % additive (post-decoding)
        sigma_map_all(valid) = max(1e-6, map_sigma_int_base + map_sigma_int_slope * z_sigma(valid));
        % tuning (φ only)
        tune_map_all(valid)  = max(0.0 , map_ext_filt_base  + map_ext_filt_slope  * z_phi(valid));
        % decay (dominant of z’s)
        decay_map_all(valid) = max(0.0 , map_ext_filt_base  + map_ext_filt_slope  * z_dom(valid));
    else
        % raw fallback: scale σ_add and φ; decay uses max(raw φ, raw σ_add)
        sigma_map_all(valid) = max(1e-6, raw_sigma_scale * PTMtbl.sigma_add(valid));
        tune_map_all(valid)  = max(0.0 , raw_phi_scale   * PTMtbl.phi(valid));
        decay_map_all(valid) = max(0.0 , raw_phi_scale   * max(PTMtbl.phi(valid), PTMtbl.sigma_add(valid)));
    end
end