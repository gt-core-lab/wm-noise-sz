function out = score_model_fit(EmpMat, PredMat, groupLabels, setSizes, variantName, overall_method, overall_r_method)
% Compute GOF between EmpMat and PredMat (nSub x nSet), per-K and overall.
% overall_method: 'sqmean' (default), 'mean', or 'median'
% overall_r_method: 'fisher' (default), 'mean', or 'median' -> for overall Pearson from per-K
%
% Returns:
%   out.perK     : table(Variant,Group,K, RMSE, MAE, r_Pearson, r_Spearman, r_Kendall, Bias, Slope_Deming)
%   out.overall  : table(Variant,Group, RMSE_overall, method)
%   out.tidy     : long table combining perK and overall rows

    if nargin < 6 || isempty(overall_method)
        overall_method = 'sqmean';  % √mean(RMSE^2) across K (weighted by #obs if K sizes differ)
    end

    if nargin < 7 || isempty(overall_r_method)
        overall_r_method = 'fisher';  % aggregate per-K Pearson via Fisher z (weighted)
    end


    % Vectorize
    [nSub, nSet] = size(EmpMat);
    x_all = EmpMat(:);
    y_all = PredMat(:);
    g_all = repmat(string(groupLabels(:)), nSet, 1);
    K_all = repmat(setSizes(:).', nSub, 1);   %  [nSub x nSet]
    K_all = K_all(:);

    % Keep finite numeric pairs
    keep = isfinite(x_all) & isfinite(y_all);
    x_all = x_all(keep);
    y_all = y_all(keep);
    g_all = g_all(keep);
    K_all = K_all(keep);

    % ---------- per-K metrics (overall + by group) ----------
    perK = table();
    Klist = unique(K_all(:),'stable');
    Glist = ["ALL"; unique(g_all,'stable')];

    for Ki = 1:numel(Klist)
        kmask = (K_all == Klist(Ki));

        % ALL
        perK = [perK; score_one_group(x_all(kmask), y_all(kmask), "ALL", variantName, Klist(Ki))]; %#ok<AGROW>

        % by group
        for gi = 2:numel(Glist)
            gmask = kmask & (g_all == Glist(gi));
            perK = [perK; score_one_group(x_all(gmask), y_all(gmask), Glist(gi), variantName, Klist(Ki))]; %#ok<AGROW>
        end
    end

    % ---------- overall aggregation from per-K RMSEs + overall Pearson r ----------
    overall = table();
    for gi = 1:numel(Glist)
        gname = Glist(gi);

        % rows for this group in the per-K table (for RMSE aggregation and per-K r's)
        rows = perK(perK.Group == gname, :);
        w_rmse = rows.Nobs;               % weights = #obs per K
        rmseK  = rows.RMSE;

        % --- aggregated RMSE across K
        switch lower(overall_method)
            case 'sqmean'   % √mean(RMSE^2) with weights
                if all(isfinite(w_rmse)) && any(w_rmse>0)
                    RMSE_overall = sqrt( sum(w_rmse .* (rmseK.^2), 'omitnan') / sum(w_rmse,'omitnan') );
                else
                    RMSE_overall = sqrt( mean(rmseK.^2,'omitnan') );
                end
            case 'mean'
                RMSE_overall = mean(rmseK,'omitnan');
            case 'median'
                RMSE_overall = median(rmseK,'omitnan');
            otherwise
                RMSE_overall = sqrt( mean(rmseK.^2,'omitnan') );
        end

        % --- overall Pearson r (pooled across K, using all pairs) for reference
        if gname == "ALL"
            mask = true(size(x_all));
        else
            mask = (g_all == gname);
        end
        Pearson_r_overall_pooled = safe_corr(x_all(mask), y_all(mask), 'Pearson');

        % --- overall Pearson r aggregated from per-K r's (preferred)
        rK = rows.Pearson_r;                   % per-K Pearson r
        wK = max(rows.Nobs - 3, 1);            % Fisher z weights (N-3), clipped at 1
        rK = rK(isfinite(rK)); wK = wK(isfinite(rK));  %#ok<NASGU> (keep rK,wK synced)

        switch lower(overall_r_method)
            case 'fisher'
                % Fisher z transform, weighted by (N-3)
                if ~isempty(rK)
                    % re-sync rK and wK after removing non-finite rK
                    rows_valid = isfinite(rows.Pearson_r);
                    rK = rows.Pearson_r(rows_valid);
                    wK = max(rows.Nobs(rows_valid) - 3, 1);
                    z  = atanh(max(min(rK, 0.999999), -0.999999));
                    zbar = sum(wK .* z, 'omitnan') / sum(wK, 'omitnan');
                    Pearson_r_overall_agg = tanh(zbar);
                else
                    Pearson_r_overall_agg = NaN;
                end
            case 'mean'
                Pearson_r_overall_agg = mean(rK,'omitnan');
            case 'median'
                Pearson_r_overall_agg = median(rK,'omitnan');
            otherwise
                Pearson_r_overall_agg = NaN;
        end

        % assemble overall row
        overall = [overall; table( ...
            string(variantName), string(gname), string(overall_method), ...
            RMSE_overall, ...
            Pearson_r_overall_pooled, ...
            Pearson_r_overall_agg, ...
            string(overall_r_method), ...
            'VariableNames', {'Variant','Group','OverallMethod','RMSE_overall', ...
                              'Pearson_r_overall_pooled','Pearson_r_overall_agg','OverallRMethod'})]; %#ok<AGROW>
    end


    % ---------- bundle outputs ----------
    out.perK    = perK;
    out.overall = overall;

    % ---------- long tidy table (per-K rows + OVERALL rows with same schema) ----------
    % Desired unified column order
    varsOrder = {'Variant','Group','K','RMSE','MAE','Pearson_r','Spearman_rho','Kendall_tau','Bias','Slope_Deming','Nobs'};

    % per-K rows -> strings for K like "K=1"
    P = perK(:, {'Variant','Group','K','RMSE','MAE','Pearson_r','Spearman_rho','Kendall_tau','Bias','Slope_Deming','Nobs'});
    P.K = "K=" + string(P.K);   % convert numeric K to string label

    % overall rows -> same columns, with K="OVERALL" and NaNs for per-K-only metrics
    K_overall = repmat("OVERALL", height(overall), 1);
    RMSE_overall = overall.RMSE_overall;
    % Nobs for overall = sum of Nobs across Ks for that group
    N_overall = nan(height(overall),1);
    for i = 1:height(overall)
        gname = overall.Group(i);
        N_overall(i) = nansum(perK.Nobs(perK.Group==gname));
    end
    tidy_overall = table(overall.Variant, overall.Group, K_overall, RMSE_overall, ...
                         nan(height(overall),1), nan(height(overall),1), nan(height(overall),1), ...
                         nan(height(overall),1), nan(height(overall),1), nan(height(overall),1), N_overall, ...
        'VariableNames', varsOrder);
    % Ensure per-K table has the same order
    P = P(:, varsOrder);

    % Final tidy table
    out.tidy = [P; tidy_overall];


end