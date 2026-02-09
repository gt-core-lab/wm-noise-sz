function T = score_one_group(x, y, gname, vname, Kval)
% Core metrics for a single (group, K). Returns a single-row table.

    if nargin < 5 || isempty(Kval), Kval = NaN; end

    m = isfinite(x) & isfinite(y);
    x = x(m); y = y(m);
    N = numel(x);

    if N < 2
        rmse = NaN; mae = NaN; bias = NaN; rP = NaN; rS = NaN; rK = NaN; slope_dem = NaN;
    else
        rmse = sqrt(mean((y - x).^2,'omitnan'));
        mae  = mean(abs(y - x),'omitnan');
        bias = mean(y - x,'omitnan');
        rP   = safe_corr(x, y, 'Pearson');
        rS   = safe_corr(x, y, 'Spearman');
        rK   = safe_corr(x, y, 'Kendall');
        slope_dem = deming_slope(x, y, 1);
    end

    T = table(string(vname), string(gname), double(Kval), rmse, mae, rP, rS, rK, bias, slope_dem, N, ...
        'VariableNames', {'Variant','Group','K','RMSE','MAE','Pearson_r','Spearman_rho','Kendall_tau','Bias','Slope_Deming','Nobs'});
end