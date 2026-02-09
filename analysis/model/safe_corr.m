function r = safe_corr(x, y, type)
    try
        r = corr(x, y, 'type', type, 'rows','pairwise');
    catch
        m = isfinite(x) & isfinite(y);
        if sum(m) < 3, r = NaN; return; end
        switch lower(type)
            case 'pearson',  c = corrcoef(x(m), y(m)); r = c(1,2);
            case 'spearman', r = corr(x(m), y(m), 'type','Spearman');
            case 'kendall',  r = corr(x(m), y(m), 'type','Kendall');
            otherwise, r = NaN;
        end
    end
end