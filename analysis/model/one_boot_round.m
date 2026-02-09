function [r_b, rmse_b] = one_boot_round(EmpAll, VarPreds, groupLabels, setSizes, variantNames, gu, grpIdx, RMSE_AGG, R_AGG)
% One bootstrap replicate: resample subjects within each group, make an "ALL" pool, score.

    groupLabels = string(groupLabels);
    nVar = numel(VarPreds);

    % Build resampled indices per group, then pool for "ALL"
    idx_all = [];
    resamp_per_group = cell(numel(gu),1);
    for gi = 1:numel(gu)
        idx = grpIdx{gi};
        if numel(idx) >= 2
            resamp = idx(randsample(numel(idx), numel(idx), true)); % with replacement
        else
            resamp = idx; % not enough to resample; just take what exists
        end
        resamp_per_group{gi} = resamp;
        idx_all = [idx_all; resamp(:)]; %#ok<AGROW>
    end

    % Group names include "ALL" first
    groupNames = ["ALL"; gu];
    nGrp = numel(groupNames);

    r_b    = nan(nVar, nGrp);
    rmse_b = nan(nVar, nGrp);

    for vi = 1:nVar
        for gj = 1:nGrp
            if gj == 1
                idx_use = idx_all;     % ALL
            else
                idx_use = resamp_per_group{gj-1};
            end
            if isempty(idx_use) || numel(idx_use) < 2
                continue;
            end

            Emp  = EmpAll(idx_use, :);
            Pred = VarPreds{vi}(idx_use, :);
            labs = groupLabels(idx_use);

            % Score with your existing function (per-K + aggregated)
            fitb = score_model_fit(Emp, Pred, labs, setSizes, variantNames(vi), RMSE_AGG, R_AGG);
            Tover = fitb.overall;

            % Extract rows for this gj ("ALL" or specific group)
            tg = groupNames(gj);
            row = strcmpi(Tover.Group, tg);
            if any(row)
                rmse_b(vi, gj) = Tover.RMSE_overall(row);
                % use the aggregated Fisher-z Pearson across K:
                if ismember('Pearson_r_overall_agg', Tover.Properties.VariableNames)
                    r_b(vi, gj) = Tover.Pearson_r_overall_agg(row);
                else
                    % fallback to pooled if agg not available
                    r_b(vi, gj) = Tover.Pearson_r_overall(row);
                end
            end
        end
    end
end
