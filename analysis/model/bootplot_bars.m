function bootplot_bars(metricMean, ciLo, ciHi, variantNames, groupNames, varargin)
% BOOTPLOT_BARS
% Draw grouped panels of bars with error bars (mean +/- CI).
% metricMean, ciLo, ciHi: [nVar x nGrp] matrices
% variantNames: [nVar x 1] string
% groupNames  : [nGrp x 1] string
% Name-value: 'YLabel', 'Title', 'YLim', 'BarColor', 'Rotate', 'Fmt'

    p = inputParser;
    addParameter(p, 'YLabel', '', @(s)ischar(s)||isstring(s));
    addParameter(p, 'Title', '', @(s)ischar(s)||isstring(s));
    addParameter(p, 'YLim', [], @(x)isnumeric(x)&&numel(x)==2 || isempty(x));
    addParameter(p, 'BarColor', [0.5 0.75 0.95], @(x)isnumeric(x)&&numel(x)==3);
    addParameter(p, 'Rotate', 25, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Fmt', '%.2f', @(s)ischar(s)||isstring(s));
    parse(p, varargin{:});
    YL = p.Results.YLabel; TTL = p.Results.Title; YLIM = p.Results.YLim;
    BC = p.Results.BarColor; ROT = p.Results.Rotate; FMT = char(p.Results.Fmt);

    [nVar, nGrp] = size(metricMean);
    figure('Color','w','Position',[160 140 280*max(3,nGrp) 420]);
    tiledlayout(1, nGrp, 'Padding','compact','TileSpacing','compact');

    for gi = 1:nGrp
        nexttile; 
        vals = metricMean(:, gi);
        errL = vals - ciLo(:, gi);
        errU = ciHi(:, gi) - vals;

        bh = bar(vals, 'FaceColor', BC); hold on; grid on; box on;
        eb = errorbar(1:nVar, vals, errL, errU, 'k.', 'LineWidth',1.1); %#ok<NASGU>
        set(gca,'XTick',1:nVar,'XTickLabel',variantNames,'XTickLabelRotation',ROT);

        if ~isempty(YLIM), ylim(YLIM); end
        ylabel(YL);
        title(string(groupNames(gi)));

        % numeric labels
        for i = 1:nVar
            if isfinite(vals(i))
                text(i, vals(i), sprintf(FMT, vals(i)), ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
            end
        end
    end
    sgtitle(string(TTL));
end