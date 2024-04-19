function pVal = stat_compare1D(dataMat, varargin)

% receives a data matrix and compares groups (columns) based on
%
% INPUT:
%   dataMat         n x m numeric where n is number of subjects and
%   flgParam        logical. perform a parametric (true) or non-paramteric test
%   flgRep          logical. repeated measures (true) or independent
%   axh             axis handle. if given will add p-value to center
%
% OUTPUT
%   pVal            p-value of statistical test
%
% 23 mar 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'flgParam', false, @islogical);
addParameter(p, 'flgRep', false, @islogical);
addParameter(p, 'axh', [])

parse(p, varargin{:})
flgParam       = p.Results.flgParam;
flgRep         = p.Results.flgRep;
axh            = p.Results.axh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ngrps = size(dataMat, 2);

if ngrps == 2

    if flgRep
        if ~flgParam
            pVal = signrank(dataMat(:, 1), dataMat(:, 2));      % Wilcoxon signed-rank test
        else
            [~, pVal] = ttest(dataMat(:, 1), dataMat(:, 2));    % Paired t-test
        end
    else
        if ~flgParam
            pVal = ranksum(dataMat(:, 1), dataMat(:, 2));       % Mann-Whitney U test
        else
            [~, pVal] = ttest2(dataMat(:, 1), dataMat(:, 2));   % Two-sample t-test
        end
    end

elseif ngrps > 2

    if flgRep
        if ~flgParam
            pVal = friedman(dataMat, 1, 'off');                 % Friedman's test
        else
            error('Repeated measures ANOVA in MATLAB requires a table setup and is not covered in this function.');
        end
    else

        % Reshape data
        dataVector = reshape(dataMat, [], 1);
        groupVector = repelem((1:ngrps), size(dataMat, 1))';

        if flgParam
            pVal = kruskalwallis(dataVector, groupVector, 'off');  % Kruskal-Wallis
        else
            pVal = anova1(dataVector, groupVector, 'off');    % ANOVA
        end
    end
else
    error('Data must contain at least 2 groups for comparison.');
end

% Annotate plot with p-value if axh is provided
if ~isempty(axh) && isgraphics(axh, 'axes')
    axes(axh);
    ylims = get(axh, 'YLim');
    xlims = get(axh, 'XLim');
    text(mean(xlims), ylims(2), sprintf('p=%.3f', pVal),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Parent', axh);
end

end