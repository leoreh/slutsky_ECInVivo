function plot_ccg(cc, ccbins, varargin)

% plots cross correlation histrogram. can mark significant bins for mono
% synaptic interactions (see monoSyn_wrapper.m)
%
% INPUT
%   cc          numeric n x 1 vector of cross correlation 
%   ccbins      numeric n x 1 of time lags 
%   pred        numeric n x 1 of the predicator for a deconvoluted ccg
%   sigbins1    numeric n x 1 of significant excitatory bins. these are
%               indices to cc and ccbins
%   sigbins2    same as sigbins 2 for significant inhibitory bins
%   clr         char. color of ccg {'k'}
%
% SEE ALSO
%   plot_monoSyn.m
% 
% 18 jun 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'pred', [], @isnumeric);
addOptional(p, 'sigbins1', [], @isnumeric);
addOptional(p, 'sigbins2', [], @isnumeric);
addOptional(p, 'clr', 'k', @ischar);

parse(p,varargin{:})
pred        = p.Results.pred;
sigbins1    = p.Results.sigbins1;
sigbins2    = p.Results.sigbins2;
clr         = p.Results.clr;

% color of significant e / i bins
clr_sig = ['br'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bh = bar(ccbins, cc, 'BarWidth', 1);
bh.FaceColor = clr;
bh.FaceAlpha = 0.4;
bh.EdgeColor = 'none';
hold on

% add significant bins
if ~isempty(sigbins1)
    bh = bar(ccbins(sigbins1), cc(sigbins1), 'BarWidth', 1);
    bh.FaceColor = clr_sig(1);
    bh.FaceAlpha = 1;
    bh.EdgeColor = 'none';
end
if ~isempty(sigbins2)
    bh = bar(ccbins(sigbins2), cc(sigbins2), 'BarWidth', 1);
    bh.FaceColor = clr_sig(2);
    bh.FaceAlpha = 1;
    bh.EdgeColor = 'none';
end
if ~isempty(pred)
    plot(ccbins, pred, '--k', 'LineWidth', 2)
end

plot([0, 0], ylim, '--k')
ylabel('Counts')
xlabel('Time [s]')
box off
axis tight

end