function result = make_bouts_per_sample(weights)
%MAKE_EPOCHS_PER_SAMPLE Given a set of weights generate the number of
% bouts per sample for each weight.
%
% result = MAKE_EPOCHS_PER_SAMPLE(weights)
%
% Parameters
% ----------
% weights: array of size (n_1_simplices, 1)
%     The weights of how much we wish to sample each 1-simplex.
%
% Note that the total number of bouts does not impact this result.
% 
% Returns
% -------
% result: array of size (n_1_simplices, 1)
%     The number of bouts per sample, one for each 1-simplex.

%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

result = -1*ones(size(weights, 1), 1);
L = weights > 0;
result(L) = max(weights)*ones(sum(L), 1)./weights(L);