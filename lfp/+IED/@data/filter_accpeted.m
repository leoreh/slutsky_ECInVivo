function [accepeted_vec] = filter_accpeted(obj, accepted_filt, accepeted_vec)
    % Match a logical vector that was calc only on accepted IEDs,
    % to a full accpeted vector. Used to fix accpeted after you
    % Used for secound - stage, to remove IEDs that were accepted
    % earlier.
    % Note that order is important - values in accepted_filt should
    % correspond with find(accpeted_vec).
    %
    % INPUT:
    %   obj - ied object calling this function.
    %   accepted_filt - logical vector, numel(accepted_filt) == sum(accpeted_vec).
    %                   true for values that should remine accpeted.
    %   accpeted_vec  - logical vector, numel(accpeted_vec) == numel(ied.pos).
    %                   what positions of IED are accepted originaly.
    %                   if not given, use ied.accepted.
    %
    % OUTPUT:
    %   accpeted_vec  - same as accpeted_vec input, but only accpeted
    %                   values that were true in accepted_filt are
    %                   still true.
    %
    % EXAMPLE:
    %   %%% General: %%%
    %   accpeted_vec = [0 1 1 0];
    %   accepted_filt = [0 1];
    %   new_accepted = filter_accpeted(obj, accepted_filt, accpeted_vec)
    %
    %   %new_accepted = [0 0 1 0]
    %
    %   %%% Assign new accpeted %%%
    %   ied.accpeted = [0 1 1 0];
    %   ied.accpeted = filter_accpeted(obj, accepted_filt)

    if ~exist("accepeted_vec","var")
        accepeted_vec = obj.accepted;
    elseif numel(accepeted_vec) ~= numel(obj.pos) || ~islogical(accepeted_vec) || ~isvector(accepeted_vec)
        error('accpeted_vec need to be a logical vector, & match size with ied.pos')
    end

    if numel(accepted_filt) ~= sum(accepeted_vec)
        error("accepted_filt must have a value for each true val in accpeted_vec, total %d",sum(accepeted_vec))
    end

    true_pos = find(accepeted_vec);
    pos2rmv = true_pos(~accepted_filt);
    accepeted_vec(pos2rmv) = false;

end