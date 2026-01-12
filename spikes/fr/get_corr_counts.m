function [clean_counts, shuffle_counts, counts, bins, dbin, mcc, vcc, min_cc, max_cc, clean_cc] = get_corr_counts(Y,nbins,binrange,is_cc, threshold)

Y = Y ./ std(Y,[],2);
if is_cc

    Yz = sum(Y,2);
    if ~sum(Yz)
        return;
    end

    Ny = size(Y,1);

    % if threshold > 0
    Yshuffle = nan(size(Y));
    for n = 1:Ny
        Yshuffle(n,:) = circshift(Y(n,:),randi(size(Y,2)));
    end
    scc = cov(Yshuffle');
    scc = triu(scc) - diag(diag(scc));
    scc = scc(scc ~= 0);
    mscc = mean(scc(:));
    vscc = std(scc(:));
    % else
    % end

    cc = cov(Y');
    I = abs(cc) >= vscc*threshold;
    clean_cc = cc.*I;
    cc = triu(cc) - diag(diag(cc));
    cc = cc(cc ~= 0);
    full_cc = cc;
    cc = cc(abs(cc) >= vscc*threshold);

    mcc = mean(cc(:));
    vcc = std(cc(:));
    max_cc = prctile(cc(:),98);
    min_cc = prctile(cc(:),2);
    % num_cells(j) = Ny;
    % frac_nonzero(j) = Nnz/Ny;

else
    full_cc = Y;
    mcc = 0;
    vcc = 0;
end

bins = linspace(binrange(1),binrange(2),nbins+1);
dbin = bins(2) - bins(1);
counts = nan(length(bins)-1,1);
shuffle_counts = nan(length(bins)-1,1);
clean_counts = nan(length(bins)-1,1);

for jj = 1:length(bins)-1
    counts(jj) = sum((full_cc > bins(jj)) .* (full_cc  <= bins(jj+1)));
    shuffle_counts(jj) = sum((scc > bins(jj)) .* (scc  <= bins(jj+1)));
    clean_counts(jj) = sum((cc > bins(jj)) .* (cc  <= bins(jj+1)));
end

counts = counts ./ sum(counts);
shuffle_counts = shuffle_counts ./ sum(shuffle_counts);
clean_counts = clean_counts ./ sum(clean_counts);
