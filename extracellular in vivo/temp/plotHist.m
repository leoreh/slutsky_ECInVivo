function plotHist(histData, clr, histBins)
if isempty(clr)
    clr = 'k';
end
if isempty(histBins)
    histBins = 200;
end
h = histogram(histData, histBins, 'Normalization', 'probability');
h.FaceColor = clr;
h.EdgeColor = 'none';
h.FaceAlpha = 0.2;
ylabel('Probability')
end