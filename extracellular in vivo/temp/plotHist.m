function plotHist(histData, clr)
if isempty(clr)
    clr = 'k';
end
histBins = 200;
h = histogram(histData, histBins, 'Normalization', 'probability');
h.FaceColor = clr;
h.EdgeColor = 'none';
h.FaceAlpha = 0.2;
ylabel('Probability')
end