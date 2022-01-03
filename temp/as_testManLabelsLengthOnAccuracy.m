

labelsMan = labels;
labelsNet = ss.labels_net;
netScores = ss.netScores;

start = 86 * 60;
jumps = [30, 60, Inf];

for ilabels = 1 : length(jumps)
    
    if jumps(ilabels) == Inf
        labelsIdx = 1 : length(labelsMan);
    else
        labelsIdx = start : start + jumps(ilabels) * 60;
    end
    labels1 = labelsMan(labelsIdx);
    labels2 = labelsNet(labelsIdx);
    
    [prc(ilabels, :), rec(ilabels, :)] = as_cm(labels1, labels2, netScores,...
        'graphics', false);
    
end

% get params from configuration file
[cfg_colors, cfg_names, ~] = as_loadConfig([]);
cfg_colors = cfg_colors(:);
nstates = length(cfg_colors) - 1;
setMatlabGraphics(false)

fh = figure;
subplot(1, 2, 1)
ph = plot(prc, 'LineWidth', 2);
set(ph, {'color'}, cfg_colors(1 : nstates))
legend(cfg_names(1 : nstates), 'Location', 'southeast');
ylabel('Precision')
xticks([1 : length(jumps)])
xticklabels(split(num2str(jumps)))
xlabel('Manual Labels Length [min]')

subplot(1, 2, 2)
ph = plot(rec, 'LineWidth', 2);
set(ph, {'color'}, cfg_colors(1 : nstates))
ylabel('Recall')
xticks([1 : length(jumps)])
xticklabels(split(num2str(jumps)))
xlabel('Manual Labels Length [min]')


