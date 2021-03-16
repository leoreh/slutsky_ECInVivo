
% removes spikes that are rvds based on fet / mDist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
basepath = 'K:\Data\lh86\lh86_210228_190000\1st_6hr\lh86_210228_190000';
cd(basepath)
[~, basename] = fileparts(basepath);

session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

psamp = [];
grps = [];

% files params
sniplength = ceil(1.6 * 10^-3 * fs);
win = [-(floor(sniplength / 2) - 1) floor(sniplength / 2)];   
precision = 'int16'; % for dat file. size of one data point in bytes
nbytes = class2bytes(precision); 

if isempty(psamp)
    psamp = round(sniplength / 2);
end
if isempty(grps)
    grps = 1 : length(spkgrp);
end
ngrps = length(grps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups and clus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1 : ngrps
    
    grp = grps(j);
    grpchans = spkgrp{j};
    
    % ---------------------------------------------------------------------
    % load data for neurosuite files
    
    % clu
    cluname = fullfile([basename '.clu.' num2str(grp)]);
    fid = fopen(cluname, 'r');
    nclu = fscanf(fid, '%d\n', 1);
    clu = fscanf(fid, '%d\n');
    rc = fclose(fid);
    if rc ~= 0 || isempty(clu)
        warning(['failed to read clu ' num2str(j)])
    end
    nspks(j) = length(clu);
    uclu = unique(clu);
    
    % res 
    resname = fullfile([basename '.res.' num2str(grp)]);
    fid = fopen(resname, 'r');
    res = fscanf(fid, '%d\n');
    rc = fclose(fid);
    if rc ~= 0 || isempty(res)
        warning(['failed to read res ' num2str(j)])
    end
    
    % spk
    spkname = fullfile([basename '.spk.' num2str(grp)]);
    fid = fopen(spkname, 'r');
    spk = fread(fid, 'int16');
    spk = reshape(spk, length(grpchans), sniplength, nspks(j)); 
    rc = fclose(fid);
    if rc ~= 0 || isempty(spk)
       warning(['failed to read spk ' num2str(j)])
    end
     
    % fet
    fetname = fullfile([basename '.fet.' num2str(grp)]);  
    fid = fopen(fetname, 'r');
    formatstring = '%d';
    for ii = 2 : nFeatures
        formatstring = [formatstring, '\t%d'];
    end
    formatstring = [formatstring, '\n'];  
    nfet = fscanf(fid, '%d\n', 1);
    fet = fscanf(fid, formatstring);
    fet = reshape(fet, nspks(j), nfet);
    npca = 3;
    nfet = npca * length(spkgrp{j});
    fet = fet(:, 1 : nfet);
    rc = fclose(fid);
    if rc ~= 0 || isempty(fet)
        warning(['failed to read fet ' num2str(j)])
    end

    for jj = 1 : nclu
        if uclu(jj) == 0
            continue
        end
        cluidx = find(clu == uclu(jj));
        fetclu = fet(cluidx, :);
        
        % mdist
        mDist = mahal(fet, fet(cluidx, :));
        mCluster = mDist(cluidx); % mahal dist of spikes in cluster        
        
        % rvd
        ref = ceil(0.002 * fs);
        rvd = find(diff([0; res(cluidx)]) < ref);
        
        % bins for mDist
        bins = linspace(0, max(mCluster), 500);
        binIntervals = [0 bins(1 : end - 1); bins]';

        % test removal
        nrmv = histcounts(mCluster(rvd), 'BinEdges', [2.3 60]);
        mu = (length(rvd) - nrmv) / (length(cluidx) - nrmv) * 100;
        
        % graphics
        figure
        [nsub] = numSubplots(nfet + 1);
        xLim = [0 25000];
        for jjj = 1 : nfet
            subplot(nsub(1), nsub(2), jjj)        
            histogram(fetclu(:, jjj), 'BinEdges', linspace(xLim(1), xLim(2), 500),...
                'FaceAlpha', 0.4, 'LineStyle', 'none')
            hold on
            sh = scatter(fetclu(rvd, jjj), mean(get(gca, 'YLim')) *...
                ones(length(rvd), 1), 'k');
            set(sh, 'Marker', 'x')
            yyaxis right
            histogram(fetclu(rvd, jjj), 'BinEdges', linspace(xLim(1), xLim(2), 500),...
                'FaceAlpha', 0.4, 'LineStyle', 'none')
            xlim(xLim)
            title(['fet #', num2str(jjj)])
        end
        subplot(nsub(1), nsub(2), jjj + 1)
        histogram(mCluster, 'BinEdges', bins, 'FaceAlpha', 0.4, 'LineStyle', 'none')
        hold on
        sh = scatter(mCluster(rvd), mean(get(gca, 'YLim')) *...
            ones(length(rvd), 1), 'k');
        yyaxis right
        histogram(mCluster(rvd), 'BinEdges', linspace(xLim(1), xLim(2), 500),...
            'FaceAlpha', 0.4, 'LineStyle', 'none')
        title('mDist')
        suptitle(['clu #', num2str(jj)])
        
    end
      
    % ---------------------------------------------------------------------
    % save files
    
 
end
     
% EOF

% -------------------------------------------------------------------------
% load and inspect rvds in parallel to mancur
j = 2;      % spkgrp
% clu
cluname = fullfile([basename '.clu.' num2str(j)]);
fid = fopen(cluname, 'r');
nclu = fscanf(fid, '%d\n', 1);
clu = fscanf(fid, '%d\n');
rc = fclose(fid);
if rc ~= 0 || isempty(clu)
    warning(['failed to read clu ' num2str(j)])
end
nspks(j) = length(clu);
uclu = unique(clu);

% res
resname = fullfile([basename '.res.' num2str(j)]);
fid = fopen(resname, 'r');
res = fscanf(fid, '%d\n');
rc = fclose(fid);
if rc ~= 0 || isempty(res)
    warning(['failed to read res ' num2str(j)])
end

ref = ceil(0.002 * fs);
for jj = 1 : nclu
        if uclu(jj) == 0 || uclu(jj) == 1
            continue
        end
        cluidx = find(clu == uclu(jj));
    
        % rvd
        rvd = find(diff([0; res(cluidx)]) < ref);
        mu(jj) = length(rvd) / length(cluidx) * 100;
 end